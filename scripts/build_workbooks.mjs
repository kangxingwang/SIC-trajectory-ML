import fs from "node:fs/promises";
import path from "node:path";
import { Workbook, SpreadsheetFile } from "@oai/artifact-tool";

const root = path.resolve(process.env.SIC_PROJECT_ROOT || path.join(import.meta.dirname, ".."));
const tableDir = path.join(root, "results", "tables");
const outDir = path.join(root, "results", "workbooks");
const previewDir = path.join(root, "results", "logs", "workbook_previews");
await fs.mkdir(outDir, { recursive: true }); await fs.mkdir(previewDir, { recursive: true });

async function csvValues(file) {
  const text = await fs.readFile(path.join(tableDir, file), "utf8");
  const temp = await Workbook.fromCSV(text, { sheetName: "Imported" });
  const rows=temp.worksheets.getItem("Imported").getUsedRange(true).values;
  const num=/^-?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?$/;
  return rows.map((r,i)=>i===0?r:r.map(v=>typeof v==="string"&&num.test(v.trim())?Number(v):v));
}
function colLetter(n) { let s=""; while(n){ n--; s=String.fromCharCode(65+n%26)+s; n=Math.floor(n/26); } return s; }
async function addSheet(wb, name, title, csvFile, note="Aggregate results only; no patient-level identifiers.") {
  const rows = await csvValues(csvFile); const cols = Math.max(...rows.map(r=>r.length)); const end = colLetter(cols);
  const sh = wb.worksheets.add(name); sh.showGridLines = false;
  sh.getRange(`A1:${end}1`).merge(); sh.getRange("A1").values=[[title]];
  sh.getRange(`A2:${end}2`).merge(); sh.getRange("A2").values=[[note]];
  sh.getRange(`A4:${end}${rows.length+3}`).values = rows;
  sh.getRange(`A1:${end}1`).format={fill:"#17365D",font:{bold:true,color:"#FFFFFF",size:15},horizontalAlignment:"center",verticalAlignment:"center"};
  sh.getRange(`A2:${end}2`).format={fill:"#D9EAF7",font:{italic:true,color:"#334155",size:9},wrapText:true};
  sh.getRange(`A4:${end}4`).format={fill:"#2B6CB7",font:{bold:true,color:"#FFFFFF"},horizontalAlignment:"center",verticalAlignment:"center",wrapText:true,
    borders:{bottom:{style:"medium",color:"#17365D"}}};
  if(rows.length>1) sh.getRange(`A5:${end}${rows.length+3}`).format={borders:{insideHorizontal:{style:"thin",color:"#E5E7EB"}},verticalAlignment:"center"};
  if(rows.length>1){for(let c=0;c<cols;c++){const vals=rows.slice(1).map(r=>r[c]).filter(v=>v!==null&&v!=="");if(vals.length&&vals.every(v=>typeof v==="number")){const h=String(rows[0][c]||"");sh.getRangeByIndexes(4,c,rows.length-1,1).format.numberFormat=/^N$|count|size|overlap/i.test(h)?"#,##0":"0.000";}}}
  sh.getRange(`A1:${end}${rows.length+3}`).format.autofitColumns(); sh.getRange(`A1:${end}${rows.length+3}`).format.autofitRows();
  for(let c=0;c<cols;c++){const rg=sh.getRangeByIndexes(0,c,rows.length+3,1);const want=Math.min(28,Math.max(10,String(rows[0][c]||"").length+3));if(rg.format.columnWidth<want)rg.format.columnWidth=want;if(rg.format.columnWidth>28)rg.format.columnWidth=28;}
  sh.getRange("A1").format.rowHeight=28; sh.getRange("A2").format.rowHeight=30; sh.freezePanes.freezeRows(4);
  return sh;
}
async function exportWorkbook(filename, specs) {
  const wb=Workbook.create();
  for(const s of specs) await addSheet(wb,s.name,s.title,s.csv,s.note);
  const inspect=await wb.inspect({kind:"sheet,table",maxChars:4000,tableMaxRows:5,tableMaxCols:10});
  await fs.writeFile(path.join(previewDir,`${filename}.inspect.ndjson`),inspect.ndjson,"utf8");
  const errors=await wb.inspect({kind:"match",searchTerm:"#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A",options:{useRegex:true,maxResults:100},summary:"formula error scan"});
  await fs.writeFile(path.join(previewDir,`${filename}.errors.ndjson`),errors.ndjson,"utf8");
  for(const s of specs){const blob=await wb.render({sheetName:s.name,autoCrop:"all",scale:1,format:"png"});await fs.writeFile(path.join(previewDir,`${filename}_${s.name}.png`),new Uint8Array(await blob.arrayBuffer()));}
  const out=await SpreadsheetFile.exportXlsx(wb);await out.save(path.join(outDir,filename));
}

await exportWorkbook("Feature_Engineering.xlsx",[
  {name:"Feature engineering",title:"Nested Boruta-LASSO feature engineering",csv:"feature_engineering_summary.csv",
   note:"Selection frequencies are estimated only within outer training folds; the three SIC components are prespecified and forced."}
]);
await exportWorkbook("Model_Performance.xlsx",[
  {name:"Model ranking",title:"Nested out-of-fold model ranking",csv:"nested_oof_model_ranking.csv"},
  {name:"Final parameters",title:"Final development tuning parameters",csv:"final_best_hyperparameters.csv"},
  {name:"Performance",title:"Temporal and external validation performance",csv:"performance_summary.csv"},
  {name:"Performance CI",title:"Performance estimates and bootstrap confidence intervals",csv:"performance_with_confidence_intervals.csv"},
  {name:"DeLong FDR",title:"Nested OOF AUROC comparisons with FDR correction",csv:"delong_comparisons_fdr.csv"},
  {name:"Thresholds",title:"Frozen development thresholds",csv:"development_frozen_thresholds.csv"}
]);
await exportWorkbook("SHAP_Interpretation.xlsx",[
  {name:"Global importance",title:"Cross-cohort global SHAP importance",csv:"shap_global_importance.csv"},
  {name:"Stability",title:"Cross-cohort SHAP rank stability",csv:"shap_cross_cohort_stability.csv"},
  {name:"Representative cases",title:"Representative local-explanation cases",csv:"shap_representative_cases.csv",
   note:"Case labels are descriptive and contain no patient identifiers. SHAP reflects predictive association, not causality."}
]);

console.log(`Created three workbooks in ${outDir}`);
