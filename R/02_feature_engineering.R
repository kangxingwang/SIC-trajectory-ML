## ============================================================================
## 1. Packages and configuration
## ============================================================================
if (!exists("SIC_CONFIG")) source(file.path(Sys.getenv("SIC_PROJECT_ROOT", unset = getwd()), "R", "00_config.R"), encoding = "UTF-8")
harmonized_file <- file.path(SIC_CONFIG$dirs$private, "harmonized_cohorts.rds")
if (!file.exists(harmonized_file)) stop("Run R/01_harmonize_data.R first.")

## ============================================================================
## 2. Import harmonized cohorts and define development split
## ============================================================================
obj <- readRDS(harmonized_file); cohorts <- obj$cohorts; FEATURES <- obj$features; FORCED_FEATURES <- obj$forced
m4 <- cohorts[["MIMIC-IV"]]
dev <- m4 %>% filter(Admission_Era != "2020-2022") %>% droplevels()
temporal <- m4 %>% filter(Admission_Era == "2020-2022") %>% droplevels()
externals <- cohorts[c("MIMIC-III","eICU","NWICU")]
if (nrow(dev) == 0 || nrow(temporal) == 0) stop("MIMIC-IV historical/temporal split is empty.")

## ============================================================================
## 3. Nested feature selection inside each outer training fold
## ============================================================================
impute_for_selection <- function(dat, features) {
  x <- dat[, c("Target",features), drop = FALSE]
  for (v in features) {
    if (is.numeric(x[[v]])) {
      x[[v]][!is.finite(x[[v]])] <- NA_real_; med <- median(x[[v]], na.rm = TRUE); if (!is.finite(med)) med <- 0
      x[[v]][is.na(x[[v]])] <- med
    } else {
      lev <- levels(x[[v]]); md <- mode_value(x[[v]]); z <- as.character(x[[v]]); z[is.na(z) | z == ""] <- md; x[[v]] <- factor(z, levels = lev)
    }
  }
  droplevels(x)
}
map_dummy_to_parent <- function(dat, features, selected_dummy) {
  trm <- terms(reformulate(features)); mm <- model.matrix(trm, dat); idx <- match(selected_dummy, colnames(mm)); idx <- idx[!is.na(idx)]
  if (!length(idx)) return(character()); unique(attr(trm,"term.labels")[attr(mm,"assign")[idx]])
}
select_features <- function(train_raw, seed) {
  set.seed(seed); x <- impute_for_selection(train_raw, FEATURES)
  bor <- Boruta::Boruta(Target ~ ., data = x, pValue = .01, mcAdj = TRUE, maxRuns = SIC_CONFIG$boruta_runs, doTrace = 0,
                        getImp = Boruta::getImpRfZ, num.threads = SIC_CONFIG$workers)
  bor <- Boruta::TentativeRoughFix(bor); bvars <- Boruta::getSelectedAttributes(bor, withTentative = FALSE)
  trm <- terms(reformulate(FEATURES)); mm <- model.matrix(trm, x)[,-1,drop=FALSE]; y <- as.integer(x$Target == "Class2")
  set.seed(seed + 1L); cv <- glmnet::cv.glmnet(mm, y, family = "binomial", alpha = 1, nfolds = max(3L,min(10L,SIC_CONFIG$inner_v)), type.measure = "auc", standardize = TRUE)
  co <- coef(cv, s = "lambda.1se"); dummies <- setdiff(rownames(co)[as.vector(co) != 0], "(Intercept)"); lvars <- map_dummy_to_parent(x, FEATURES, dummies)
  raw_intersection <- intersect(bvars, lvars); selected <- unique(c(FORCED_FEATURES, raw_intersection))
  list(boruta=bvars,lasso=lvars,intersection=raw_intersection,selected=selected,lambda_1se=cv$lambda.1se,
       boruta_decision=data.frame(Feature=names(bor$finalDecision),Decision=as.character(bor$finalDecision)))
}

set.seed(SIC_CONFIG$seed); outer <- rsample::vfold_cv(dev, v = SIC_CONFIG$outer_v, strata = Target)
checkpoint <- file.path(SIC_CONFIG$dirs$private, ifelse(SIC_CONFIG$fast_test,"feature_selection_debug.rds","feature_selection_checkpoint.rds"))
selection_list <- if (SIC_CONFIG$resume && file.exists(checkpoint)) readRDS(checkpoint) else list()
start_fold <- length(selection_list) + 1L
for (i in if (start_fold <= nrow(outer)) seq.int(start_fold,nrow(outer)) else integer()) {
  cat(sprintf("Feature selection outer fold %d/%d\n", i, nrow(outer)))
  fs <- select_features(rsample::analysis(outer$splits[[i]]), SIC_CONFIG$seed + i*100L)
  selection_list[[i]] <- tibble(Outer_Fold=i,Feature=FEATURES,Boruta=FEATURES%in%fs$boruta,LASSO=FEATURES%in%fs$lasso,
    Intersection=FEATURES%in%fs$intersection,Final_Selected=FEATURES%in%fs$selected,Forced=FEATURES%in%FORCED_FEATURES,Lambda_1se=fs$lambda_1se)
  saveRDS(selection_list, checkpoint)
}
selection_long <- bind_rows(selection_list)
frequency <- selection_long %>% group_by(Feature) %>% summarise(Boruta_Frequency=mean(Boruta),LASSO_Frequency=mean(LASSO),
  Intersection_Frequency=mean(Intersection),Forced=dplyr::first(Forced),.groups="drop") %>% arrange(desc(Intersection_Frequency),desc(Boruta_Frequency),desc(LASSO_Frequency))
final_features <- unique(c(FORCED_FEATURES, frequency %>% filter(Intersection_Frequency >= .70) %>% pull(Feature)))
if (length(final_features) < 5L) final_features <- unique(c(FORCED_FEATURES, head(frequency$Feature, 10L)))

## ============================================================================
## 4. Feature-selection figure (all candidate predictors)
## ============================================================================
display_map <- c(Age="Age",Gender="Gender",Race="Race",ICU_Type="ICU type",SOFA="SOFA",OASIS="OASIS",Charlson="Charlson",
  Hypertension="Hypertension",Diabetes="Diabetes",COPD="COPD",HF="Heart failure",Stroke="Stroke",Malignancy="Malignancy",
  Lung_Infection="Lung infection",GI_Infection="GI infection",GU_Infection="GU infection",Other_Infection="Other infection",
  MAP="MAP",Heart_Rate="Heart rate",Resp_Rate="Respiratory rate",Temperature="Temperature",Hb="Hemoglobin",WBC="WBC",
  Plt="Platelet count",Cre="Creatinine",AST="AST",TBil="Total bilirubin",Na="Sodium",K="Potassium",Ca="Calcium",Cl="Chloride",
  Plt_S="SIC platelet component",INR_S="SIC coagulation component",SOFA_S="SIC SOFA component")
plot_dat <- frequency %>% mutate(Display=unname(display_map[Feature]),Display=factor(Display,levels=rev(unname(display_map[frequency$Feature])))) %>%
  pivot_longer(c(Boruta_Frequency,LASSO_Frequency,Intersection_Frequency),names_to="Method",values_to="Frequency") %>%
  mutate(Method=recode(Method,Boruta_Frequency="Boruta",LASSO_Frequency="LASSO",Intersection_Frequency="Intersection"))
p_feature <- ggplot(plot_dat,aes(Frequency,Display,color=Method,shape=Method))+geom_vline(xintercept=.70,linetype=2,color="grey45",linewidth=.65)+
  geom_point(size=2.5,alpha=.95)+scale_x_continuous(limits=c(0,1),breaks=seq(0,1,.2),labels=percent_format(accuracy=1))+
  scale_color_manual(values=c(Boruta="#2CA25F",LASSO="#756BB1",Intersection="#DE2D26"))+
  labs(title="Nested feature-selection stability",x="Selection frequency across outer folds",y=NULL,color=NULL,shape=NULL)+
  theme_manuscript(10)+theme(legend.position="top",axis.text.y=element_text(size=7.8),plot.title=element_text(size=14))
save_png(p_feature,"Figure_ML1_Feature_Selection.png",8.3,8.8)

## ============================================================================
## 5. Export aggregate feature-engineering tables and private object
## ============================================================================
domain_map <- c(Age="Demographics",Gender="Demographics",Race="Demographics",ICU_Type="ICU context",SOFA="Severity",OASIS="Severity",Charlson="Comorbidity",
  Hypertension="Comorbidity",Diabetes="Comorbidity",COPD="Comorbidity",HF="Comorbidity",Stroke="Comorbidity",Malignancy="Comorbidity",
  Lung_Infection="Infection site",GI_Infection="Infection site",GU_Infection="Infection site",Other_Infection="Infection site",
  MAP="Vital signs",Heart_Rate="Vital signs",Resp_Rate="Vital signs",Temperature="Vital signs",Hb="Laboratory",WBC="Laboratory",Plt="Laboratory",
  Cre="Laboratory",AST="Laboratory",TBil="Laboratory",Na="Laboratory",K="Laboratory",Ca="Laboratory",Cl="Laboratory",
  Plt_S="SIC component",INR_S="SIC component",SOFA_S="SIC component")
feature_summary <- frequency %>% mutate(Clinical_Name=unname(display_map[Feature]),Domain=unname(domain_map[Feature]),
  Final_Model=Feature%in%final_features,Decision=case_when(Forced~"Prespecified and forced",Final_Model~"Stable intersection (>=70%)",TRUE~"Not retained")) %>%
  select(Feature,Clinical_Name,Domain,Boruta_Frequency,LASSO_Frequency,Intersection_Frequency,Forced,Final_Model,Decision)
fwrite(selection_long,file.path(SIC_CONFIG$dirs$tables,"feature_selection_by_fold.csv"),bom=TRUE)
fwrite(feature_summary,file.path(SIC_CONFIG$dirs$tables,"feature_engineering_summary.csv"),bom=TRUE)
fwrite(data.frame(Feature=final_features),file.path(SIC_CONFIG$dirs$tables,"final_features.csv"),bom=TRUE)
saveRDS(list(dev=dev,temporal=temporal,externals=externals,outer=outer,selection=selection_long,frequency=frequency,
             final_features=final_features,features=FEATURES,forced=FORCED_FEATURES),file.path(SIC_CONFIG$dirs$private,"feature_engineering.rds"))
cat("Feature engineering complete | final features:",paste(final_features,collapse=", "),"\n")
