## ============================================================================
## 1. Packages and configuration
## ============================================================================
if (!exists("SIC_CONFIG")) source(file.path(Sys.getenv("SIC_PROJECT_ROOT", unset=getwd()),"R","00_config.R"),encoding="UTF-8")
core_file <- file.path(SIC_CONFIG$dirs$private,"modeling_core.rds"); if(!file.exists(core_file)) stop("Run R/03_nested_modeling.R first.")

## ============================================================================
## 2. Import OOF and locked validation predictions
## ============================================================================
core <- readRDS(core_file); best_model <- core$best_model
oof <- core$oof %>% mutate(Validation_Set="Development nested OOF")
pred <- core$predictions
all_pred <- bind_rows(oof %>% select(stay_id,Cohort,Target,Model,Validation_Set,.pred_Class2),
                      pred %>% select(stay_id,Cohort,Target,Model,Validation_Set,.pred_Class2))

## ============================================================================
## 3. Performance, bootstrap CIs, DeLong tests, calibration, and DCA
## ============================================================================
safe_auc <- function(y,p) tryCatch(as.numeric(pROC::auc(pROC::roc(y,p,levels=c(0,1),direction="<",quiet=TRUE))),error=function(e) NA_real_)
safe_prauc <- function(y,p) tryCatch(PRROC::pr.curve(scores.class0=p[y==1],scores.class1=p[y==0],curve=FALSE)$auc.integral,error=function(e) NA_real_)
calibration_stats <- function(y,p) {
  p <- pmin(pmax(p,1e-6),1-1e-6); lp <- qlogis(p)
  ci <- tryCatch(coef(glm(y~offset(lp),family=binomial()))[1],error=function(e) NA_real_)
  cs <- tryCatch(coef(glm(y~lp,family=binomial()))[2],error=function(e) NA_real_)
  bins <- cut(p,breaks=unique(quantile(p,seq(0,1,.1),na.rm=TRUE)),include.lowest=TRUE)
  dd <- tibble(y=y,p=p,b=as.character(bins)) %>% filter(!is.na(b)) %>% group_by(b) %>% summarise(n=n(),obs=mean(y),pred=mean(p),.groups="drop")
  c(Calibration_Intercept=unname(ci),Calibration_Slope=unname(cs),Brier=mean((y-p)^2),
    ECE=if(nrow(dd)) weighted.mean(abs(dd$obs-dd$pred),dd$n) else NA_real_,MCE=if(nrow(dd)) max(abs(dd$obs-dd$pred)) else NA_real_)
}
class_stats <- function(y,p,t) {
  z<-as.integer(p>=t);tp<-as.numeric(sum(y==1&z==1));tn<-as.numeric(sum(y==0&z==0));fp<-as.numeric(sum(y==0&z==1));fn<-as.numeric(sum(y==1&z==0));div<-function(a,b)ifelse(b>0,a/b,NA_real_)
  den<-sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn));c(Sensitivity=div(tp,tp+fn),Specificity=div(tn,tn+fp),PPV=div(tp,tp+fp),NPV=div(tn,tn+fn),
    Accuracy=div(tp+tn,tp+tn+fp+fn),F1=div(2*tp,2*tp+fp+fn),MCC=ifelse(den>0,(tp*tn-fp*fn)/den,NA_real_))
}
thresholds <- oof %>% group_by(Model) %>% group_modify(~{
  y<-as.integer(.x$Target=="Class2");r<-pROC::roc(y,.x$.pred_Class2,levels=c(0,1),direction="<",quiet=TRUE)
  tibble(Threshold=as.numeric(pROC::coords(r,"best",best.method="youden",ret="threshold",transpose=FALSE)[1]))
}) %>% ungroup()
metric_one <- function(df,t,seed) {
  y<-as.integer(df$Target=="Class2");p<-df$.pred_Class2;point<-c(AUROC=safe_auc(y,p),AUPRC=safe_prauc(y,p),calibration_stats(y,p),class_stats(y,p,t))
  set.seed(seed);n<-length(y);boots<-replicate(SIC_CONFIG$bootstrap_n,{ii<-sample.int(n,n,replace=TRUE);c(AUROC=safe_auc(y[ii],p[ii]),AUPRC=safe_prauc(y[ii],p[ii]),Brier=mean((y[ii]-p[ii])^2))})
  ci<-apply(boots,1,quantile,c(.025,.975),na.rm=TRUE);lower<-upper<-setNames(rep(NA_real_,length(point)),names(point));idx<-match(names(point),colnames(ci));ok<-!is.na(idx)
  lower[ok]<-ci[1,idx[ok]];upper[ok]<-ci[2,idx[ok]];tibble(Metric=names(point),Estimate=as.numeric(point),Lower=as.numeric(lower),Upper=as.numeric(upper))
}
performance <- all_pred %>% left_join(thresholds,by="Model") %>% group_by(Model,Validation_Set) %>%
  group_modify(~metric_one(.x,unique(.x$Threshold),SIC_CONFIG$seed+sum(utf8ToInt(paste(.y$Model,.y$Validation_Set))))) %>% ungroup()
performance_wide <- performance %>% select(Model,Validation_Set,Metric,Estimate) %>% pivot_wider(names_from=Metric,values_from=Estimate)

base <- oof %>% filter(Model==best_model) %>% select(stay_id,Target,p_best=.pred_Class2)
delong <- bind_rows(lapply(setdiff(unique(oof$Model),best_model),function(nm){
  z<-oof%>%filter(Model==nm)%>%select(stay_id,p_other=.pred_Class2)%>%inner_join(base,by="stay_id");y<-as.integer(z$Target=="Class2")
  r1<-pROC::roc(y,z$p_best,quiet=TRUE);r2<-pROC::roc(y,z$p_other,quiet=TRUE)
  tibble(Reference=best_model,Comparator=nm,Delta_AUROC=as.numeric(pROC::auc(r1)-pROC::auc(r2)),P_value=pROC::roc.test(r1,r2,paired=TRUE,method="delong")$p.value)
})) %>% mutate(P_FDR=p.adjust(P_value,"BH"))

best_pred <- pred %>% filter(Model==best_model)
cal_data <- best_pred %>% group_by(Validation_Set) %>% group_modify(~{
  q<-unique(quantile(.x$.pred_Class2,seq(0,1,.1),na.rm=TRUE));if(length(q)<3)q<-seq(0,1,length.out=11)
  .x%>%mutate(Bin=cut(.pred_Class2,breaks=q,include.lowest=TRUE))%>%group_by(Bin)%>%summarise(N=n(),Predicted=mean(.pred_Class2),Observed=mean(Target=="Class2"),.groups="drop")
}) %>% ungroup()
dca_data <- best_pred %>% group_by(Validation_Set) %>% group_modify(~{
  y<-as.integer(.x$Target=="Class2");p<-.x$.pred_Class2;n<-length(y);prev<-mean(y)
  map_dfr(seq(.05,.80,.01),function(t){z<-p>=t;tibble(Threshold=t,Model_NB=sum(y==1&z)/n-sum(y==0&z)/n*t/(1-t),Treat_All_NB=prev-(1-prev)*t/(1-t),Treat_None_NB=0)})
}) %>% ungroup() %>% pivot_longer(c(Model_NB,Treat_All_NB,Treat_None_NB),names_to="Strategy",values_to="Net_Benefit") %>%
  mutate(Strategy=recode(Strategy,Model_NB=best_model,Treat_All_NB="Treat all",Treat_None_NB="Treat none"))

## ============================================================================
## 4. Publication figures: 4 ROC panels, DCA, and calibration
## ============================================================================
model_colors <- setNames(viridisLite::plasma(length(unique(pred$Model)),end=.86),sort(unique(pred$Model)))
cohort_colors <- c("MIMIC-IV temporal"="#377EB8","MIMIC-III"="#E41A1C","eICU"="#FF7F00","NWICU"="#4DAF4A")
make_roc <- function(ds,file,title) {
  d <- pred %>% filter(Validation_Set==ds)
  curves <- d %>% group_by(Model) %>% group_modify(~{y<-as.integer(.x$Target=="Class2");r<-pROC::roc(y,.x$.pred_Class2,quiet=TRUE);tibble(FPR=1-r$specificities,TPR=r$sensitivities,AUC=as.numeric(pROC::auc(r)))}) %>% ungroup() %>%
    mutate(Label=sprintf("%s (AUC = %.3f)",Model,AUC))
  labs <- curves %>% distinct(Model,Label,AUC) %>% arrange(desc(AUC))
  p <- ggplot(curves,aes(FPR,TPR,color=Model))+geom_abline(slope=1,intercept=0,linetype=2,color="grey60")+geom_line(linewidth=1.05)+coord_equal(xlim=c(0,1),ylim=c(0,1))+
    scale_color_manual(values=model_colors,breaks=labs$Model,labels=labs$Label)+labs(title=title,x="1 - Specificity",y="Sensitivity",color=NULL)+theme_manuscript(10)+
    theme(legend.position=c(.68,.25),legend.justification=c(.5,.5),legend.background=element_rect(fill=alpha("white",.90),color="grey70"),legend.text=element_text(size=7.5),plot.title=element_text(size=13))
  save_png(p,file,6.6,5.8)
}
make_roc("MIMIC-IV temporal","Figure_ML2A_ROC_MIMICIV_Temporal.png","MIMIC-IV temporal validation")
make_roc("MIMIC-III","Figure_ML2B_ROC_MIMICIII.png","MIMIC-III external validation")
make_roc("eICU","Figure_ML2C_ROC_eICU.png","eICU external validation")
make_roc("NWICU","Figure_ML2D_ROC_NWICU.png","NWICU external validation")

p_dca <- ggplot(dca_data,aes(Threshold,Net_Benefit,color=Strategy,linetype=Strategy))+geom_line(linewidth=.85)+facet_wrap(~Validation_Set,ncol=2,scales="free_y")+
  scale_color_manual(values=c(setNames("#D62728",best_model),`Treat all`="grey35",`Treat none`="grey70"))+scale_linetype_manual(values=c(setNames("solid",best_model),`Treat all`="dashed",`Treat none`="dotted"))+
  labs(title=paste0(best_model," decision-curve analysis"),x="Threshold probability",y="Net benefit",color=NULL,linetype=NULL)+theme_manuscript(9)+
  theme(legend.position=c(.80,.13),legend.background=element_rect(fill=alpha("white",.90),color="grey75"),strip.text=element_text(face="bold"),plot.title=element_text(size=12))
save_png(p_dca,"Figure_ML2E_XGBoost_DCA.png",8.3,6.4)

p_cal <- ggplot(cal_data,aes(Predicted,Observed,color=Validation_Set))+geom_abline(slope=1,intercept=0,linetype=2,color="grey60")+geom_line(linewidth=.9)+geom_point(size=1.8)+
  scale_color_manual(values=cohort_colors)+coord_equal(xlim=c(0,1),ylim=c(0,1))+labs(title=paste0(best_model," calibration"),x="Mean predicted probability",y="Observed Class 2 proportion",color=NULL)+theme_manuscript(10)+
  theme(legend.position=c(.30,.80),legend.background=element_rect(fill=alpha("white",.90),color="grey75"),legend.text=element_text(size=8),plot.title=element_text(size=12))
save_png(p_cal,"Figure_ML2F_XGBoost_Calibration.png",6.6,5.8)

## ============================================================================
## 5. Export aggregate evaluation tables
## ============================================================================
fwrite(thresholds,file.path(SIC_CONFIG$dirs$tables,"development_frozen_thresholds.csv"),bom=TRUE)
fwrite(performance,file.path(SIC_CONFIG$dirs$tables,"performance_with_confidence_intervals.csv"),bom=TRUE)
fwrite(performance_wide,file.path(SIC_CONFIG$dirs$tables,"performance_summary.csv"),bom=TRUE)
fwrite(delong,file.path(SIC_CONFIG$dirs$tables,"delong_comparisons_fdr.csv"),bom=TRUE)
fwrite(cal_data,file.path(SIC_CONFIG$dirs$tables,"calibration_plot_data.csv"),bom=TRUE)
fwrite(dca_data,file.path(SIC_CONFIG$dirs$tables,"dca_plot_data.csv"),bom=TRUE)
saveRDS(list(performance=performance,performance_wide=performance_wide,thresholds=thresholds,delong=delong,calibration=cal_data,dca=dca_data),file.path(SIC_CONFIG$dirs$private,"evaluation_results.rds"))
cat("Model evaluation complete | locked model:",best_model,"\n")
