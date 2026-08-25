## ============================================================================
## 1. Packages and configuration
## ============================================================================
if (!exists("SIC_CONFIG")) source(file.path(Sys.getenv("SIC_PROJECT_ROOT", unset=getwd()),"R","00_config.R"),encoding="UTF-8")
core_file <- file.path(SIC_CONFIG$dirs$private,"modeling_core.rds"); eval_file <- file.path(SIC_CONFIG$dirs$private,"evaluation_results.rds")
if(!file.exists(core_file)||!file.exists(eval_file)) stop("Run model evaluation before SHAP analysis.")

## ============================================================================
## 2. Import locked model and validation cohorts
## ============================================================================
core <- readRDS(core_file); evaluation <- readRDS(eval_file); best_model <- core$best_model; fit_best <- core$final_fits[[best_model]]; features <- core$final_features
sets <- c(list("MIMIC-IV temporal"=core$temporal),core$externals)
threshold <- evaluation$thresholds$Threshold[match(best_model,evaluation$thresholds$Model)]
pred_wrapper <- function(object,newdata) as.numeric(predict(object,new_data=as.data.frame(newdata),type="prob")$.pred_Class2)
set.seed(SIC_CONFIG$seed+300000L); background <- core$dev %>% slice_sample(n=min(200,nrow(core$dev))) %>% select(all_of(features))

## ============================================================================
## 3. Cross-cohort global and representative local SHAP analysis
## ============================================================================
shap_objects <- list(); importance <- list(); sampled_rows <- list(); shap_matrices <- list()
for(i in seq_along(sets)) {
  nm <- names(sets)[i]; dat <- sets[[i]]; set.seed(SIC_CONFIG$seed+300100L+i); samp <- dat %>% slice_sample(n=min(SIC_CONFIG$shap_sample,nrow(dat)))
  X <- samp %>% select(all_of(features)); cat("SHAP:",nm,"n =",nrow(X),"\n");set.seed(SIC_CONFIG$seed+300200L+i)
  S <- fastshap::explain(fit_best,X=background,pred_wrapper=pred_wrapper,newdata=X,nsim=SIC_CONFIG$shap_nsim,adjust=TRUE,shap_only=TRUE)
  shap_matrices[[nm]] <- S; shap_objects[[nm]] <- shapviz::shapviz(S,X=X); sampled_rows[[nm]] <- samp
  importance[[nm]] <- tibble(Cohort=nm,Feature=colnames(S),MeanAbsSHAP=colMeans(abs(S),na.rm=TRUE),MeanSHAP=colMeans(S,na.rm=TRUE)) %>% arrange(desc(MeanAbsSHAP)) %>% mutate(Rank=row_number())
}
importance_all <- bind_rows(importance); ref <- importance_all %>% filter(Cohort=="MIMIC-IV temporal") %>% select(Feature,Ref_Rank=Rank)
stability <- importance_all %>% left_join(ref,by="Feature") %>% group_by(Cohort) %>% summarise(Spearman_Rank=cor(Rank,Ref_Rank,method="spearman",use="complete.obs"),
  Top10_Overlap=sum(Feature[Rank<=10]%in%Feature[Ref_Rank<=10]),.groups="drop")

tmp <- sampled_rows[["MIMIC-IV temporal"]]; tmp$Probability <- pred_wrapper(fit_best,tmp[,features,drop=FALSE]); tmp$Predicted <- tmp$Probability>=threshold
types <- list(TP=which(tmp$Target=="Class2"&tmp$Predicted),FP=which(tmp$Target=="Other"&tmp$Predicted),FN=which(tmp$Target=="Class2"&!tmp$Predicted))
local_public <- list(); local_private <- list()
for(tp in names(types)) {
  ids<-types[[tp]];if(!length(ids))next;med<-median(tmp$Probability[ids]);id<-ids[which.min(abs(tmp$Probability[ids]-med))]
  local_public[[tp]]<-tibble(Case_ID=paste0("Representative_",tp),Type=tp,Observed=as.character(tmp$Target[id]),Probability=tmp$Probability[id],Threshold=threshold)
  local_private[[tp]]<-tibble(Case_ID=paste0("Representative_",tp),stay_id=tmp$stay_id[id],Type=tp,Observed=as.character(tmp$Target[id]),Probability=tmp$Probability[id],Threshold=threshold)
}

## ============================================================================
## 4. Publication SHAP importance and beeswarm figures
## ============================================================================
cohort_colors <- c("MIMIC-IV temporal"="#377EB8","MIMIC-III"="#E41A1C","eICU"="#FF7F00","NWICU"="#4DAF4A")
top15 <- importance_all %>% filter(Cohort=="MIMIC-IV temporal") %>% slice_head(n=min(15,length(features))) %>% pull(Feature)
p_importance <- importance_all %>% filter(Feature%in%top15) %>% mutate(Feature=factor(Feature,levels=rev(top15))) %>%
  ggplot(aes(MeanAbsSHAP,Feature,fill=Cohort))+geom_col(position=position_dodge(width=.80),width=.72)+scale_fill_manual(values=cohort_colors)+
  labs(title="Cross-cohort SHAP feature importance",x="Mean absolute SHAP value",y=NULL,fill=NULL)+theme_manuscript(9)+
  theme(legend.position=c(.72,.25),legend.background=element_rect(fill=alpha("white",.90),color="grey75"),legend.text=element_text(size=7.5),
        axis.text.y=element_text(face="bold",size=7.6),plot.title=element_text(size=11.5))
save_png(p_importance,"Figure_ML2G_XGBoost_SHAP_Importance.png",7.2,5.8)

p_beeswarm <- shapviz::sv_importance(shap_objects[["MIMIC-IV temporal"]],kind="beeswarm",max_display=min(20,length(features)))+
  labs(title=paste0(best_model," SHAP summary"),x=paste0("SHAP value for ",best_model," prediction"),y=NULL)+theme_manuscript(9)+
  theme(axis.text.y=element_text(face="bold",size=7.3),axis.text.x=element_text(size=8),plot.title=element_text(size=11.5),
        legend.position=c(.82,.18),legend.background=element_rect(fill=alpha("white",.90),color="grey75"),legend.title=element_text(size=7.5),legend.text=element_text(size=7))
save_png(p_beeswarm,"Figure_ML2H_XGBoost_SHAP_Beeswarm.png",7.2,5.8)

## ============================================================================
## 5. Export aggregate SHAP tables and private calculator bundle
## ============================================================================
fwrite(importance_all,file.path(SIC_CONFIG$dirs$tables,"shap_global_importance.csv"),bom=TRUE)
fwrite(stability,file.path(SIC_CONFIG$dirs$tables,"shap_cross_cohort_stability.csv"),bom=TRUE)
fwrite(bind_rows(local_public),file.path(SIC_CONFIG$dirs$tables,"shap_representative_cases.csv"),bom=TRUE)
fwrite(bind_rows(local_private),file.path(SIC_CONFIG$dirs$private,"shap_representative_cases_with_ids.csv"),bom=TRUE)
defaults <- lapply(features,function(v){x<-core$dev[[v]];if(is.numeric(x))list(type="numeric",value=median(x,na.rm=TRUE),min=quantile(x,.01,na.rm=TRUE),max=quantile(x,.99,na.rm=TRUE))else list(type="categorical",value=mode_value(x),choices=levels(x))});names(defaults)<-features
bundle <- list(model=fit_best,model_name=best_model,features=features,background=background,threshold=threshold,defaults=defaults,
  global_importance=importance[["MIMIC-IV temporal"]],created=as.character(Sys.time()),seed=core$seed)
saveRDS(bundle,file.path(SIC_CONFIG$dirs$private,"model_bundle.rds"));saveRDS(list(shap=shap_objects,matrices=shap_matrices),file.path(SIC_CONFIG$dirs$private,"shap_objects.rds"))
cat("SHAP analysis complete | public tables contain no patient identifiers.\n")
