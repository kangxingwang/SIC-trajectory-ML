## ============================================================================
## 1. Packages and configuration
## ============================================================================
if (!exists("SIC_CONFIG")) source(file.path(Sys.getenv("SIC_PROJECT_ROOT", unset = getwd()), "R", "00_config.R"), encoding = "UTF-8")
feature_file <- file.path(SIC_CONFIG$dirs$private,"feature_engineering.rds")
if (!file.exists(feature_file)) stop("Run R/02_feature_engineering.R first.")

## ============================================================================
## 2. Import fold-specific feature sets and model data
## ============================================================================
fe <- readRDS(feature_file); dev <- fe$dev; temporal <- fe$temporal; externals <- fe$externals; outer <- fe$outer
final_features <- fe$final_features
MODELS_TO_RUN <- if (SIC_CONFIG$fast_test) c("ElasticNet","XGBoost") else c("ElasticNet","RandomForest","XGBoost","SVM_RBF","KNN")

make_recipe <- function(dat, features, normalize = FALSE) {
  rec <- recipe(Target ~ ., data = dat[,c("Target",features),drop=FALSE]) %>%
    step_novel(all_nominal_predictors(),new_level="Novel") %>% step_unknown(all_nominal_predictors(),new_level="Unknown") %>%
    step_impute_median(all_numeric_predictors()) %>% step_zv(all_predictors())
  if (normalize) rec <- rec %>% step_YeoJohnson(all_numeric_predictors(),num_unique=5) %>% step_normalize(all_numeric_predictors())
  rec %>% step_dummy(all_nominal_predictors(),one_hot=TRUE)
}
model_spec <- function(name) switch(name,
  ElasticNet=logistic_reg(penalty=tune(),mixture=tune()) %>% set_engine("glmnet") %>% set_mode("classification"),
  RandomForest=rand_forest(mtry=tune(),min_n=tune(),trees=800) %>% set_engine("ranger",importance="permutation",probability=TRUE,num.threads=1) %>% set_mode("classification"),
  XGBoost=boost_tree(trees=tune(),tree_depth=tune(),learn_rate=tune(),loss_reduction=tune(),min_n=tune(),sample_size=tune(),mtry=tune()) %>% set_engine("xgboost",nthread=1,verbosity=0) %>% set_mode("classification"),
  SVM_RBF=svm_rbf(cost=tune(),rbf_sigma=tune()) %>% set_engine("kernlab") %>% set_mode("classification"),
  KNN=nearest_neighbor(neighbors=tune(),dist_power=tune(),weight_func=tune()) %>% set_engine("kknn") %>% set_mode("classification"))
needs_norm <- function(name) name %in% c("ElasticNet","SVM_RBF","KNN")
tune_fit_predict <- function(name, train_raw, val_raw, features, seed) {
  rec <- make_recipe(train_raw,features,needs_norm(name)); wf <- workflow() %>% add_recipe(rec) %>% add_model(model_spec(name))
  set.seed(seed); inner <- vfold_cv(train_raw,v=SIC_CONFIG$inner_v,strata=Target)
  ctrl <- control_grid(save_pred=FALSE,save_workflow=TRUE,parallel_over="everything",verbose=FALSE)
  set.seed(seed+1L); tuned <- tune_grid(wf,resamples=inner,grid=SIC_CONFIG$grid_n,metrics=metric_set(roc_auc,pr_auc),control=ctrl)
  best <- select_best(tuned,metric="roc_auc"); final_wf <- finalize_workflow(wf,best); fitted <- fit(final_wf,data=train_raw)
  pred <- predict(fitted,new_data=val_raw,type="prob") %>% bind_cols(val_raw %>% select(stay_id,Cohort,Target))
  list(pred=pred,best=best,fit=fitted,cv=collect_metrics(tuned))
}

## ============================================================================
## 3. Nested tuning, OOF model ranking, and locked external validation
## ============================================================================
future::plan(future::multisession,workers=SIC_CONFIG$workers); on.exit(future::plan(future::sequential),add=TRUE)
checkpoint_file <- file.path(SIC_CONFIG$dirs$private,ifelse(SIC_CONFIG$fast_test,"nested_modeling_debug.rds","nested_modeling_checkpoint.rds"))
if (SIC_CONFIG$resume && file.exists(checkpoint_file)) {
  cp <- readRDS(checkpoint_file); oof_list <- cp$oof; best_list <- cp$best; cv_list <- cp$cv; start_outer <- cp$completed + 1L
} else { oof_list <- list(); best_list <- list(); cv_list <- list(); start_outer <- 1L }

for (i in if(start_outer<=nrow(outer)) seq.int(start_outer,nrow(outer)) else integer()) {
  tr <- analysis(outer$splits[[i]]); va <- assessment(outer$splits[[i]])
  fold_features <- fe$selection %>% filter(Outer_Fold==i,Final_Selected) %>% pull(Feature)
  if (!length(fold_features)) stop("No selected features in outer fold ",i)
  for (j in seq_along(MODELS_TO_RUN)) {
    nm <- MODELS_TO_RUN[j]; cat(sprintf("Outer %d/%d | tuning %s\n",i,nrow(outer),nm))
    ans <- tune_fit_predict(nm,tr,va,fold_features,SIC_CONFIG$seed+i*1000L+j*10L)
    key <- paste(i,nm,sep="__"); oof_list[[key]] <- ans$pred %>% mutate(Model=nm,Outer_Fold=i)
    best_list[[key]] <- ans$best %>% mutate(Model=nm,Outer_Fold=i); cv_list[[key]] <- ans$cv %>% mutate(Model=nm,Outer_Fold=i)
  }
  saveRDS(list(oof=oof_list,best=best_list,cv=cv_list,completed=i),checkpoint_file)
}
oof <- bind_rows(oof_list); best_outer <- bind_rows(best_list); cv_outer <- bind_rows(cv_list)
model_rank <- oof %>% group_by(Model) %>% summarise(AUROC=roc_auc_vec(Target,.pred_Class2,event_level="first"),
  AUPRC=pr_auc_vec(Target,.pred_Class2,event_level="first"),.groups="drop") %>% arrange(desc(AUROC),desc(AUPRC))
best_model <- model_rank$Model[1]

set.seed(SIC_CONFIG$seed+90000L); inner_final <- vfold_cv(dev,v=SIC_CONFIG$inner_v,strata=Target)
final_fits <- list(); final_tuning <- list(); final_best <- list(); predictions <- list()
for (j in seq_along(MODELS_TO_RUN)) {
  nm <- MODELS_TO_RUN[j]; cat("Final development tuning:",nm,"\n")
  rec <- make_recipe(dev,final_features,needs_norm(nm)); wf <- workflow() %>% add_recipe(rec) %>% add_model(model_spec(nm))
  set.seed(SIC_CONFIG$seed+91000L+j); tuned <- tune_grid(wf,resamples=inner_final,grid=SIC_CONFIG$grid_n,metrics=metric_set(roc_auc,pr_auc),
    control=control_grid(save_pred=TRUE,save_workflow=TRUE,parallel_over="everything",verbose=FALSE))
  bst <- select_best(tuned,metric="roc_auc"); fitted <- fit(finalize_workflow(wf,bst),data=dev)
  final_fits[[nm]] <- fitted; final_tuning[[nm]] <- collect_metrics(tuned); final_best[[nm]] <- bst
  datasets <- c(list("MIMIC-IV temporal"=temporal),externals)
  for (ds in names(datasets)) predictions[[paste(nm,ds,sep="__")]] <- predict(fitted,new_data=datasets[[ds]],type="prob") %>%
    bind_cols(datasets[[ds]] %>% select(stay_id,Cohort,Target)) %>% mutate(Model=nm,Validation_Set=ds)
}
pred_all <- bind_rows(predictions)
final_best_tbl <- bind_rows(lapply(names(final_best),function(nm) final_best[[nm]] %>% mutate(Model=nm,.before=1)))

## ============================================================================
## 4. Model-ranking plot
## ============================================================================
p_rank <- model_rank %>% mutate(Model=factor(Model,levels=rev(Model))) %>%
  ggplot(aes(AUROC,Model,color=Model))+geom_segment(aes(x=.5,xend=AUROC,y=Model,yend=Model),color="grey80",linewidth=1)+geom_point(size=3)+
  geom_text(aes(label=sprintf("%.3f",AUROC)),hjust=-.25,size=3.4,color="black")+scale_x_continuous(limits=c(.5,max(.9,max(model_rank$AUROC)+.06)))+
  scale_color_viridis_d(option="plasma",end=.85)+labs(title="Nested out-of-fold model performance",x="AUROC",y=NULL)+theme_manuscript(10)+theme(legend.position="none")
save_png(p_rank,"Supplement_Model_Ranking.png",7.2,4.8)

## ============================================================================
## 5. Export aggregate tuning tables and private model objects
## ============================================================================
fwrite(model_rank,file.path(SIC_CONFIG$dirs$tables,"nested_oof_model_ranking.csv"),bom=TRUE)
fwrite(best_outer,file.path(SIC_CONFIG$dirs$tables,"outer_fold_best_hyperparameters.csv"),bom=TRUE)
fwrite(final_best_tbl,file.path(SIC_CONFIG$dirs$tables,"final_best_hyperparameters.csv"),bom=TRUE)
fwrite(cv_outer,file.path(SIC_CONFIG$dirs$tables,"inner_cv_metrics.csv"),bom=TRUE)
fwrite(pred_all,file.path(SIC_CONFIG$dirs$private,"locked_validation_predictions.csv"),bom=TRUE)
saveRDS(list(seed=SIC_CONFIG$seed,outer_v=SIC_CONFIG$outer_v,inner_v=SIC_CONFIG$inner_v,models=MODELS_TO_RUN,dev=dev,temporal=temporal,
  externals=externals,final_features=final_features,oof=oof,model_rank=model_rank,best_model=best_model,final_fits=final_fits,
  final_tuning=final_tuning,final_best=final_best_tbl,predictions=pred_all),file.path(SIC_CONFIG$dirs$private,"modeling_core.rds"))
cat("Nested modeling complete | locked model:",best_model,"| features:",paste(final_features,collapse=", "),"\n")

