rm(list=ls())
library(lubridate)
library(corrplot)
library(pheatmap)
library(epiR)
library(readxl)
library(tidyr)
library(dplyr)
library(data.table)
library(stringr)
library(table1)
library(ggpubr)
library(purrr)
library(VennDiagram)
library(RColorBrewer)
library(scales)
library(Seurat)
library(cowplot)
library(patchwork)
library(psych)
library(PerformanceAnalytics)
library(rcompanion)
library(lmtest)
library("survival")
library("survminer")
library(gbm)
library(scales)
library(EnhancedVolcano)
library(umap)
library("FactoMineR")
library("factoextra")
library(gridExtra)
library(gplots)
library(caret)
library(plotROC)
library(MLmetrics)
library(MLeval)
library(forestplot)
library(scatterplot3d)
data_concater <- function(x){
  x<- levels(factor(x))
  paste(x, collapse = "+")
}
my.render.cont <- function(x) {
  with(stats.apply.rounding(stats.default(x), digits=2), 
       c("", "Mean (SD)"=sprintf("%s (&plusmn; %s)", MEAN, SD)))
}
my.render.cat <- function(x) {
  c("", sapply(stats.default(x), function(y) with(y, sprintf("%d (%0.0f %%)", FREQ, PCT))))
}
rndr.strat <- function(label, n, ...) {
  ifelse(n==0, label, render.strat.default(label, n, ...))
}
collapse<-function(x){
  paste(unique(x),collapse = ",")
}
cox.stats<-function(x){ 
  if(!is(x,"try-error")){
    x <- summary(x)
    p.value<-signif(x$wald["pvalue"], digits=2)
    wald.test<-signif(x$wald["test"], digits=2)
    beta<-signif(x$coef[1], digits=2);#coeficient beta
    HR <-signif(x$coef[2], digits=2);#exp(beta)
    HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
    HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
    HR <- paste0(HR, " (", 
                 HR.confint.lower, "-", HR.confint.upper, ")")
    res<-c(beta, HR, wald.test, p.value)
    names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                  "p.value")
    return(res)
    #return(exp(cbind(coef(x),confint(x))))
  }
}
plot.umap<-function(x, labels,
         main="UMAP",
         colors=c("#ff7f00", "#e377c2", "#17becf"),
         pad=0.1, cex=0.6, pch=19, add=FALSE, legend.suffix="",
         cex.main=1, cex.legend=0.85) {
  layout = x
  if (is(x, "umap")) {
    layout = x$layout
  }
  xylim = range(layout)
  xylim = xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)
  if (!add) {
    par(mar=c(0.2,0.7,1.2,0.7), ps=10)
    plot(xylim, xylim, type="n", axes=F, frame=F)
    rect(xylim[1], xylim[1], xylim[2], xylim[2], border="#aaaaaa", lwd=0.25)
  }
  points(layout[,1], layout[,2], col=colors[as.integer(labels)],
         cex=cex, pch=pch)
  mtext(side=3, main, cex=cex.main)
  labels.u = unique(labels)
  legend.pos = "topleft"
  legend.text = as.character(labels.u)
  if (add) {
    legend.pos = "bottomleft"
    legend.text = paste(as.character(labels.u), legend.suffix)
  }
  legend(legend.pos, legend=legend.text, inset=0.03,
         col=colors[as.integer(labels.u)],
         bty="n", pch=pch, cex=cex.legend)
}
save_pheatmap_png <- function(x, filename, width,height, res = 400) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

dir.create("stats")

#load data
df.test<-read_excel("deID.clem.xlsx",sheet="data")

#####ANALYSIS

#define variables of interest
df<-df.test
outcome.num<-unique(c("tempdurhome","ptc_ox_sat","ptc_gcs","ptc_supp_liters",
               "vs_o2vol","vs_fio2","oxygen_saturation_spo2_v2","ventdays","icu_days",
               "pd_calc_icu_days","icu.count","ordscore","length_stay","time_to_fu","enc.count",
               grep(paste(c("floor","min","max","avg"),collapse="|"),colnames(df),value=T)))
vars.num<-c("dem_age","a1c","avg.bmi","age_year","number_immuno_drugs",grep(paste(c("min","max","avg"),collapse="|"),colnames(df),value=T))
outcome.cat<-unique(c("outcome","ordscore","hospvent","patmanicudt","vvecmo","corangio","dschstat","ever_hosp","oc_newstatus","ptc_entry","ptc_present_intubate_yn",
               "ptc_ox_supp_yn","ptc_intubate_yn","causdth","loc","death","vap","ecmo","rrt","were_infiltrates_present",
               "disch_disp_full","poe_o2","poe_vent","inout_cd","ibax_disch_disp desc","death_yn",
               grep(paste(c("shk","discharge","oxygenation","comp_","cxr_pattern","ct_chest","lung_ultrasound"),collapse="|"),colnames(df),value=T)))
vars.cat<-unique(c("ptc_immuno_drug","ptc_immuno_drug_core","ckd_stages","dem_sex","dem_gender","race.x","ethnicity.x","sub_etoh","sub_drug","sub_cig","sub_ends","sub_mj",
            "ptc_insteroid_dur","ptc_orsteroid_dur","prone","if_yes_steroid_route_of_ad_v2","zipcode.x","agegroup","med.immunomod.drugs_yn","Tobacco use",
            colnames(comorb)[-1],"pat_blood_type","pat_rh","gender","race.y","ethnicity.y","ethnic_cd desc","hispanic_ind","zipcode.y","gender_MF","dem_preg",
            grep(paste(c("ptc_tx","pmh","ptc_sx__","path","blood_products","medications__","antiviral__","lab.abnormal","_yn"),collapse="|"),colnames(df),value=T)))

vars.select<-unique(c("ptc_ox_sat","oxygen_saturation_spo2_v2","ventdays","icu_days","length_stay","tmax","enc.count",
                      grep(paste(c("min_","max_","avg_"),collapse="|"),colnames(df),value=T),"number_immuno_drugs",
                      "pd_calc_icu_days","icu.count","ordscore","vvecmo","dschstat","ever_hosp","oc_newstatus","loc","death","vap","ecmo",
                      "disch_disp_full","poe_o2","poe_vent","inout_cd","ibax_disch_disp desc","death_yn",
                      "agegroup","gender_MF","race.y","zipcode.y","dem_age","a1c","avg.bmi","dem_age",
                      colnames(comorb)[-1],grep(paste(c("immun"),collapse="|"),colnames(df),value=T),
                      grep("immun",colnames(pmh),value=T),colnames(med)))


#demographic and outcomes tables (with stats) 
i="death_yn"
df2<-df[!is.na(df[[i]]),]
rndr <- function(x, name, ...) {
  if (length(x) == 0) {
    y <- df2[[name]]
    s <- rep("", length(render.default(x=y, name=name, ...)))
    if (is.numeric(y)) {
      p <- wilcox.test(y ~ df2[[i]])$p.value #t.test
    } else {
      p <- chisq.test(table(y, droplevels(df2[[i]])))$p.value
    }
    s[2] <- sub("<", "&lt;", format.pval(p, digits=3, eps=0.001))
    s
  } else {
    render.default(x=x, name=name, ...)
  }
}
df2[[i]]<- factor(df2[[i]], levels=c("no","yes","p"), labels=c("No", "Yes", "P-value"))
table1(~ gender_MF+dem_age+race.y+pat_blood_type+avg.bmi| df2[[i]],
       data=df2, droplevels=F, render=rndr, render.strat=rndr.strat, overall=T,
       render.continuous=my.render.cont2, render.categorical=my.render.cat)%>%
  cat(., file =paste0(i,".dem.pval.table1.html"))
df3<-df2[!is.na(df2$htn_yn),]
table1(~ htn_yn+CKD_yn+diabetes_yn+obesity_yn+
         rheumatologic_yn+autoimmune_yn+malignancy_yn+immunosup_yn+copd_yn+asthma_yn+cad_yn+cvd_yn| df3[[i]],
       data=df3, droplevels=F, render=rndr, render.strat=rndr.strat, overall=T,
       render.continuous=my.render.cont2, render.categorical=my.render.cat)%>%
  cat(., file =paste0(i,".comorb.pval.table1.html"))
table1(~ `mean.Absolute Lymphocyte Count`+`mean.C-Reactive Protein`+mean.Creatinine+mean.Ferritin+`mean.D-Dimer`+
         `mean.Creatine Kinase (CK)`+`mean.INR(PT)`+`mean.Lactate Dehydrogenase (LD)`+mean.pH+`mean.Platelet Count`+
         mean.PT+mean.PTT+`mean.Absolute Neutrophil Count`+a1c | df2[[i]],
       data=df2, droplevels=F, render=rndr, render.strat=rndr.strat, overall=T,
       render.continuous=my.render.cont2, render.categorical=my.render.cat)%>%
  cat(., file =paste0(i,".labs.pval.table1.html"))
table1(~ avg_RR+avg_HR+tmax+min_SBP+min_DBP | df2[[i]],
       data=df2, droplevels=F, render=rndr, render.strat=rndr.strat, overall=T,
       render.continuous=my.render.cont2, render.categorical=my.render.cat)%>%
  cat(., file =paste0(i,".vitals.pval.table1.html"))
df3<-df2[!is.na(df2$med.corticosteroid_yn),]
table1(~ med.corticosteroid_yn+med.calcineurin_yn+med.antirheumatic_yn+immunosup.drug_yn+med.nsaid.drugs_yn+#med.card.drugs_yn+
         med.onc.drugs_yn+med.antiglycemic.drugs_yn+med.asthma.drugs_yn+
         med.biologics.drug_yn+med.bone.drugs_yn+med.opioid.drugs_yn+med.antihistamine.drugs_yn+med.gout.drugs_yn+med.htn.drugs_yn | df3[[i]],
       data=df3, droplevels=F, render=rndr, render.strat=rndr.strat, overall=T,
       render.continuous=my.render.cont2, render.categorical=my.render.cat)%>%
  cat(., file =paste0(i,".meds.pval.table1.html"))
df2$patmanicudt<-as.factor(df2$patmanicudt)
df2$icu_yn<-"no"
df2$icu_yn[df2$patmanicudt==1]<-"yes"
table1(~ tmax + comp_sum___rrt+poe_o2+poe_vent+enc.count+inout_cd+length_stay+ordscore+patmanicudt+icu_yn+icu_days| df2[[i]],
       data=df2, droplevels=F, render=rndr, render.strat=rndr.strat, overall=T,
       render.continuous=my.render.cont2, render.categorical=my.render.cat)%>%
  cat(., file =paste0(i,".outcomes.pval.table1.html"))
df3<-df2[df2$inout_cd=="I",]
table1(~ length_stay| df3[[i]],
       data=df3, droplevels=F, render=rndr, render.strat=rndr.strat, overall=T,
       render.continuous=my.render.cont2, render.categorical=my.render.cat)%>%
  cat(., file =paste0(i,".outcomes.inpat.pval.table1.html"))
df3<-df2[df2$icu_days>0 &!is.na(df2$icu_days),]
table1(~ icu_days| df3[[i]],
       data=df3, droplevels=F, render=rndr, render.strat=rndr.strat, overall=T,
       render.continuous=my.render.cont2, render.categorical=my.render.cat)%>%
  cat(., file =paste0(i,".outcomes.icu.pval.table1.html"))

##Kaplan-Meyer Survival analysis overall
df3<-df[c("time_to_fu","death_yn")]
#convert categorical to numerical
must_convert<-(sapply(df3,is.character)|sapply(df3,is.factor)|sapply(df3,is.POSIXt) )# logical vector telling if a variable needs to be displayed as numeric
df3[must_convert]<-lapply(df3[must_convert],factor)
df3[must_convert]<-lapply(df3[must_convert],unclass)
df3[must_convert]<-lapply(df3[must_convert],as.numeric)
#model
fit<-survfit(Surv(df3$time_to_fu,df3$death_yn) ~1,df3)
fit
#plot
p<-ggsurvplot(fit,
              pval = T, conf.int = TRUE,
              surv.median.line = "hv", # Specify median survival
              xlab="Time (days after positive COVID-19 result)",
              xlim=c(0,50),
              break.time.by = 10,
              fontsize=2,color="black",
              ggtheme=theme_classic(),
)
png("overall.km.png", width = 1000, height = 800, res = 300)
print(p)
dev.off()


##Kaplan-Meyer Survival analysis per cohort
for(i in c("lab.abnormal.pH","agegroup","bm_failure_yn","white_yn","race.y","immunosup.any_yn","immunosup.drug_yn","lab.abnormal.Creatinine",
           "gender_MF","CKD_yn","immunosup_yn","pat_blood_type","lab.abnormal.Albumin","core.immunosup.drug_yn","lab.abnormal.Phosphate",
           "htn_yn","diabetes_yn","cvd_yn","malignancy_yn")){
  df2<-df[!is.na(df[[i]]),]
  df2<-df2[df2[[i]] %in% names(which(table(df[[i]])>4, arr.ind = FALSE, useNames = TRUE)),] #only keep groups with at least 5 pts
  df3<-df2[c("time_to_fu","death_yn",i)]
  #convert categorical to numerical
  must_convert<-(sapply(df3,is.character)|sapply(df3,is.factor)|sapply(df3,is.POSIXt) )# logical vector telling if a variable needs to be displayed as numeric
  df3[must_convert]<-lapply(df3[must_convert],factor)
  df3[must_convert]<-lapply(df3[must_convert],unclass)
  df3[must_convert]<-lapply(df3[must_convert],as.numeric)
  #model
  fit<-survfit(Surv(time_to_fu,death_yn) ~df3[[i]],df3)
  fit
  #plot
  png(paste0(i,".KM.png"), width = 1000, height = 1000, res = 300)
  print(ggsurvplot(fit,
               pval = T, conf.int = TRUE,
               linetype = "strata", # Change line type by groups
               surv.median.line = "hv", # Specify median survival
               ggtheme = theme_bw(), # Change ggplot2 theme
               xlab="Time (days after positive COVID-19 result)",
               xlim=c(0,50),
               break.time.by = 10,
               fontsize=5,
               ))
  dev.off()
  #cox regression
  res.cox <- coxph(Surv(time_to_fu,death_yn) ~df3[[i]],df3)
  print(summary(res.cox))
}
survfit(formula = Surv(df3$time_to_fu,df3$death_yn) ~ 1)


####Unbiased analyses
##Univariate dataframes
chisq.pval.df<-data.frame(test=vars.cat)
or.pval.df<-data.frame(test=vars.cat)
ttest.pval.df<-data.frame(test=vars.num)
wilcox.pval.df<-data.frame(test=vars.num)
glm.pval.df<-data.frame(test=vars.num)
aov.pval.df<-data.frame(test=vars.num)
kruskal.pval.df<-data.frame(test=vars.num)
for(j in outcome.cat){
  chisq.pval.df[[j]]<-NA
  or.pval.df[[j]]<-NA
  or.list<-list()
  for(i in vars.cat){
    #chi squared
    tbl<-table(df[[i]],df[[j]])
    if(nrow(tbl)>1 & ncol(tbl)>1 & sum(tbl>0)!=0){
      chisq.pval.df[[j]][chisq.pval.df$test==i]<-chisq.test(tbl)$p.value
    }
    #odds ratio
    if(nrow(tbl)==2 & ncol(tbl)==2 & sum(tbl>0)!=0){
      test<-try(epi.2by2(tbl,method="cohort.count"))
      if(!is(test,"try-error")){
        or.list[[i]]<-cbind(test$massoc$OR.strata.score,test$massoc$chisq.strata)
        or.pval.df[[j]][or.pval.df$test==i]<-test$massoc$chisq.strata$p.value
      }
    }
  }
  if(length(or.list)>1){
    or.df<-bind_rows(or.list, .id = "test")
    or.df$p.adj<-p.adjust(or.df$p.value,method="BH")
    or.df<-or.df[order(or.df$p.adj),]
    rownames(or.df)<-or.df$test
    or.df$test<-NULL
    assign(paste0(j,".or.df"),or.df)
    write.csv(or.df,paste0("./stats/",j,".or.df.csv"))
  }
  if(dim(table(df[[j]]))>2){
    aov.pval.df[[j]]<-NA
    kruskal.pval.df[[j]]<-NA
    for(i in vars.num){
      #one way ANOVA
      test<-try(aov(df[[i]] ~ df[[j]], data = df))
      if(!is(test,"try-error")){
        if(!is.null(summary(test)[[1]][["Pr(>F)"]])){
          aov.pval.df[[j]][aov.pval.df$test==i]<-summary(test)[[1]][["Pr(>F)"]][1]
        }
      }
      test<-try(kruskal.test(df[[i]] ~ df[[j]], data = df))
      if(!is(test,"try-error")){
        kruskal.pval.df[[j]][kruskal.pval.df$test==i]<-test$p.value
      }
    }
  }
  if(dim(table(df[[j]]))==2){
    ttest.pval.df[[j]]<-NA
    ttest.list<-list()
    wilcox.pval.df[[j]]<-NA
    df3<-df
    df3[[j]]<-as.factor(df3[[j]])
    glm.pval.df[[j]]<-NA
    for(i in vars.num){
      #ttest
      test<-try(t.test(df[[i]] ~ df[[j]], data = df))
      if(!is(test,"try-error")){
        ttest.list[[i]]<-c(p.val=test$p.value,test$estimate)
        ttest.pval.df[[j]][ttest.pval.df$test==i]<-test$p.value
      }
      #wilcox
      test<-try(wilcox.test(df[[i]] ~ df[[j]], data = df))
      if(!is(test,"try-error")){
        wilcox.pval.df[[j]][wilcox.pval.df$test==i]<-test$p.value
      }
      #logistic regression
      test<-try(glm(df3[[j]] ~ df3[[i]], data = df3,family=binomial))
      if(!is(test,"try-error")){
        if(!is.na(test$coefficients[2])){
          glm.pval.df[[j]][glm.pval.df$test==i]<-coef(summary(test))[2,4]
        }
      }
    }
    if(length(ttest.list)>1){
      ttest.df<-as.data.frame(map_dfr(lapply(ttest.list,unlist), ~as_tibble(t(.))))
      rownames(ttest.df)<-names(ttest.list)
      ttest.df$p.adj<-p.adjust(ttest.df$p.val, method = "BH")
      ttest.df<-ttest.df[order(ttest.df$p.adj),]
      assign(paste0(j,".ttest.df"),ttest.df)
      write.csv(ttest.df,paste0("./stats/",j,".ttest.csv"))
    }
  }
}


#plot univariate pvals
for(i in ls(pattern=".pval.df")){
  df2<-get(i)
  if("death_yn" %in% colnames(df2)){df2<-df2[order(df2$death_yn),]}
  if("ordscore" %in% colnames(df2)){df2<-df2[order(df2$ordscore),]}
  df2<-df2[!grepl("order_meds|death_yn|ordscore|disch|age_year",df2$test),] #exclude from variables
  df2<-df2[,!grepl("comp_sum|disch|causdth|shk|entry|vvecmo|corangio|hosp|ptc",colnames(df2))] #exclude from outcomes
  rownames(df2)<-df2$test
  df2$test<-NULL
  df2<-df2[rowSums(is.na(df2)) != ncol(df2), ] #delete rows that entirely empty
  df2<-df2[,colSums(is.na(df2)) != nrow(df2)] #delete rows that entirely empty
  df3<-as.data.frame(sapply(df2,p.adjust,method="BH"))#p.adj
  rownames(df3)<-rownames(df2)
  if("ordscore" %in% colnames(df2) & !"death_yn" %in% colnames(df2)){
    df3<-df3[df3$ordscore <0.05 & !is.na(df3$ordscore),]
  }
  if(!"ordscore" %in% colnames(df2) & "death_yn" %in% colnames(df2)){
    df3<-df3[df3$death_yn <0.05 & !is.na(df3$death_yn),]
  }
  if("ordscore" %in% colnames(df2) & "death_yn" %in% colnames(df2)){
    df3<-df3[(df3$death_yn <0.05 & !is.na(df3$death_yn))|(df3$ordscore <0.05 & !is.na(df3$ordscore)),]
  }
  write.csv(df3,paste0("./stats/",i,".csv"))
  df3<- -log(df3)
  df4<-do.call(data.frame,lapply(df3, function(x) replace(x, is.infinite(x),max(x[is.finite(x)], na.rm=T))))
  rownames(df4)<-rownames(df3)
  p<-pheatmap(df4,cluster_rows=F,cluster_cols=F)
  save_pheatmap_png(p, width=1500+ncol(df3)*100, height=1000+nrow(df3)*50,paste0(i,".heatmap.png"))
}


#cox regression (survival)
df2<-df[!is.na(df$death_yn),]
df3<-df2[c(vars.cat,"time_to_fu","death_yn")]
#convert categorical to numerical
must_convert<-(sapply(df3,is.character)|sapply(df3,is.factor)|sapply(df3,is.POSIXt) )# logical vector telling if a variable needs to be displayed as numeric
df3[must_convert]<-lapply(df3[must_convert],factor)
df3[must_convert]<-lapply(df3[must_convert],unclass)
df3[must_convert]<-lapply(df3[must_convert],as.numeric)
colnames(df3)<-gsub(" |,|-","_",colnames(df3))
colnames(df3)<-gsub("\\(|\\)|\\?|/","",colnames(df3))
#compute cox regression
univ_formulas <- sapply(colnames(df3),function(x) as.formula(paste("Surv(time_to_fu,death_yn) ~", x)))
univ_models <- lapply( univ_formulas, function(x){try(coxph(x, data = df3))})
univ_results <- lapply(univ_models,cox.stats)
univ_list<-univ_results[lengths(univ_results)==4]
res <- t(as.data.frame(univ_list, check.names = FALSE))
stats<-as.data.frame(res)
stats$p.value<-as.numeric(as.character(stats$p.value))
coxph.df<-stats
write.csv(stats[order(stats$p.value),],"coxph.csv")

##Spearman correlation plot
corr.var<-c(outcome.num,vars.num,outcome.cat,vars.cat)
corr.var<-corr.var[!str_detect(corr.var,pattern="_dt")]
corr.var.select<-vars.select
corr.var.select2<-c("dem_age","gender_MF","race.y","pat_blood_type","death_yn","ordscore","avg.bmi","icu_days","length_stay","enc.count","core.immunosup.drug_yn",
                    "mean.Absolute Lymphocyte Count","mean.C-Reactive Protein","mean.Creatinine","mean.Ferritin","mean.D-Dimer","mean.Creatine Kinase (CK)",
                    "mean.INR(PT)","mean.Lactate Dehydrogenase (LD)","mean.pH","mean.Platelet Count","mean.PT","mean.PTT","mean.Absolute Neutrophil Count",
                    "htn_yn","diabetes_yn","CKD_yn","immunosup_yn","rheumatologic_yn","autoimmune_yn","cancer_yn","obesity_yn",
                    "med.antiglycemic.drugs_yn","med.htn.drugs_yn","immunosup.any_yn","immunosup.drug_yn")
corr.var.select3<-c("death_yn","ordscore","dem_age","a1c","avg.bmi","agegroup","gender_MF","race.y","zipcode.y","dem_age","a1c",
                    "avg.bmi","dem_age","number_immuno_drugs","enc.count",
                    colnames(comorb)[-1],grep(paste(c("immun"),collapse="|"),colnames(df),value=T),
                    grep("immun",colnames(pmh),value=T),colnames(med))
corr.var.select4<-c("htn_yn","CKD_yn","diabetes_yn","obesity_yn","death_yn","dem_age","gender_MF","race.y","pat_blood_type",
                      "rheumatologic_yn","autoimmune_yn","malignancy_yn","immunosup_yn","mean.Absolute Lymphocyte Count","mean.C-Reactive Protein","mean.Creatinine","mean.Ferritin","mean.D-Dimer",
                      "mean.Creatine Kinase (CK)","mean.Lactate Dehydrogenase (LD)","mean.pH","mean.Platelet Count",
                      "mean.PT","mean.PTT","mean.Absolute Neutrophil Count","a1c","avg_RR","avg_HR","tmax","min_SBP","min_DBP","med.corticosteroid_yn","med.calcineurin_yn","med.antirheumatic_yn","med.nsaid.drugs_yn",
                      "med.onc.drugs_yn","med.antiglycemic.drugs_yn","med.asthma.drugs_yn","immunosup.drug_yn",
                      "med.biologics.drug_yn","med.bone.drugs_yn","med.opioid.drugs_yn","med.antihistamine.drugs_yn","med.gout.drugs_yn","med.htn.drugs_yn" ,
                    "comp_sum___rrt","poe_o2","poe_vent","enc.count","inout_cd","length_stay","ordscore","icu_days")
corr.var.select4<-corr.var.select4[!grepl("a1c|biologics|rrt",corr.var.select4)]
var.rename4<-c("Hypertension","CKD","Diabetes","Obesity","Mortality","Age","Gender","Race","ABO type","Rheumatologic Disease","Autoimmune Disease",
               "Cancer","Immunosuppression","Absolute Lymphocyte Count","C-Reactive Protein","Creatinine","Ferritin","D-Dimer","Creatine Kinase",
               "Lactate Dehydrogenase","pH","Platelet Count","PT","PTT","Absolute Neutrophil Count","RR","HR","Tmax","SBP","DBP","Corticosteroids",
               "Calcineurin Inhibitors","Antirheumatic Therapy","NSAIDs","Chemotherapy","Antiglycemic Therapy","Asthma Therapy",
               "Immunosuppressive Therapy","Osteoporosis Therapy","Opioids","Antihistamines","Gout Therapy","Antihypertensive Therapy","Supplemental O2",
               "Mechanical Ventilation","Encounters","Hospitalization","Length Admission","Ordinal Score","ICU days")
for(i in c("corr.var","corr.var.select","corr.var.select2","corr.var.select3","corr.var.select4")){
  var.list<-get(i)
  df2<-df[,unique(var.list)]
  if(i=="corr.var.select4"){
    for(j in 1:length(corr.var.select4)){
      names(df2)[names(df2) == corr.var.select4[j]] <- var.rename4[j]
    }
  }
  df2<-df2[rowSums(is.na(df2)) != ncol(df2), ] #delete rows that entirely empty
  df2<-df2[,colSums(is.na(df2)) != nrow(df2)] #delete cols that entirely empty
  df2<-df2[,colSums(!is.na(df2)) >50] #delete cols that have less than 50 pts
  df2<-df2[, sapply(df2, function(col) sum(!is.na(unique(col)))) > 1] #delete columns with only one type of factor
  #convert categorical to numerical
  must_convert<-(sapply(df2,is.character)|sapply(df2,is.factor)|sapply(df2,is.POSIXt) )# logical vector telling if a variable needs to be displayed as numeric
  df2[must_convert]<-lapply(df2[must_convert],factor)
  df2[must_convert]<-lapply(df2[must_convert],unclass)
  df2[must_convert]<-lapply(df2[must_convert],as.numeric)
  must_convert<-sapply(df2,is.logical) # pull logical columns
  df2[must_convert]<-lapply(df2[must_convert],as.numeric)
  #calculate correlation
  M<-cor(df2,use="pairwise.complete.obs",method="spearman")
  if(i=="corr.var"){
    w=15000
    h=15000}else{
      w=5000
      h=5000
    }
  #plot
  png(paste0(i,".corrplot.png"), width=w, height=h, res = 300)
  print(corrplot(M)) #na.labe=" "
  dev.off()
  test<-try(corrplot(M,order="hclust"))
  if(!is(test,"try-error")){
    res1 <- cor.mtest(df2, conf.level = .95)
    png(paste0(i,".corrplot.upper.png"), width=w, height=h, res = 300)
    print(corrplot(M, p.mat = res1$p, method = "color", type = "upper",
             sig.level = c(.001, .01, .05), pch.cex = .9,
             insig = "label_sig", pch.col = "white", order = "AOE"))
    dev.off()
    png(paste0(i,".corrplot.hclust.png"), width=w, height=h, res = 300)
    print(corrplot(M, p.mat = res1$p, sig.level = .05,insig="blank",na.label=".", order = "hclust", addrect = 4))
    dev.off()
  }
}

#cohort comparisons of Spearman, odds ratios, and FC plot
df2<-df[,unique(corr.var)]
df2<-df2[rowSums(is.na(df2)) != ncol(df2), ] #delete rows that entirely empty
df2<-df2[,colSums(is.na(df2)) != nrow(df2)] #delete cols that entirely empty
df2<-df2[,colSums(!is.na(df2)) >50] #delete cols that have less than 50 pts
df2<-df2[, sapply(df2, function(col) sum(!is.na(unique(col)))) > 1] #delete columns with only one type of factor
#convert categorical to numerical
must_convert<-(sapply(df2,is.character)|sapply(df2,is.factor)|sapply(df2,is.POSIXt) )# logical vector telling if a variable needs to be displayed as numeric
df2[must_convert]<-lapply(df2[must_convert],factor)
df2[must_convert]<-lapply(df2[must_convert],unclass)
df2[must_convert]<-lapply(df2[must_convert],as.numeric)
must_convert<-sapply(df2,is.logical) # pull logical columns
df2[must_convert]<-lapply(df2[must_convert],as.numeric)
for(k in c("ordscore","length_stay","death_yn","inout_cd","icu_days","enc.count")){
  for(i in c("immunosup_yn","core.immunosup.drug_yn","core.immunomod.drug_yn","core.biologic.drug_yn","ptc_immuno_drug",
             "ptc_immuno_drug_core","med.immunosup.drug_yn","immunosup.any_yn","immunosup.drug_yn")){
    cor.df<-list()
    or.df<-list()
    ttest.df<-list()
    df4<-df2[!is.na(df2[[i]]),]
    for(l in unique(df4[[i]])){
      df3<-df4[df4[[i]]==l,]
      cor.scores<-data.frame(metric=rep(NA,ncol(df3)),pval=rep(NA,ncol(df3)),rho=rep(NA,ncol(df3)),
                             z.score=rep(NA,ncol(df3)),stdev=rep(NA,ncol(df3)),std.error=rep(NA,ncol(df3)))
      or.list<-list()
      ttest.list<-list()
      for(j in 1:ncol(df3)){
        cor<-try(cor.test(df3[[j]],df3[[k]],method="spearman",use="pairwise.complete.obs"))
        if(!is(cor,"try-error")){
          cor.scores$pval[j]<-cor$p.value
          cor.scores$rho[j]<-cor$estimate
          cor.scores$metric[j]<-colnames(df3)[j]
          cor.scores$z.score[j]<-qnorm(1-cor$p.value)
          cor.scores$stdev[j]<-abs(cor$estimate)/abs(cor.scores$z.score[j])
          cor.scores$std.error[j]<-cor.scores$stdev[j]/sqrt(sum(!is.na(df3[[j]])))
        }
        tbl<-table(df3[[j]],df3[[k]])
        if(nrow(tbl)==2 & ncol(tbl)==2 & sum(tbl>0)!=0){
          or.test<-try(epi.2by2(tbl,method="cohort.count"))
          if(!is(or.test,"try-error")){
            or.result<-log(or.test$massoc$OR.strata.score)
            colnames(or.result)<-paste0("log.",colnames(or.result))
            or.result<-cbind(or.result,or.test$massoc$chisq.strata,or.test$massoc$OR.strata.score)
            or.result$n<-sum(!is.na(df3[[j]]))
            or.result$std.error<-sqrt((1/tbl[1,1])+(1/tbl[1,2])+(1/tbl[2,1])+(1/tbl[2,2]))
            or.list[[colnames(df3)[j]]]<-or.result
          }
        }
        if(dim(table(df3[[k]]))==2){
          test<-try(t.test(df3[[j]] ~ df3[[k]], data = df3))
          if(!is(test,"try-error")){
            ttest.result<-c(p.val=test$p.value,test$estimate,std.error=test$stderr,FC=as.numeric(test$estimate[2]/test$estimate[1]),
                            logFC=log(as.numeric(test$estimate[2]/test$estimate[1])),n=sum(!is.na(df3[[j]])))
            ttest.list[[colnames(df3)[j]]]<-ttest.result
          }
        }
      }
      # cor.scores[is.na(cor.scores)]<-0
      if(l==1){name<-"no"}
      if(l==2){name<-"yes"}
      cor.scores<-cor.scores[cor.scores$metric!=k&!is.na(cor.scores$rho),]
      colnames(cor.scores)<-paste0(name,".",colnames(cor.scores))
      colnames(cor.scores)[grepl("metric",colnames(cor.scores))]<-"metric"
      cor.df[[name]]<-cor.scores
      if(length(or.list)>1){
        or.scores<-bind_rows(or.list, .id = "or.test")
        colnames(or.scores)<-paste0(name,".",colnames(or.scores))
        colnames(or.scores)[grepl("or.test",colnames(or.scores))]<-"metric"
        or.df[[name]]<-or.scores
      }
      if(length(ttest.list)>1){
        ttest.scores<-as.data.frame(map_dfr(lapply(ttest.list,unlist), ~as_tibble(t(.))))
        ttest.scores$metric<-names(ttest.list)
        colnames(ttest.scores)<-paste0(name,".",colnames(ttest.scores))
        colnames(ttest.scores)[grepl("metric",colnames(ttest.scores))]<-"metric"
        ttest.df[[name]]<-ttest.scores
      }
    }
    cor.compare<-merge(cor.df[[1]], cor.df[[2]], by="metric")
    assign(paste0(k,".",i,".cor.compare"),cor.compare)
    if(length(or.df)>1){
      or.compare<-merge(or.df[[1]], or.df[[2]], by="metric")
      assign(paste0(k,".",i,".or.compare"),or.compare)
    }
    if(length(ttest.df)>1){
      ttest.compare<-merge(ttest.df[[1]], ttest.df[[2]], by="metric")
      assign(paste0(k,".",i,".ttest.compare"),ttest.compare)
    }
  }
  df3<-df2
  cor.scores<-data.frame(metric=rep(NA,ncol(df3)),pval=rep(NA,ncol(df3)),rho=rep(NA,ncol(df3)))
  for(j in 1:ncol(df3)){
    cor<-try(cor.test(df3[[j]],df3[[k]],method="spearman",use="pairwise.complete.obs"))
    if(!is(cor,"try-error")){
      cor.scores$pval[j]<-cor$p.value
      cor.scores$rho[j]<-cor$estimate
      cor.scores$metric[j]<-colnames(df3)[j]
    }else{
      cor.scores$pval[j]<-NA
      cor.scores$rho[j]<-NA
      cor.scores$metric[j]<-NA
    }
  }
  cor.scores<-cor.scores[cor.scores$metric!=k &!is.na(cor.scores$rho),]
  cor.scores$p.adj<-p.adjust(cor.scores$pval)
  write.csv(cor.scores[order(cor.scores$p.adj),],file=paste0(k,".cor.scores.csv"))
  assign(paste0(k,".cor.scores"),cor.scores)
  ranks <- cor.scores$rho
  names(ranks) <- cor.scores$metric
  png(paste0(k,".ranks.png"),width=5,height=4,units="in",res=200)
  print(barplot(sort(ranks, decreasing = T),cex.names=0.7,ylab=expression(paste("Correlation (",rho,")"))))
  dev.off()
}



##cohort vs cohort Volcano plots
coxph.df$`HR (95% CI for HR)`<-as.character(coxph.df$`HR (95% CI for HR)`)
coxph.df$HR<-as.numeric(sub(" .*","",coxph.df$`HR (95% CI for HR)`))
coxph.df$logFC<-log(coxph.df$HR)
coxph.df$p.val<-coxph.df$p.value
for(i in ls(pattern="ttest.df")){
  df2<-get(i)
  if(length(df2)>1){
    df2$logFC<-log(as.numeric(df2[[2]]/df2[[1]]))
    assign(i,df2)
  }
}
for(i in ls(pattern="\\.or.df")){
  df2<-get(i)
  if(length(df2)>1){
    df2$logFC<-log(df2$est)
    df2$p.val<-df2$p.value
    assign(i,df2)
  }
}
for(i in ls(pattern="\\.or.compare")){
  df2<-get(i)
  df2$logFC<-df2$yes.log.est-df2$no.log.est
  df2$std.dev<-sqrt((df2$yes.std.error)^2+(df2$no.std.error)^2)
  df2$z.score<-abs(df2$logFC)/df2$std.dev
  df2$p.val<-pnorm(df2$z.score, mean = 0, sd = 1, lower.tail = F)
  df2$p.adj<-p.adjust(df2$p.val,method="BH")
  # head(df2[order(df2$p.adj),])
  assign(i,df2)
}
for(i in ls(pattern=".cor.compare")){
  df2<-get(i)
  df2$logFC<-df2$yes.rho-df2$no.rho
  df2$std.dev<-sqrt((df2$yes.std.error)^2+(df2$no.std.error)^2)
  df2$z.score<-abs(df2$logFC)/df2$std.dev
  df2$p.val<-pnorm(df2$z.score, mean = 0, sd = 1, lower.tail = F)
  df2$p.adj<-p.adjust(df2$p.val,method="BH")
  # head(df2[order(df2$p.adj),])
  assign(i,df2)
}
for(i in ls(pattern=".ttest.compare")){
  df2<-get(i)
  df2$logFC<-df2$yes.logFC-df2$no.logFC
  df2$std.dev<-sqrt((df2$yes.std.error)^2+(df2$no.std.error)^2)
  df2$z.score<-abs(df2$logFC)/df2$std.dev
  df2$p.val<-pnorm(df2$z.score, mean = 0, sd = 1, lower.tail = F)
  df2$p.adj<-p.adjust(df2$p.val,method="BH")
  # head(df2[order(df2$p.adj),])
  assign(i,df2)
}

plots<-c("coxph.df",ls(pattern="or.df|ttest.df|.compare"))#.cor.scores|
plots<-plots[!grepl("comp_sum|ptc|shk|enc.count|drug",plots)]
for(i in plots){
  df2<-get(i)
  if(length(df2)>1&is(df2,"data.frame")){
    if("metric" %in% colnames(df2)){rownames(df2)<-df2$metric}
    df2<-df2[!(is.na(df2$logFC)|is.na(df2$p.val)),]
    df2<-df2[!(is.infinite(df2$logFC)),]
    df2<-df2[!grepl("order_meds|death_yn|ordscore|disch|minus|time_to_fu",rownames(df2)),]
    if(grepl("ttest",i)){
      df2<-df2[!grepl("comp_sum|pmh_|ptc_|abnormal|path|shk",rownames(df2)),]
    }
    df2$p.adj<-p.adjust(df2$p.val,method="BH")
    df2<-df2[df2$p.adj !=1,]
    x.label<-bquote(~Log~ (frac("Yes","No")))
    if(grepl("or.df",i)){x.label<-"Log(OR)"}
    if(i=="coxph.df"){x.label<-"Log(HR)"}
    if(grepl(".cor.scores",i)){x.label<-expression(paste("Log(",rho,")"))}
    if(grepl("inout_cd|patmanicudt",i)){
      df2$logFC<- -df2$logFC
    } #if(i %in% c("inout_cd.or.df","patmanicudt.or.df")){
    if(nrow(df2)>10){
      write.csv(df2[order(df2$p.adj),],file=paste0(i,".volcano.csv"))
      metric.list<-head(rownames(df2[order(df2$p.adj),]),n=10)
      EnhancedVolcano(df2,lab = rownames(df2),
                      x = 'logFC',y = 'p.adj',title = i,col=c("black","black","black","red3"),
                      selectLab=metric.list,xlab=x.label,
                      pCutoff = 0.05,FCcutoff = 0.2,pointSize = 0.5,labSize = 3,axisLabSize=10,colAlpha = 1, #transparencyxlim = c(-1.5, 1.5),
                      legendPosition="none",#drawConnectors = TRUE,widthConnectors = 0.2,colConnectors = 'grey30',
                      subtitle="", caption="",border="full",cutoffLineWidth=0,
                      gridlines.major=F,gridlines.minor=F,titleLabSize=10
      )
      ggsave2(paste0(i,".volcano.png"),width=4, height=4,device="png")
    }
  }
}

# violin plots for continuous vs death
df2<-df[!is.na(df$death_yn),]
df2<-df2[,colSums(!is.na(df2)) >100] #delete cols that have less than 600 pts
# boxplot.list<-list()
# violin.list<-list()
violin.list2<-list()
violin.list3<-list()
for(i in grep("mean",colnames(df2),value=T)){
  violin.list2[[i]]<-ggplot(data=df2,aes(x=death_yn,y=.data[[i]]))+ #,fill=death_yn
    geom_violin(trim=F)+labs(x="Death",y=i)+theme_classic()+
    stat_compare_means(method = "wilcox.test",label="p.format",label.x=1.5,label.y=1.3*max(df2[[i]],na.rm=T))+ 
    geom_jitter(shape=16, position=position_jitter(0.2),color="black",size=0.5) + 
    theme(legend.position = "none")
  violin.list3[[i]]<-ggplot(data=df2,aes(x=death_yn,y=.data[[i]],fill=ptc_immuno_drug))+
    geom_violin(trim=F)+labs(x="Death",y=i)+theme_classic()+
    stat_compare_means(method = "wilcox.test",label="p.format",label.x=1.5,label.y=1.3*max(df2[[i]],na.rm=T))+ 
    geom_point(position=position_jitterdodge(jitter.width=0.1,dodge.width=0.9),color="black",size=0.5)
}
for(i in c("violin.list2","violin.list3")){
  list<-get(i)
  wrap_plots(list,guides="collect")&theme(legend.position = "bottom")
  ggsave2(paste0(i,".mortality.png"),width=15, height=15,device="png")
}

# 
#forest plot of Odds ratios
var<-c("gender_MF","htn_yn","CKD_yn","cvd_yn","diabetes_yn","cad_yn","asthma_yn","copd_yn",
       "obesity_yn","rheumatologic_yn","autoimmune_yn","malignancy_yn","immunosup_yn")#,"lab.abnormal.pH"
lab<-c("Male","Hypertension","Chronic Kidney Disease","Cardiovascular Disease","Diabetes","Coronary Artery Disease",
       "Asthma","COPD","Obesity","Rheumatologic Disease",
       "Autoimmune Disease","Cancer","Immunosuppression")
for(i in c("death_yn.or.df","inout_cd.or.df","patmanicudt.or.df","poe_vent.or.df","poe_o2.or.df")){
  df2<-get(i)
  df2<-df2[var,]
  fplot<-cbind(metric=rownames(df2),df2[1:3])
  if(i %in% c("inout_cd.or.df","patmanicudt.or.df")){
    fplot[2:4]<-1/fplot[2:4]
  }
  fplot$metric<-factor(fplot$metric,levels=var,labels=lab)
  write.csv(fplot,paste0(i,".fplot.csv"))
  ggplot(fplot,aes(x=metric,y=est,ymin=lower,ymax=upper))+
    geom_pointrange()+geom_hline(yintercept=1, lty=2) + coord_flip() +
    labs(x="",y="OR (95% CI)")+
    theme(panel.border = element_rect(colour = "black", size=1,fill=NA),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  ggsave2(paste0(i,".fplot.png"),width=3, height=3,device="png")
}
for(i in ls(pattern=".ptc_immuno_drug.or.compare")){
  df2<-get(i)
  df2<-df2[df2$metric %in% var,]
  df3<-as.data.frame(rbind(cbind(metric=df2$metric,est=df2$no.est,lower=df2$no.lower,upper=df2$no.upper,immmun_drug="no"),
                           cbind(metric=df2$metric,est=df2$yes.est,lower=df2$yes.lower,upper=df2$yes.upper,immmun_drug="yes")),
                     stringsAsFactors = F)
  df3[c("est","lower","upper")]<-lapply(df3[c("est","lower","upper")],as.numeric)
  ggplot(df3,aes(x=metric,y=est,ymin=lower,ymax=upper,color=immmun_drug))+
    geom_pointrange(position=position_dodge(width=0.5))+
    geom_hline(yintercept=1, lty=2) + coord_flip() +
    labs(x="",y="OR (95% CI)")+
    theme(panel.border = element_rect(colour = "black", size=1,fill=NA),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  ggsave2(paste0(i,".fplot.png"),width=5, height=7,device="png")
}

#scatter plot for continuous variables w/ correlation coefficient
df.names<-df
colnames(df.names)<-make.names(colnames(df.names))
scatter.vars<-c("max.Absolute Neutrophil Count","max.C-Reactive Protein","max.D-Dimer","mean.Albumin",
                "max.pH","max.Ferritin","ptc_ox_sat","max.Phosphate","min.Hemoglobin","min.Albumin",
                "dem_age","number_immuno_drugs")
scatter.vars<-c("max_RR","max.Absolute Neutrophil Count","max.C-Reactive Protein","min.Albumin",
                "max.Neutrophils","min.Lymphocytes","max.Phosphate","max_HR",
                "tmax")#"max.D-Dimer",,"min.Hemoglobin","min_SBP"
for(j in c("ordscore","length_stay")){ #,"enc.count","icu_days"
  plot.list<-list()
  for(k in make.names(scatter.vars)){
    df4<-df.names[!is.na(df.names[[k]])&!is.na(df.names[[j]]),]
    if(is.character(df4[[k]])){df4[[k]]<-as.numeric(as.character(df4[[k]]))}
    ylabel<-max(df4[[k]])+IQR(df4[[k]])
    ylimit<-max(df4[[k]])+IQR(df4[[k]])*1.5
    p<-ggscatter(df4, x = j, y = k, add = "reg.line", conf.int = TRUE, cor.method = "spearman",size=0.5,
                 add.params = list(color = "blue",fill = "grey20"),cor.coef=T,ylim=c(NA,ylimit),
                 cor.coeff.args = list(method = "spearman", label.sep = "\n",label.y=ylabel)) #label.x = 3, 
    if(j=="length_stay"|j=="icu_days"){
      p<-p+scale_x_continuous(trans="log10")
    }
    plot.list[[k]]<-p
  }
  wrap_plots(plot.list)
  ggsave2(paste0(j,".scatter.png"),width=10, height=10,device="png")
}


##PCA/UMAP
#select variables to input
outcome.num<-unique(c("tempdurhome","ptc_ox_sat","ptc_gcs","ptc_supp_liters",
                      "vs_o2vol","vs_fio2","oxygen_saturation_spo2_v2","ventdays","icu_days",
                      "pd_calc_icu_days","icu.count","ordscore","length_stay","time_to_fu","enc.count",
                      grep(paste(c("floor","min","max","avg"),collapse="|"),colnames(df),value=T)))
vars.num<-c("dem_age","a1c","avg.bmi","age_year","number_immuno_drugs",grep(paste(c("min","max","avg"),collapse="|"),colnames(df),value=T))
outcome.cat<-unique(c("outcome","ordscore","hospvent","patmanicudt","vvecmo","corangio","dschstat","ever_hosp","oc_newstatus","ptc_entry","ptc_present_intubate_yn",
                      "ptc_ox_supp_yn","ptc_intubate_yn","causdth","loc","death","vap","ecmo","rrt","were_infiltrates_present",
                      "disch_disp_full","poe_o2","poe_vent","inout_cd","ibax_disch_disp desc","death_yn",
                      grep(paste(c("shk","discharge","oxygenation","comp_","cxr_pattern","ct_chest","lung_ultrasound"),collapse="|"),colnames(df),value=T)))
vars.cat<-unique(c("ptc_immuno_drug","ptc_immuno_drug_core","ckd_stages","dem_sex","dem_gender","race.x","ethnicity.x","sub_etoh","sub_drug","sub_cig","sub_ends","sub_mj",
                   "ptc_insteroid_dur","ptc_orsteroid_dur","prone","if_yes_steroid_route_of_ad_v2","zipcode.x","agegroup","med.immunomod.drugs_yn","Tobacco use",
                   colnames(comorb)[-1],"pat_blood_type","pat_rh","gender","race.y","ethnicity.y","ethnic_cd desc","hispanic_ind","zipcode.y","gender_MF","dem_preg",
                   grep(paste(c("ptc_tx","pmh","ptc_sx__","path","blood_products","medications__","antiviral__","lab.abnormal","_yn"),collapse="|"),colnames(df),value=T)))
pca.var<-c("deID_key",outcome.num,vars.num,outcome.cat,vars.cat)
pca.var<-subset(pca.var, !(pca.var %in% c("ordscore","length_stay","icu_days","ventdays","pd_calc_icu_days","time_to_fu","outcome",
                                       "hospvent","patmanicudt","vvecmo","corangio","dschstat","oc_newstatus","causdth","death","death_yn",
                                       "disch_disp_full","ibax_disch_disp desc","enc.count","time_to_fu","outcome","poe_vent",
                                       grep("discharge|order_meds|ptc_tx|minus|path|abnormal",colnames(df),value=T)))) 
pca.var<-subset(pca.var, !(pca.var %in% c(grep("path|ptc_tx|comp_sum|pmh|core|ptc_sx",colnames(df),value=T)))) 
pca.var<-c("deID_key","dem_age","a1c","avg.bmi",colnames(lab.val)[-1],"race.y","zipcode.y","gender_MF",
           grep(paste(c("_yn"),collapse="|"),colnames(df),value=T))
pca.var<-subset(pca.var, !(pca.var %in% c(grep("order_meds|death|immunomod|core",colnames(df),value=T),
                                          "immunosup.any_yn"))) #exclude listed variables, |med.
# pca.var<-corr.var
df2<-df[,unique(pca.var)]#unique(vars.num,vars.cat)
# df2<-df
df2<-df2[,colSums(is.na(df2)) != nrow(df2)] #delete cols that entirely empty
df2<-df2[, sapply(df2, function(col) sum(!is.na(unique(col)))) > 1] #delete columns with only one type of factor
df2<-df2[,colSums(!is.na(df2)) >500] #delete cols that have less than 600 pts (use 500 if excluding meds)
# df2<-df2[rowSums(!is.na(df2))>100, ] #delete rows that have less than 400 variables
df2<-df2[rowSums(is.na(df2))==0, ] #delete cols that have missing data
df2<-df2[, colSums(is.na(df2))==0] #delete cols that have missing data
df2<-df2[, sapply(df2, function(col) sum(!is.na(unique(col)))) > 1] #delete columns with only one type of factor
dim(df2)
#convert categorical to numerical
df3<-df2 #store df for factor reference
df4<-left_join(df2[,"deID_key"],df,by="deID_key") #store df for removed var reference
df2$deID_key<-NULL
must_convert<-(sapply(df2,is.character)|sapply(df2,is.factor)|sapply(df2,is.POSIXt) )# logical vector telling if a variable needs to be displayed as numeric
df2[must_convert]<-lapply(df2[must_convert],factor)
df2[must_convert]<-lapply(df2[must_convert],unclass)
df2[must_convert]<-lapply(df2[must_convert],as.numeric)
must_convert<-sapply(df2,is.logical) # pull logical columns
df2[must_convert]<-lapply(df2[must_convert],as.numeric)
###run PCA
res.pca <- PCA(df2,  graph = FALSE)
# get_eig(res.pca)
p1<-fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 50))
# var <- get_pca_var(res.pca)
p2<-fviz_contrib(res.pca, choice = "var", axes = 1, top = 10,xtickslab.rt=90,y.text.angle=90)
p2
ggsave("PC1.png",width=2.5, height=3,device="png") #p1,
p3<-fviz_contrib(res.pca, choice = "var", axes = 2, top = 10,xtickslab.rt=90,y.text.angle=90)
p3
ggsave("PC2.png",width=2.5, height=3.5,device="png") #p1,
grid.arrange(p1,p2,p3,ncol=3)
ggsave("PCs.png",arrangeGrob(p2,p3,ncol=2),width=6, height=4,device="png") #p1,
# ind <- get_pca_ind(res.pca)
fviz_pca_var(res.pca, col.var = "black",repel=TRUE,select.var=list(contrib=5))
ggsave("PCA.vars.png",width=6, height=6,device="png")
# fviz_pca_biplot(res.pca, repel = TRUE,habillage = as.factor(df2$death_yn),select.var=list(contrib=10),label="var")
fviz_pca_biplot(res.pca, select.var=list(contrib=5),label="var")
fviz_pca_biplot(res.pca, select.var=list(
  name=c("immunosup_yn","max.Neutrophils","mean.Asparate Aminotransferase (AST)","max.Ferritin")),label="var")
ggsave("PCA.biplot.png",width=4, height=4,device="png")
#PCA plot
plot.metrics<-c("death_yn","ptc_immuno_drug","gender_MF","race.y","length_stay","dem_age","htn_yn","enc.count",
                "number_immuno_drugs","icu_days","max.Absolute Neutrophil Count",
                "max_RR","max.Creatinine","max_HR","min.Albumin","inout_cd","poe_vent",
                "ordscore","immunosup_yn","core.immunosup.drug_yn")
plot.metrics<-c("death_yn","gender_MF","race.y","length_stay","dem_age",
                "ordscore")
# plot.metrics.present<-plot.metrics[plot.metrics %in% colnames(df2)]
pca.plot<-cbind(as.data.frame(res.pca$ind$coord)[1:2],df4[,c(plot.metrics)])
factor.metrics<-c("death_yn","gender_MF","race.y")#"ptc_immuno_drug","inout_cd","poe_vent"
pca.plot[factor.metrics]<-lapply(pca.plot[factor.metrics],factor)
plot.list<-list()
for(i in plot.metrics){
  p<-ggplot(data=pca.plot,aes(x=Dim.1,y=Dim.2,color=.data[[i]]))+
    geom_point(size=0.5)+theme(panel.border = element_rect(colour = "black", size=1,fill=NA),
                               panel.background = element_blank(),
                               panel.grid.major = element_blank(), 
                               panel.grid.minor = element_blank())
  if(! i %in% factor.metrics){p<-p+scale_color_gradient(low="blue", high="red")}
  plot.list[[i]]<-p
}
wrap_plots(plot.list)
ggsave2("pca.png",width=11, height=5,device="png")

###run UMAP
custom.config = umap.defaults
# custom.config$random_state = 123
umap<-umap(df2, config=custom.config)
# umap.plot<-cbind(as.data.frame(umap$layout),df3[,c(plot.metrics.present)])
umap.plot<-cbind(as.data.frame(umap$layout),df4[,c(plot.metrics)])
umap.plot[factor.metrics]<-lapply(umap.plot[factor.metrics],factor)
plot.list<-list()
for(i in plot.metrics){
  plot.list[[i]]<-ggplot(data=umap.plot,aes(x=V1,y=V2,color=.data[[i]]))+
    geom_point(size=1)+theme(panel.border = element_rect(colour = "black", size=1,fill=NA),
                       panel.background = element_blank(),
                       panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank())+
    labs(x="UMAP_1",y="UMAP_2")
}
wrap_plots(plot.list)
ggsave2("umap.png",width=12, height=5,device="png")

##Machine learning for mortality
#choose variables to use
outcome.num<-unique(c("tempdurhome","ptc_ox_sat","ptc_gcs","ptc_supp_liters",
                      "vs_o2vol","vs_fio2","oxygen_saturation_spo2_v2","ventdays","icu_days",
                      "pd_calc_icu_days","icu.count","ordscore","length_stay","time_to_fu","enc.count",
                      grep(paste(c("floor","min","max","avg"),collapse="|"),colnames(df),value=T)))
vars.num<-c("dem_age","a1c","avg.bmi","age_year","number_immuno_drugs",grep(paste(c("min","max","avg"),collapse="|"),colnames(df),value=T))
outcome.cat<-unique(c("outcome","ordscore","hospvent","patmanicudt","vvecmo","corangio","dschstat","ever_hosp","oc_newstatus","ptc_entry","ptc_present_intubate_yn",
                      "ptc_ox_supp_yn","ptc_intubate_yn","causdth","loc","death","vap","ecmo","rrt","were_infiltrates_present",
                      "disch_disp_full","poe_o2","poe_vent","inout_cd","ibax_disch_disp desc","death_yn",
                      grep(paste(c("shk","discharge","oxygenation","comp_","cxr_pattern","ct_chest","lung_ultrasound"),collapse="|"),colnames(df),value=T)))
vars.cat<-unique(c("ptc_immuno_drug","ptc_immuno_drug_core","ckd_stages","dem_sex","dem_gender","race.x","ethnicity.x","sub_etoh","sub_drug","sub_cig","sub_ends","sub_mj",
                   "ptc_insteroid_dur","ptc_orsteroid_dur","prone","if_yes_steroid_route_of_ad_v2","zipcode.x","agegroup","med.immunomod.drugs_yn","Tobacco use",
                   colnames(comorb)[-1],"pat_blood_type","pat_rh","gender","race.y","ethnicity.y","ethnic_cd desc","hispanic_ind","zipcode.y","gender_MF","dem_preg",
                   grep(paste(c("ptc_tx","pmh","ptc_sx__","path","blood_products","medications__","antiviral__","lab.abnormal","_yn"),collapse="|"),colnames(df),value=T)))
ML.var<-c(outcome.num,vars.num,outcome.cat,vars.cat)
# ML.var<-pca.var
ML.var<-subset(ML.var, !(ML.var %in% c("ordscore","length_stay","icu_days","ventdays","pd_calc_icu_days","time_to_fu","outcome",
                                           "hospvent","patmanicudt","vvecmo","corangio","dschstat","oc_newstatus","causdth","death",
                                           "disch_disp_full","ibax_disch_disp desc","enc.count","time_to_fu","outcome","poe_vent",
                                           grep("discharge|order_meds|ptc_tx|minus|path",colnames(df),value=T)))) #exclude outcome variables (except one outcome modeling)
ML.var<-c("deID_key","dem_age","a1c","avg.bmi",colnames(lab.val)[-1],"race.y","zipcode.y","gender_MF",
           grep(paste(c("_yn"),collapse="|"),colnames(df),value=T))
ML.var<-subset(ML.var, !(ML.var %in% c(grep("order_meds|immunomod|core",colnames(df),value=T),
                                       "immunosup.any_yn"))) #exclude listed variables, |med.
df2<-df[,unique(ML.var)]#unique(vars.num,vars.cat)
# df2<-df
df2<-df2[,colSums(is.na(df2)) != nrow(df2)] #delete cols that entirely empty
df2<-df2[, sapply(df2, function(col) sum(!is.na(unique(col)))) > 1] #delete columns with only one type of factor
df2<-df2[,colSums(!is.na(df2)) >500] #delete cols that have less than 600 pts (use 500 if excluding meds)
df2<-df2[rowSums(is.na(df2))==0, ] #delete rows that have missing data
df2<-df2[, colSums(is.na(df2))==0] #delete cols that have missing data
df2<-df2[, sapply(df2, function(col) sum(!is.na(unique(col)))) > 1] #delete columns with only one type of factor
dim(df2)
#convert categorical to numerical
must_convert<-(sapply(df2,is.character)|sapply(df2,is.factor)|sapply(df2,is.POSIXt) )# logical vector telling if a variable needs to be displayed as numeric
df2[must_convert]<-lapply(df2[must_convert],factor)
colnames(df2)<-make.names(colnames(df2))
#split into training and test data sets
index<-createDataPartition(df2$death_yn,p=0.75,list=F)
df.training<-df2[index,]
df.test<-df2[-index,]
set.seed(123)
ctrl=trainControl(method="repeatedcv",repeats=10,classProbs=TRUE,summaryFunction=twoClassSummary,savePredictions=T)
lift_results<-data.frame(obs=df.test[["death_yn"]])
mat.list<-list()
performance.list<-list()
ROC.list<-list()
ROC.list2<-list()
#train ML
i="gbm"
model<- train(as.factor(death_yn) ~.,data=df.training,#df.training[, !names(df.training) %in% c("death_yn")], df.training[["death_yn"]], #
              method=i,preProcess=c("center","scale"),metric="ROC",trControl=ctrl) 
assign(paste0(i,".model"),model)
#visualize performance
model
p<-try(ggplot(model))
if(!is(p,"try-error")){
  ggsave2(paste0(i,".training.png"),plot=p,width=4, height=4,device="png")
}
#variable importance
var.imp<-try(varImp(model,scale=F))
if(!is(var.imp,"try-error")){
  var.imp
  write.csv(var.imp$importance,paste0(i,".var.imp.csv"))
  png(paste0(i,".var.imp.png"), width = 1200, height = 1500, res = 300)
  print(plot(var.imp,top=20))
  dev.off()
}
roc_imp<-try(filterVarImp(x=df.training[, !names(df.training) %in% c("death_yn")],y=df.training[["death_yn"]]))
if(!is(roc_imp,"try-error")){
  roc_imp<-roc_imp[order(-roc_imp$no),]
  write.csv(roc_imp,paste0(i,".var.roc.imp.csv"))
  roc_imp2<-head(roc_imp,n=20)
  ggplot(roc_imp2,aes(x=reorder(rownames(roc_imp2),yes), y=yes)) +
    geom_point( color="blue", size=4, alpha=0.6)+
    geom_segment( aes(x=rownames(roc_imp2), xend=rownames(roc_imp2), y=0, yend=yes)) +#,color='skyblue'
    labs(x="",y="Importance")+ theme(panel.border = element_rect(colour = "black", size=1,fill=NA), panel.background = element_blank())+
    coord_flip() 
  ggsave2(paste0(i,".roc.imp.png"),width=4, height=7,device="png")
}
# Predict the labels of the test set
pred<-predict.train(object=model,df.test[, !names(df.test) %in% c("death_yn")], type="raw")
predictions.prob<-predict.train(object=model,df.test[, !names(df.test) %in% c("death_yn")], type="prob")
mat.list[[i]]<-confusionMatrix(pred,df.test[["death_yn"]])
#performance metrics
performance.list[[i]]<-postResample(pred,df.test[["death_yn"]])
test_set<-data.frame(obs=df.test[["death_yn"]],pred=pred)
test_set<-cbind(test_set,predictions.prob)
ROC.list[[i]]<-twoClassSummary(test_set, lev = levels(df.test[["death_yn"]]))
ROC.list2[[i]]<-prSummary(test_set, lev = levels(df.test[["death_yn"]]))
#ROC
res<-evalm(model)
res$roc
ggsave2(paste0(i,".ROC.png"),width=10, height=10,device="png")
#lift curve
lift_results[[i]]<-predictions.prob[,"yes"]


