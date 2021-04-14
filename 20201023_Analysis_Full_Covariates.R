library(OlinkAnalyze)
library(corrplot)
library(heatmap3)
library(Cairo)


#load image rather than running lines 12:292
rm(list=ls())
load("20201023_Analysis_Full_Covariates.RData")


####### setup #######
#set working directory to file/project location
##Be sure to load workspace generated at the end of 20201023_Analysis.R

#### To do

# 1. Correlation heatmaps      - DONE 17JULY2020
# 2. Sample heatmaps           - DONE 17JULY2020
# 3. LMER's, control for age   - DONE 17JULY2020


####### 24JUL20 - Add full covariates to all LMM models
# 1. Correlation heatmaps      -
# 2. Sample heatmaps           -
# 3. LMER's, control for age   -


lmm.covariates <- c("heart.disease","Diabetes","HTN","HLD","Pulmonary","kidney.disease","Immuno.comprimised","sex","ethnicity","age.decile")

load("../../NPX data/20201023_clean_merged_data")
source("../heatmap_functions.R")

npx.info <- npx.info %>% 
  filter(COVID=="Positive") %>% 
  filter(Timepoint != "DE")


# LMM models with all covariates ---------------------------------------------------------------
# 
# severity.lmer.covars <- olink_lmer(df=npx.info,
#                                 variable=c("Timepoint","WHO.max"),
#                                 outcome="NPX",
#                                 covariates = lmm.covariates,
#                                 random="Patient")
# 
# sig.interaction.covars <- severity.lmer.covars %>%
#   filter(Threshold=="Significant" & term=="Timepoint:WHO.max") %>%
#   pull(OlinkID)
# 
# ##Pvalue distributions
# ggplot(severity.lmer.covars,aes(x=Adjusted_pval))+
#   geom_histogram(color="white",bins=20)+
#   set_plot_theme()+
#   facet_wrap(~term)+
#   labs(title="BH Adjusted p-values")
# ggplot(severity.lmer.covars,aes(x=p.value))+
#   geom_histogram(color="white",bins=20)+
#   set_plot_theme()+
#   facet_wrap(~term)+
#   labs(title="Unadjusted p-values")
# 
# 
# ##tally significances
# severity.lmer.covars %>% select(OlinkID,term,Threshold) %>%
#   pivot_wider(names_from="term",values_from="Threshold") %>%
#   group_by(Timepoint,WHO.max,`Timepoint:WHO.max`) %>%
#   tally() 
# 
# %>%
# #   left_join(severity.lmer %>% select(OlinkID,term,Threshold) %>%
# #               pivot_wider(names_from="term",values_from="Threshold") %>%
# #               group_by(Timepoint,WHO.max,`Timepoint:WHO.max`) %>%
# #               tally(),
# #             by=c("Timepoint","WHO.max","Timepoint:WHO.max"),suffix=c(".Covariates",".None"))
# 
# 
# 
# 
##all sig interaction post hocs
# all.ph <- olink_lmer_posthoc(df=npx.info,    ##This is available from the loaded environment
#                                     variable=c("Timepoint","WHO.max"),
#                                     olinkid_list = unique(npx.info$OlinkID),
#                                     effect=c("Timepoint","WHO.max"),
#                                     outcome="NPX",
#                                     # covariates = lmm.covariates,
#                                     random="Patient")
# 
# all.ph.covars <- olink_lmer_posthoc(df=npx.info,
#                                  variable=c("Timepoint","WHO.max"),
#                                  olinkid_list = unique(npx.info$OlinkID),
#                                  effect=c("Timepoint","WHO.max"),
#                                  outcome="NPX",
#                                  covariates = lmm.covariates,
#                                  random="Patient")
# 
# severity.lmer.covars.ph <- list()
# 
# severity.lmer.covars.ph[["Significant Interaction"]] <- all.ph.covars %>%
#   filter(OlinkID %in% sig.interaction.covars) %>%
#   arrange(Adjusted_pval)%>%
#   separate(contrast,c("day.left","who.left","day.right","who.right"),",| - ",remove=F) %>%
#   mutate(Day_Comparison = paste(day.left,day.right,sep="-"),
#          WHO_Comparison= paste(who.left,who.right,sep="-")) %>%
#   select(-c(day.left,day.right,who.left,who.right))
# 
# ##Due to missingness structure, model has to be adjusted to Timepoint + WHO.max. Indicator included for significant interaction
# severity.lmer.covars.ph[["Significant Timepoint"]] <- olink_lmer_posthoc(df=npx.info,
#                                                                       variable=c("Timepoint"),
#                                                                       olinkid_list = severity.lmer.covars %>%
#                                                                         filter(!OlinkID %in% sig.interaction.covars) %>%
#                                                                         filter(Threshold=="Significant" & term=="Timepoint") %>%
#                                                                         pull(OlinkID),
#                                                                       effect=c("Timepoint"),
#                                                                       covariates = c(lmm.covariates,"WHO.max"),
#                                                                       outcome="NPX",
#                                                                       random="Patient") %>%
#   mutate(significant_interaction = ifelse(OlinkID %in% sig.interaction.covars,"Yes","No")) %>%
#   arrange(Adjusted_pval) %>%
#   filter(significant_interaction=="No")
# 
# severity.lmer.covars.ph[["Significant WHO.max"]] <- olink_lmer_posthoc(df=npx.info,
#                                                                     variable=c("WHO.max"),
#                                                                     olinkid_list = severity.lmer.covars %>%
#                                                                       filter(Threshold=="Significant" & term=="WHO.max") %>%
#                                                                       filter(!OlinkID %in% sig.interaction.covars) %>%
#                                                                       pull(OlinkID),
#                                                                     effect=c("WHO.max"),
#                                                                     covariates=c(lmm.covariates,"Timepoint"),
#                                                                     outcome="NPX",
#                                                                     random="Patient") %>%
#   mutate(significant_interaction = ifelse(OlinkID %in% sig.interaction.covars,"Yes","No")) %>%
#   arrange(Adjusted_pval) %>%
#   filter(significant_interaction=="No")
# 
# severity.lmer.covars.ph[["Global F-test"]] <- severity.lmer.covars
# 
# #Write results to excel
# olinkr::output_dfs_to_excel(severity.lmer.covars.ph,
#                             filename = "lmm/20201023_LMER_all_covariates_results.xlsx")
# 
# 
# 
# 
# 
# # volcano plots all covariates-----------------------------------------------------------
# 
# volcano.data.covars <- all.ph.covars %>%
#   separate(contrast,c("day.left","who.left","day.right","who.right"),",| - ",remove=F) %>%
#   filter(day.left==day.right) %>%
#   mutate(significant = ifelse(Threshold=="Significant" &
#                                 OlinkID %in% (severity.lmer.covars %>%
#                                                 filter(Threshold=="Significant" & term %in% c("WHO.max","Timepoint:WHO.max")) %>%
#                                                 pull(OlinkID)),
#                               "Significant","Non-significant")) %>%
#   mutate(Day_Comparison = paste(day.left,day.right,sep="-"),
#          WHO_Comparison= paste(who.left,who.right,sep="-")) %>%
#   filter(!is.na(estimate))
# 
# 
# volc.plots.covars <- volcano.data.covars %>%
#   mutate(Adjusted_pval = ifelse(Adjusted_pval==0,min(Adjusted_pval[Adjusted_pval!=0]),Adjusted_pval),
#          Day_Comparison=factor(Day_Comparison,levels=c("DE-DE","D0-D0","D3-D3","D7-D7")),
#          WHO_Comparison=paste0("WHO: ",WHO_Comparison)) %>%
#   ggplot(aes(x=estimate,y=-log10(Adjusted_pval))) +
#   geom_point(aes(color=significant),alpha=0.4)+
#   geom_hline(yintercept=-log10(.05),linetype="dashed")+
#   geom_vline(xintercept=0,linetype="dashed")+
#   facet_grid(WHO_Comparison~Day_Comparison)+
#   theme_bw()+
#   theme(strip.background = element_rect(fill = NA),
#         legend.position="bottom")+
#   scale_color_manual(values=c("lightgray","#00559EFF"))+
#   labs(x="NPX Difference",y="-log10(p-value)",color="")+
#   guides(colour = guide_legend(override.aes = list(alpha=1)))
# 
# ggsave(filename="lmm/volcano_plots_all_covariates.pdf",
#        plot=volc.plots.covars,
#        height=20,
#        width=10,
#        units = "in",
#        device = cairo_pdf)
# 
# 
# 
# 
# 
# 
# 
# # post hoc with and without covariates scatter -----------------------------------
# 
# all.means <- olink_lmer_posthoc(df=npx.info,
#                                 variable=c("Timepoint","WHO.max"),
#                                 olinkid_list = unique(npx.info$OlinkID),
#                                 effect=c("Timepoint","WHO.max"),
#                                 outcome="NPX",
#                                 random="Patient",
#                                 mean_return = T)
# all.means.covars <- olink_lmer_posthoc(df=npx.info,
#                                     variable=c("Timepoint","WHO.max"),
#                                     olinkid_list = unique(npx.info$OlinkID),
#                                     effect=c("Timepoint","WHO.max"),
#                                     outcome="NPX",
#                                     covariates = lmm.covariates,
#                                     random="Patient",
#                                     mean_return = T)
# 
# means.comp <- all.means %>% select(OlinkID,Timepoint,WHO.max,emmean) %>% filter(!is.na(emmean)) %>%
#   left_join(all.means.covars %>% select(OlinkID,Timepoint,WHO.max,emmean) %>% filter(!is.na(emmean)),
#             by=c("OlinkID","Timepoint","WHO.max"),
#             suffix=c("WHO and Covariates","WHO Only")) %>%
#   mutate(WHO.max=paste("WHO:",WHO.max)) %>%
#   ggplot(aes(y=`emmeanWHO and Covariates`,x=`emmeanWHO Only`))+
#   geom_point(alpha=.2,color="#00559EFF")+
#   set_plot_theme()+
#   facet_wrap(Timepoint~WHO.max,scales='free')+
#   geom_abline(intercept=0,slope=1)+
#   labs(x="Mean Estimate, WHO Only",
#        y="Mean Estimate, WHO and Covariates")
# 
# all.means %>% select(OlinkID,Timepoint,WHO.max,emmean) %>% filter(!is.na(emmean)) %>%
#   left_join(all.means.covars %>% select(OlinkID,Timepoint,WHO.max,emmean) %>% filter(!is.na(emmean)),
#             by=c("OlinkID","Timepoint","WHO.max"),
#             suffix=c("WHO and Covariates","WHO Only")) %>%
#   mutate(WHO.max=paste("WHO:",WHO.max)) %>%
#   mutate(diff=`emmeanWHO Only`-`emmeanWHO and Covariates`) %>%
#   ggplot(aes(x=diff))+
#   geom_density()+
#   facet_wrap(Timepoint~WHO.max)+
#   set_plot_theme()+
#   geom_vline(xintercept=0)
# 
# ggsave(filename="lmm/lmer_means_scatter_compare_all_covars.pdf",
#        plot=means.comp,
#        height=15,
#        width=10,
#        units = "in",
#        device = cairo_pdf)
# 
# #PH
# 
# ph.comp <- all.ph %>% select(OlinkID,contrast,estimate) %>% filter(!is.na(estimate)) %>%
#   left_join(all.ph.covars %>% select(OlinkID,contrast,estimate) %>% filter(!is.na(estimate)),
#             by=c("OlinkID","contrast"),
#             suffix=c("WHO and Age","WHO Only")) %>%
#   separate(contrast,c("day.left","who.left","day.right","who.right"),",| - ",remove=F) %>%
#   filter(day.left==day.right) %>%
#   mutate(Day_Comparison = paste(day.left,day.right,sep="-"),
#          WHO_Comparison= paste(who.left,who.right,sep="-")) %>%
#   mutate(Day_Comparison=factor(Day_Comparison,levels=c("DE-DE","D0-D0","D3-D3","D7-D7")),
#          WHO_Comparison=paste0("WHO: ",WHO_Comparison)) %>%
#   ggplot(aes(y=`estimateWHO and Age`,x=`estimateWHO Only`))+
#   geom_point(alpha=.2,color="#00559EFF")+
#   set_plot_theme()+
#   facet_wrap(WHO_Comparison~Day_Comparison,scales="free")+
#   geom_vline(xintercept=0)+geom_hline(yintercept=0)+
#   geom_abline(intercept=0,slope=1)+
#   labs(x="Contrast Estimate, WHO Only",
#        y="Contrast Estimate, WHO and Age")
# 
# 
# all.ph %>% select(OlinkID,contrast,estimate) %>% filter(!is.na(estimate)) %>%
#   left_join(all.ph.covars %>% select(OlinkID,contrast,estimate) %>% filter(!is.na(estimate)),
#             by=c("OlinkID","contrast"),
#             suffix=c("WHO and Age","WHO Only")) %>%
#   separate(contrast,c("day.left","who.left","day.right","who.right"),",| - ",remove=F) %>%
#   filter(day.left==day.right) %>%
#   mutate(Day_Comparison = paste(day.left,day.right,sep="-"),
#          WHO_Comparison= paste(who.left,who.right,sep="-")) %>%
#   mutate(Day_Comparison=factor(Day_Comparison,levels=c("DE-DE","D0-D0","D3-D3","D7-D7")),
#          WHO_Comparison=paste0("WHO: ",WHO_Comparison)) %>%
#   mutate(diff=`estimateWHO Only`-`estimateWHO and Age`) %>%
#   ggplot(aes(x=diff))+
#   geom_density()+
#   facet_wrap(Day_Comparison~WHO_Comparison)+
#   set_plot_theme()+
#   geom_vline(xintercept=0)
# 
# 
# 
# ggsave(filename="lmm/lmer_contrasts_scatter_compare_all_covars.pdf",
#        plot=ph.comp,
#        height=15,
#        width=10,
#        units = "in",
#        device = cairo_pdf)
# 
# 
# 
# 
# save.image("20201023_Analysis_Full_Covariates.RData")
# 
# assay Correlation matrices by significance group ------------------------------

cut.seq <- seq(.5,1.9,.01)
pref.clusters <- 8

cluster.list <- list()

oid.list <- list(interaction=severity.lmer.covars$OlinkID[severity.lmer.covars$Threshold=="Significant" & severity.lmer.covars$term=="Timepoint:WHO.max"],
                 timepoint=severity.lmer.covars$OlinkID[severity.lmer.covars$Threshold=="Significant" & severity.lmer.covars$term=="Timepoint"],
                 who=severity.lmer.covars$OlinkID[severity.lmer.covars$Threshold=="Significant" & severity.lmer.covars$term=="WHO.max"])

for(i in names(oid.list)){
  
  cluster.list[[i]] <- makeAssayCorrHeatmap(data=npx.info,
                                            OIDs=oid.list[[i]],
                                            filename = paste0("assay correlations/",i,"_proteins_corr_heatmap_allCovars.pdf"),
                                            out.width=20,out.height=20,
                                            cut.search=cut.seq,
                                            cut.n=pref.clusters,
                                            method="pearson",
                                            return.cluster=T)
  
}

olinkr::output_dfs_to_excel(cluster.list,
                            filename="assay correlations/clusters_Proteins_allCovars.xlsx",
                            df_names = names(cluster.list))


# assay correlation heatmaps by timepoint ---------------------------------



oid.list.time <- list(D0 = unique(c(severity.lmer.covars.ph[["Significant Interaction"]] %>% 
                                 filter(Day_Comparison == "D0-D0",Threshold=="Significant") %>% pull(OlinkID),
                               severity.lmer.covars.ph[["Significant WHO.max"]] %>% 
                                 filter(Threshold=="Significant") %>% pull(OlinkID))),
                 D3 = unique(c(severity.lmer.covars.ph[["Significant Interaction"]] %>% 
                                 filter(Day_Comparison == "D3-D3",Threshold=="Significant") %>% pull(OlinkID),
                               severity.lmer.covars.ph[["Significant WHO.max"]] %>% 
                                 filter(Threshold=="Significant") %>% pull(OlinkID))),
                 D7 = unique(c(severity.lmer.covars.ph[["Significant Interaction"]] %>% 
                                 filter(Day_Comparison == "D7-D7",Threshold=="Significant") %>% pull(OlinkID),
                               severity.lmer.covars.ph[["Significant WHO.max"]] %>% 
                                 filter(Threshold=="Significant") %>% pull(OlinkID)))
)


cluster.list <- list()

for(i in names(oid.list.time)){
  
  cluster.list[[i]] <- makeAssayCorrHeatmap(data=npx.info %>% filter(Timepoint==i),
                                            OIDs=oid.list.time[[i]],
                                            filename = paste0("assay correlations/by_time/",i,"_proteins_corr_heatmap.pdf"),
                                            out.width=20,out.height=20,
                                            cut.search=cut.seq,
                                            cut.n=pref.clusters,
                                            method="pearson",
                                            return.cluster=T)
  
  
  for(j in names(oid.list)){
    cluster.list[[paste(i,j)]] <- makeAssayCorrHeatmap(data=npx.info %>% filter(Timepoint==i),
                                              OIDs=oid.list[[j]],
                                              filename = paste0("assay correlations/by_time/",j,"_proteins_corr_heatmap_",i,".pdf"),
                                              out.width=20,out.height=20,
                                              cut.search=cut.seq,
                                              cut.n=pref.clusters,
                                              method="pearson",
                                              return.cluster=T)
    
  }
  
  
  
}

olinkr::output_dfs_to_excel(cluster.list,
                            filename="assay correlations/by_time/clusters_Proteins_allCovars.xlsx",
                            df_names = names(cluster.list))




# Sample Heatmaps ---------------------------------------

cut.seq <- seq(.5,1.95,.01)
pref.clusters <- 2

top50.label <- severity.lmer.covars.ph[["Global F-test"]] %>% filter(term %in% c("WHO.max","Timepoint:WHO.max"),Threshold=="Significant") %>% 
  mutate(OlinkID=paste0(Assay,"_",OlinkID)) %>% 
  group_by(OlinkID) %>% summarize(statistic=max(statistic)) %>% arrange(desc(statistic)) %>% ungroup() %>% 
  slice(1:50) %>% pull(OlinkID)


cluster.list <- list()


for(useN in  c(25,50,100,200)){
  oid.list[[paste0("top",useN)]] <- severity.lmer.covars.ph[["Global F-test"]] %>% filter(term %in% c("WHO.max","Timepoint:WHO.max")) %>% 
    group_by(OlinkID) %>% summarize(statistic=max(statistic)) %>% arrange(desc(statistic)) %>% ungroup() %>% 
    slice(1:useN) %>% pull(OlinkID)
}

for(i in names(oid.list)){
  
  
  
  cluster.list[[i]] <- makeHeatmap(data=npx.info,color.vars=c("WHO.max","WHO2","Timepoint","age.group"),color.labels=c("WHO Max","WHO2","Time","Age Group"),
                                   OIDs=oid.list[[i]],
                                   filename = paste0("sample heatmaps/sample_heatmap_",i,"_proteins.pdf"),
                                   out.width=20,out.height=20,
                                   cut.search=cut.seq,
                                   cut.n=pref.clusters,
                                   return.cluster=T)
                                   # row.labels = top50.label)
  
}

olinkr::output_dfs_to_excel(cluster.list,
                            filename="sample heatmaps/sample_heatmap_clusters_Proteins.xlsx",
                            df_names = names(cluster.list))


# Sample Heatmaps by timepoint ---------------------------------------



cluster.list <- list()

for(i in names(oid.list.time)){
  
  cluster.list[[i]] <- makeHeatmap(data=npx.info %>% filter(Timepoint==i),
                                   color.vars=c("WHO.max","WHO2","age.group"),color.labels=c("WHO Max","WHO2","Age Group"),
                                   OIDs=oid.list.time[[i]],
                                   filename = paste0("sample heatmaps/by_time_sample_heatmap_",i,"_proteins.pdf"),
                                   out.width=20,out.height=20,
                                   cut.search=cut.seq,
                                   cut.n=pref.clusters,
                                   return.cluster=T,
                                   row.labels = top50.label)
  
  for(j in names(oid.list)){
    
    cluster.list[[paste(i,j)]] <- makeHeatmap(data=npx.info %>% filter(Timepoint==i),
                                              color.vars=c("WHO.max","WHO2","age.group"),color.labels=c("WHO Max","WHO2","Age Group"),
                                              OIDs=oid.list[[j]],
                                              filename = paste0("sample heatmaps/by_time/sample_heatmap_",i,"_proteins_",j,".pdf"),
                                              out.width=20,out.height=20,
                                              cut.search=cut.seq,
                                              cut.n=pref.clusters,
                                              return.cluster=T,
                                              row.labels = top50.label)
    
  }
  
}

olinkr::output_dfs_to_excel(cluster.list,
                            filename="sample heatmaps/by_time/sample_heatmap_clusters_Proteins.xlsx",
                            df_names = names(cluster.list))







# Use Residuals -----------------------------------------------------------

single_lmm_return_residuals <- function(formula,data){
  mod <- lmer(as.formula(formula),data=data)
  data$lm.residual <- residuals(mod)
  return(data)
}

npx.info.resids <- npx.info %>% 
  group_by(OlinkID) %>% 
  mutate(NPX=ifelse(is.na(NPX),median(NPX,na.rm=T),NPX)) %>% 
  group_modify(~single_lmm_return_residuals(formula=paste0("NPX~",paste(lmm.covariates,collapse="+"),"+(1|Patient)"),data=.x)) %>% 
  ungroup()

npx.info.resids <- npx.info.resids %>% 
  mutate(NPX=lm.residual)


# Use Residuals assay Correlation matrices by significance group ------------------------------

cut.seq <- seq(.5,1.9,.01)
pref.clusters <- 8

cluster.list <- list()

oid.list <- list(interaction=severity.lmer.covars$OlinkID[severity.lmer.covars$Threshold=="Significant" & severity.lmer.covars$term=="Timepoint:WHO.max"],
                 timepoint=severity.lmer.covars$OlinkID[severity.lmer.covars$Threshold=="Significant" & severity.lmer.covars$term=="Timepoint"],
                 who=severity.lmer.covars$OlinkID[severity.lmer.covars$Threshold=="Significant" & severity.lmer.covars$term=="WHO.max"])

for(i in names(oid.list)){
  
  cluster.list[[i]] <- makeAssayCorrHeatmap(data=npx.info.resids,
                                            OIDs=oid.list[[i]],
                                            filename = paste0("assay correlations/on_residuals/",i,"_proteins_corr_heatmap_allCovars_useResids.pdf"),
                                            out.width=20,out.height=20,
                                            cut.search=cut.seq,
                                            cut.n=pref.clusters,
                                            method="pearson",
                                            return.cluster=T)
  
}

olinkr::output_dfs_to_excel(cluster.list,
                            filename="assay correlations/on_residuals/clusters_Proteins_allCovars_useResids.xlsx",
                            df_names = names(cluster.list))


# Use Residuals assay correlation heatmaps by timepoint ---------------------------------



oid.list.time <- list(D0 = unique(c(severity.lmer.covars.ph[["Significant Interaction"]] %>% 
                                      filter(Day_Comparison == "D0-D0",Threshold=="Significant") %>% pull(OlinkID),
                                    severity.lmer.covars.ph[["Significant WHO.max"]] %>% 
                                      filter(Threshold=="Significant") %>% pull(OlinkID))),
                      D3 = unique(c(severity.lmer.covars.ph[["Significant Interaction"]] %>% 
                                      filter(Day_Comparison == "D3-D3",Threshold=="Significant") %>% pull(OlinkID),
                                    severity.lmer.covars.ph[["Significant WHO.max"]] %>% 
                                      filter(Threshold=="Significant") %>% pull(OlinkID))),
                      D7 = unique(c(severity.lmer.covars.ph[["Significant Interaction"]] %>% 
                                      filter(Day_Comparison == "D7-D7",Threshold=="Significant") %>% pull(OlinkID),
                                    severity.lmer.covars.ph[["Significant WHO.max"]] %>% 
                                      filter(Threshold=="Significant") %>% pull(OlinkID)))
)


cluster.list <- list()

for(i in names(oid.list.time)){
  
  cluster.list[[i]] <- makeAssayCorrHeatmap(data=npx.info.resids %>% filter(Timepoint==i),
                                            OIDs=oid.list.time[[i]],
                                            filename = paste0("assay correlations/on_residuals/by_time/",i,"_proteins_corr_heatmap_useResids.pdf"),
                                            out.width=20,out.height=20,
                                            cut.search=cut.seq,
                                            cut.n=pref.clusters,
                                            method="pearson",
                                            return.cluster=T)
  
  
  for(j in names(oid.list)){
    cluster.list[[paste(i,j)]] <- makeAssayCorrHeatmap(data=npx.info.resids %>% filter(Timepoint==i),
                                                       OIDs=oid.list[[j]],
                                                       filename = paste0("assay correlations/on_residuals/by_time/",j,"_proteins_corr_heatmap_",i,"_useResids.pdf"),
                                                       out.width=20,out.height=20,
                                                       cut.search=cut.seq,
                                                       cut.n=pref.clusters,
                                                       method="pearson",
                                                       return.cluster=T)
    
  }
  
  
  
}

olinkr::output_dfs_to_excel(cluster.list,
                            filename="assay correlations/on_residuals/by_time/clusters_Proteins_allCovars_useResids.xlsx",
                            df_names = names(cluster.list))




# Use Residuals Sample Heatmaps ---------------------------------------

cut.seq <- seq(.5,1.95,.01)
pref.clusters <- 2
top50.label <- severity.lmer.covars.ph[["Global F-test"]] %>% filter(term %in% c("WHO.max","Timepoint:WHO.max"),Threshold=="Significant") %>% 
  mutate(OlinkID=paste0(Assay,"_",OlinkID)) %>% 
  group_by(OlinkID) %>% summarize(statistic=max(statistic)) %>% arrange(desc(statistic)) %>% ungroup() %>% 
  slice(1:50) %>% pull(OlinkID)


cluster.list <- list()


for(useN in  c(25,50,100,200)){
  oid.list[[paste0("top",useN)]] <- severity.lmer.covars.ph[["Global F-test"]] %>% filter(term %in% c("WHO.max","Timepoint:WHO.max")) %>% 
    group_by(OlinkID) %>% summarize(statistic=max(statistic)) %>% arrange(desc(statistic)) %>% ungroup() %>% 
    slice(1:useN) %>% pull(OlinkID)
}

for(i in names(oid.list)){
  
  cluster.list[[i]] <- makeHeatmap(data=npx.info.resids,color.vars=c("WHO.max","WHO2","Timepoint","age.group"),color.labels=c("WHO Max","WHO2","Time","Age Group"),
                                   OIDs=oid.list[[i]],
                                   filename = paste0("sample heatmaps/on_residuals/sample_heatmap_",i,"_proteins_useResids.pdf"),
                                   out.width=20,out.height=20,
                                   cut.search=cut.seq,
                                   cut.n=pref.clusters,
                                   return.cluster=T,
                                   row.labels = top50.label)
  
}

olinkr::output_dfs_to_excel(cluster.list,
                            filename="sample heatmaps/on_residuals/sample_heatmap_clusters_Proteins_useResids.xlsx",
                            df_names = names(cluster.list))


# Use Residuals Sample Heatmaps by timepoint ---------------------------------------



cluster.list <- list()

for(i in names(oid.list.time)){
  
  cluster.list[[i]] <- makeHeatmap(data=npx.info.resids %>% filter(Timepoint==i),
                                   color.vars=c("WHO.max","WHO2","age.group"),color.labels=c("WHO Max","WHO2","Age Group"),
                                   OIDs=oid.list.time[[i]],
                                   filename = paste0("sample heatmaps/on_residuals/by_time/sample_heatmap_",i,"_proteins_useResids.pdf"),
                                   out.width=20,out.height=20,
                                   cut.search=cut.seq,
                                   cut.n=pref.clusters,
                                   return.cluster=T,
                                   row.labels = top50.label)
  
  for(j in names(oid.list)){
    
    cluster.list[[paste(i,j)]] <- makeHeatmap(data=npx.info.resids %>% filter(Timepoint==i),
                                              color.vars=c("WHO.max","WHO2","age.group"),color.labels=c("WHO Max","WHO2","Age Group"),
                                              OIDs=oid.list[[j]],
                                              filename = paste0("sample heatmaps/on_residuals/by_time/sample_heatmap_",i,"_proteins_",j,"_useResids.pdf"),
                                              out.width=20,out.height=20,
                                              cut.search=cut.seq,
                                              cut.n=pref.clusters,
                                              return.cluster=T,
                                              row.labels = top50.label)
    
  }
  
}

olinkr::output_dfs_to_excel(cluster.list,
                            filename="sample heatmaps/on_residuals/by_time/sample_heatmap_clusters_Proteins_useResids.xlsx",
                            df_names = names(cluster.list))


