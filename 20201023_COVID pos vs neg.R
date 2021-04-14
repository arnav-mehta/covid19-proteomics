library(OlinkAnalyze)
library(corrplot)
library(heatmap3)
library(Cairo)
# library(factoextra)
#set working directory to file/project location


#### To do

# 1. 

# load data ---------------------------------------------------------------
heatmap.color.scale <- rev(colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                                              "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                              "#4393C3", "#2166AC", "#053061"))(1024))

lmm.covariates <- c("heart.disease","Diabetes","HTN","HLD","Pulmonary","kidney.disease","Immuno.comprimised","sex","ethnicity","age.decile")

getCutHeight <- function(mod,search.seq, n.clusters){
  
  cuth <-  sapply(search.seq,function(x){
    length(unique(cutree(mod,h=x)))
  })
  cuth <- cuth[cuth>=n.clusters]
  hld <- which.min(abs(n.clusters-cuth))
  cat("Cutting into ",cuth[hld]," clusters\n")
  cuth <- search.seq[hld]
  
  return(cuth)
}



single_lm <- function(formula,data){
  mod <- lm(as.formula(formula),data=data)
  out <- as.data.frame(Anova(mod,type=3))
  out$term <- rownames(out)
  rownames(out) <- NULL
  return(out[!grepl("(Intercept)|Residuals",out$term),c(5,1:4)])
}
single_lm_posthoc <- function(formula,data){
  mod <- lm(as.formula(formula),data=data)
  out <- emmeans(mod,pairwise~COVID,data=data)$contrasts
  return(as.data.frame(out))
}


load("../../NPX data/20201023_clean_merged_data")



# COVID+ vs COVID- , all covars --------------------------------------------------------

##No COVID- samples at DE and D7. Too few at D3. Remove from analysis. Remove replicate samples, use first
npx.info.covid <- npx.info %>% filter(Timepoint %in% c("D0")) %>% filter(Rep==1)

#######check all covars

covid.covars <- npx.info.covid %>%
  group_by(OlinkID,Assay) %>%
  group_modify(~single_lm(as.formula(paste0("NPX~COVID + ",paste(lmm.covariates,collapse="+"))),
                          data=.x)) %>%
  ungroup() %>% rename(p.value=`Pr(>F)`) %>%
  group_by(term) %>%
  mutate(Adjusted_pval=p.adjust(p.value,method="BH")) %>%
  ungroup()

##pvalue histograms; kidney disease has alot of significances
covid.covars %>%  ##red=BH pvalue; blue=nominal pvalue
  ggplot(aes(x=p.value))+
  geom_histogram(color="black",alpha=.4,fill="blue",bins=40)+
  geom_histogram(aes(x=Adjusted_pval),color="black",alpha=.4,fill="red",bins=40)+
  facet_wrap(~term)+
  set_plot_theme()


covid.lm.results <- npx.info.covid %>%
  group_by(OlinkID,Assay) %>%
  group_modify(~single_lm_posthoc(as.formula(paste0("NPX~COVID + ",paste(lmm.covariates,collapse="+"))),
                          data=.x)) %>% ungroup()

covid.lm.results <- covid.lm.results %>%
  mutate(Adjusted_pval=p.adjust(p.value,method="BH"),
         Threshold=ifelse(Adjusted_pval<.05,"Significant","Non-significant")) %>%
  arrange(desc(Threshold),desc(abs(estimate)))
# 
# 
# olinkr::output_dfs_to_excel(list(covid.covars,covid.lm.results),
#                             df_names=c("ANOVA","posthoc"),
#                             filename="lm/20201023_COVID_lm_results.xlsx")

##volcano plot
volc.plot <- covid.lm.results %>% 
  mutate(add.label = rank(-abs(estimate))<=40) %>% {
    ggplot(data=.,aes(x=estimate,y=-log10(Adjusted_pval),color=Threshold))+
      geom_point(alpha=.4)+
      geom_label_repel(aes(label=Assay),data=subset(.,add.label),show.legend=F,size=3)+
      labs(x="NPX Difference\nCOVID POS - NEG",y="-log10(pvalue)",color="")+
      set_plot_theme()+
      scale_color_manual(values=c("lightgray","#00559EFF"))
    }
ggsave(filename="lm/volcano_plots_all_covariates.pdf",
       plot=volc.plot+theme(legend.position = "bottom"),
       height=8,
       width=10,
       units = "in",
       device = cairo_pdf)

##### for ida
# look.oids <- npx.info %>% select(Assay,OlinkID) %>% unique() %>% 
#   filter(grepl("ACE2|IL6|MCP|KRT19|CXCL10|AREG|IFN|TLR3",Assay))
# 
# covid.lm.results %>% filter(OlinkID %in% look.oids$OlinkID)
# 
# 
# boxplots.look <- olink_boxplot(npx.info.covid,"COVID",number_of_proteins_per_plot = 15,
#                                olinkid_list = look.oids$OlinkID)




# assay correlations ------------------------------------------------------
cut.seq <- seq(.5,1.8,.01)
pref.clusters <- 8

cluster.list <- list()


##all proteins
tmp <- npx.info.covid %>% 
  mutate(Assay_OlinkID=paste0(Assay,"_",OlinkID)) %>% 
  select(Assay_OlinkID,SampleID,NPX) %>% 
  group_by(Assay_OlinkID) %>% 
  mutate(NPX=(NPX-mean(NPX,na.rm=T))/sd(NPX,na.rm=T)) %>% 
  mutate(NPX=ifelse(is.na(NPX)|!is.finite(NPX),0,NPX)) %>% ##18 values are missing. replace with mean (scaled data, mean=0)
  ungroup() %>% 
  pivot_wider(names_from="Assay_OlinkID",values_from="NPX") %>% 
  select(-SampleID) %>% 
  cor(method="pearson")
min.col <- which.min(min(tmp) > seq(-1,1,length.out = 1024))
max.col <- which.min(max(tmp) > seq(-1,1,length.out = 1024))

write.csv(tmp,
          file = "assay correlations/all_proteins_pearson_cor.csv")



##### All proteins

CairoPDF("assay correlations/all_proteins_corr_heatmap.pdf",
         width=20,
         height=20)


tmp.mod <- hclust(as.dist(1-as.matrix(tmp)))
cuth <- getCutHeight(tmp.mod,search.seq = cut.seq,n.clusters=pref.clusters)
out.clusters <- cutree(tmp.mod,h=cuth)
cluster.list[["All Proteins"]] <- data.frame(Assay = names(out.clusters),
                                Cluster=as.character(as.numeric(out.clusters))) %>% arrange(Cluster)

heatmap3(tmp,
         scale="none",
         method="complete",
         revC=T,
         ColSideCut = cuth,
         distfun = function(x) as.dist(1 - (x)),
         col = heatmap.color.scale[min.col:max.col],
         legendfun =  function() {showLegend(legend = "",col="white")})
plotfunctions::gradientLegend(valRange=c(-1,1),
                              color = heatmap.color.scale,
                              pos=c(0,.85,.025,.99), # xleft, ybottom, xright, ytop
                              fit.margin = F,cex=2)
dev.off()




##Significant proteins
tmp <- npx.info.covid %>% 
  filter(OlinkID %in% covid.lm.results$OlinkID[covid.lm.results$Threshold=="Significant"]) %>% 
  mutate(Assay_OlinkID=paste0(Assay,"_",OlinkID)) %>% 
  select(Assay_OlinkID,SampleID,NPX) %>% 
  group_by(Assay_OlinkID) %>% 
  mutate(NPX=(NPX-mean(NPX,na.rm=T))/sd(NPX,na.rm=T)) %>% 
  mutate(NPX=ifelse(is.na(NPX)|!is.finite(NPX),0,NPX)) %>% ##18 values are missing. replace with mean (scaled data, mean=0)
  ungroup() %>% 
  pivot_wider(names_from="Assay_OlinkID",values_from="NPX") %>% 
  select(-SampleID) %>% 
  cor(method="pearson")
min.col <- which.min(min(tmp) > seq(-1,1,length.out = 1024))
max.col <- which.min(max(tmp) > seq(-1,1,length.out = 1024))

write.csv(tmp,
          file = "assay correlations/significant_proteins_pearson_cor.csv")



##### All significant proteins

CairoPDF("assay correlations/significant_proteins_corr_heatmap.pdf",
         width=20,
         height=20)


tmp.mod <- hclust(as.dist(1-as.matrix(tmp)))
cuth <- getCutHeight(tmp.mod,search.seq = cut.seq,n.clusters=pref.clusters)
out.clusters <- cutree(tmp.mod,h=cuth)
cluster.list[["Significant Proteins"]] <- data.frame(Assay = names(out.clusters),
                                             Cluster=as.character(as.numeric(out.clusters))) %>% arrange(Cluster)

heatmap3(tmp,
         scale="none",
         method="complete",
         revC=T,
         ColSideCut = cuth,
         distfun = function(x) as.dist(1 - (x)),
         col = heatmap.color.scale[min.col:max.col],
         legendfun =  function() {showLegend(legend = "",col="white")})
plotfunctions::gradientLegend(valRange=c(-1,1),
                              color = heatmap.color.scale,
                              pos=c(0,.85,.025,.99), # xleft, ybottom, xright, ytop
                              fit.margin = F,cex=2)
dev.off()


olinkr::output_dfs_to_excel(cluster.list,
                            filename="assay correlations/clusters_Proteins.xlsx",
                            df_names = names(cluster.list))



#### Upregulated in COVID
tmp <- npx.info.covid %>%
  filter(OlinkID %in% covid.lm.results$OlinkID[covid.lm.results$Threshold=="Significant" &
                                                 covid.lm.results$estimate>0]) %>%
  mutate(Assay_OlinkID=paste0(Assay,"_",OlinkID)) %>%
  select(Assay_OlinkID,SampleID,NPX) %>%
  group_by(Assay_OlinkID) %>%
  mutate(NPX=(NPX-mean(NPX,na.rm=T))/sd(NPX,na.rm=T)) %>%
  mutate(NPX=ifelse(is.na(NPX)|!is.finite(NPX),0,NPX)) %>% ##18 values are missing. replace with mean (scaled data, mean=0)
  ungroup() %>%
  pivot_wider(names_from="Assay_OlinkID",values_from="NPX") %>%
  select(-SampleID) %>%
  cor(method="pearson")
min.col <- which.min(min(tmp) > seq(-1,1,length.out = 1024))
max.col <- which.min(max(tmp) > seq(-1,1,length.out = 1024))
tmp.mod <- hclust(as.dist(1-as.matrix(tmp)))
cuth <- getCutHeight(tmp.mod,search.seq = cut.seq,n.clusters=8)

CairoPDF("assay correlations/upregulated_proteins_corr_heatmap.pdf",
         width=20,
         height=20)
heatmap3(tmp,
         scale="none",
         method="complete",
         revC=T,
         ColSideCut = cuth,
         distfun = function(x) as.dist(1 - (x)),
         col = heatmap.color.scale[min.col:max.col],
         legendfun =  function() {showLegend(legend = "",col="white")})
plotfunctions::gradientLegend(valRange=c(-1,1),
                              color = heatmap.color.scale,
                              pos=c(0,.85,.025,.99), # xleft, ybottom, xright, ytop
                              fit.margin = F,cex=2)
dev.off()
#### Down regulated in COVID
tmp <- npx.info.covid %>%
  filter(OlinkID %in% covid.lm.results$OlinkID[covid.lm.results$Threshold=="Significant" &
                                                 covid.lm.results$estimate<0]) %>%
  mutate(Assay_OlinkID=paste0(Assay,"_",OlinkID)) %>%
  select(Assay_OlinkID,SampleID,NPX) %>%
  group_by(Assay_OlinkID) %>%
  mutate(NPX=(NPX-mean(NPX,na.rm=T))/sd(NPX,na.rm=T)) %>%
  mutate(NPX=ifelse(is.na(NPX)|!is.finite(NPX),0,NPX)) %>% ##18 values are missing. replace with mean (scaled data, mean=0)
  ungroup() %>%
  pivot_wider(names_from="Assay_OlinkID",values_from="NPX") %>%
  select(-SampleID) %>%
  cor(method="pearson")
min.col <- which.min(min(tmp) > seq(-1,1,length.out = 1024))
max.col <- which.min(max(tmp) > seq(-1,1,length.out = 1024))
tmp.mod <- hclust(as.dist(1-as.matrix(tmp)))
cuth <- getCutHeight(tmp.mod,search.seq = cut.seq,n.clusters=8)

CairoPDF("assay correlations/downregulated_proteins_corr_heatmap.pdf",
         width=20,
         height=20)
heatmap3(tmp,
         scale="none",
         method="complete",
         revC=T,
         ColSideCut = cuth,
         distfun = function(x) as.dist(1 - (x)),
         col = heatmap.color.scale[min.col:max.col],
         legendfun =  function() {showLegend(legend = "",col="white")})
plotfunctions::gradientLegend(valRange=c(-1,1),
                              color = heatmap.color.scale,
                              pos=c(0,.85,.025,.99), # xleft, ybottom, xright, ytop
                              fit.margin = F,cex=2)
dev.off()


# sample heatmaps ----------------------------------------------------------------

cut.seq <- seq(.5,1.95,.01)
pref.clusters <- 2

### All proteins

tmp.wide <- npx.info.covid %>% 
  # filter(OlinkID %in% sig.oids ) %>% 
  mutate(Assay_OlinkID=paste0(Assay,"_",OlinkID)) %>% 
  select(Assay_OlinkID,SampleID,NPX,COVID,age.group) %>% 
  group_by(Assay_OlinkID) %>% 
  mutate(NPX=(NPX-mean(NPX,na.rm=T))/sd(NPX,na.rm=T)) %>% 
  mutate(NPX= ifelse(abs(NPX)>4,4*sign(NPX),NPX)) %>% ##Replace any values beyond 4sd with 4
  mutate(NPX=ifelse(is.na(NPX)|!is.finite(NPX),0,NPX)) %>% ##if missing, replace with mean/0
  ungroup() %>% 
  pivot_wider(names_from="Assay_OlinkID",values_from="NPX") %>% 
  mutate(COVID=factor(COVID),
         age.group=factor(age.group))  

tmp.wide <- tmp.wide %>% 
  mutate(covid.color=c("#FE1F04FF","#00C7E1FF")[as.numeric(tmp.wide$COVID)],
         age.color=c("#00559EFF","#FFC700FF")[as.numeric(tmp.wide$age.group)])

rowside.cols <- cbind(COVID=tmp.wide$covid.color,
                      Age_Group = tmp.wide$age.color)

x <- as.matrix(tmp.wide %>% select(-c(SampleID,COVID,covid.color,age.group,age.color)))
rownames(x) <- tmp.wide$SampleID


CairoPDF(paste0("sample heatmaps/samples_heatmap_all_proteins.pdf"),
         width=20,
         height=20)

tmp.covid <- tmp.wide %>% select(COVID,covid.color) %>% unique() %>% arrange(COVID)
tmp.age.group <- tmp.wide %>% select(age.group,age.color) %>% unique() %>% arrange(age.group)

tmp.mod <- hclust(as.dist(1-cor(t(x),method="pearson")))
cuth <- getCutHeight(tmp.mod,search.seq = cut.seq,n.clusters=pref.clusters)


all.prots.sample.clusters <- cutree(tmp.mod,k=2)

min.col <- which.min(min(x) > seq(range(x)[1],range(x)[2],length.out = 1024))
max.col <- which.min(max(x) > seq(range(x)[1],range(x)[2],length.out = 1024))


heatmap3(t(x),
         scale="none",
         ColSideColors=(rowside.cols),
         method="complete",
         ColSideCut = cuth,
         distfunC = function(x) as.dist(1-cor(t(x),method="pearson")),
         distfunR = function(x) as.dist(1-cor(t(x),method="pearson")),
         col = heatmap.color.scale[min.col:max.col],
         legendfun =  function() {showLegend(legend = c("COVID",as.character(tmp.covid$COVID),
                                                        "Age Group",as.character(tmp.age.group$age.group)),
                                             col=c("white",tmp.covid$covid.color,
                                                   "white",tmp.age.group$age.color),
                                             cex=2,
                                             lwd=5)}
)
##Color bar
plotfunctions::gradientLegend(valRange=c(-4,4),
                              color = heatmap.color.scale,
                              pos=c(0.11,.85,0.13,1.02), # xleft, ybottom, xright, ytop
                              fit.margin = F,cex=1.75)

dev.off()


### significant proteins 

sig.oids <- covid.lm.results %>% filter(Threshold=="Significant") %>% pull(OlinkID)

tmp.wide <- npx.info.covid %>% 
  filter(OlinkID %in% sig.oids ) %>% 
  mutate(Assay_OlinkID=paste0(Assay,"_",OlinkID)) %>% 
  select(Assay_OlinkID,SampleID,NPX,COVID,age.group) %>% 
  group_by(Assay_OlinkID) %>% 
  mutate(NPX=(NPX-mean(NPX,na.rm=T))/sd(NPX,na.rm=T)) %>% 
  mutate(NPX= ifelse(abs(NPX)>4,4*sign(NPX),NPX)) %>% ##Replace any values beyond 4sd with 4
  mutate(NPX=ifelse(is.na(NPX)|!is.finite(NPX),0,NPX)) %>% ##if missing, replace with mean/0
  ungroup() %>% 
  pivot_wider(names_from="Assay_OlinkID",values_from="NPX") %>% 
  mutate(COVID=factor(COVID),
         age.group=factor(age.group))  

tmp.wide <- tmp.wide %>% 
  mutate(covid.color=c("#FE1F04FF","#00C7E1FF")[as.numeric(tmp.wide$COVID)],
         age.color=c("#00559EFF","#FFC700FF")[as.numeric(tmp.wide$age.group)])

rowside.cols <- cbind(COVID=tmp.wide$covid.color,
                      Age_Group = tmp.wide$age.color)

x <- as.matrix(tmp.wide %>% select(-c(SampleID,COVID,covid.color,age.group,age.color)))
rownames(x) <- tmp.wide$SampleID

CairoPDF(paste0("sample heatmaps/samples_heatmap_significant_proteins.pdf"),
         width=20,
         height=20)

tmp.covid <- tmp.wide %>% select(COVID,covid.color) %>% unique() %>% arrange(COVID)
tmp.age.group <- tmp.wide %>% select(age.group,age.color) %>% unique() %>% arrange(age.group)

tmp.mod <- hclust(as.dist(1-cor(t(x),method="pearson")))
cuth <- getCutHeight(tmp.mod,search.seq = cut.seq,n.clusters=pref.clusters)
sig.prots.sample.clusters <- cutree(tmp.mod,k=2)

min.col <- which.min(min(x) > seq(range(x)[1],range(x)[2],length.out = 1024))
max.col <- which.min(max(x) > seq(range(x)[1],range(x)[2],length.out = 1024))

heatmap3(t(x),
         scale="none",
         ColSideColors=(rowside.cols),
         method="complete",
         ColSideCut = cuth,
         distfunC = function(x) as.dist(1-cor(t(x),method="pearson")),
         distfunR = function(x) as.dist(1-cor(t(x),method="pearson")),
         col = heatmap.color.scale[min.col:max.col],
         legendfun =  function() {showLegend(legend = c("COVID",as.character(tmp.covid$COVID),
                                                        "Age Group",as.character(tmp.age.group$age.group)),
                                             col=c("white",tmp.covid$covid.color,
                                                   "white",tmp.age.group$age.color),
                                             cex=2,
                                             lwd=5)}
)
##Color bar
plotfunctions::gradientLegend(valRange=c(-4,4),
                              color = heatmap.color.scale,
                              pos=c(0.11,.85,0.13,1.02), # xleft, ybottom, xright, ytop
                              fit.margin = F,cex=1.75)

dev.off()






### Extract patient clusters
all(names(sig.prots.sample.clusters ) == names(all.prots.sample.clusters))

cluster.out <- tibble(all_proteins_cluster=all.prots.sample.clusters,
                     sig_proteins_cluster=sig.prots.sample.clusters,
                     SampleID=names(sig.prots.sample.clusters))
cluster.out <- tmp.wide %>% 
  select(SampleID,age.group,COVID) %>% 
  left_join(cluster.out,"SampleID")






### Top N significant proteins 

topn <- c(25,50,100,200)

for(useN in topn){
  
  sig.oids <- covid.lm.results %>% filter(Threshold=="Significant") %>% arrange(Adjusted_pval) %>% 
    slice(1:useN) %>% pull(OlinkID)
  
  tmp.wide <- npx.info.covid %>% 
    filter(OlinkID %in% sig.oids ) %>% 
    mutate(Assay_OlinkID=paste0(Assay,"_",OlinkID)) %>% 
    select(Assay_OlinkID,SampleID,NPX,COVID,age.group) %>% 
    group_by(Assay_OlinkID) %>% 
    mutate(NPX=(NPX-mean(NPX,na.rm=T))/sd(NPX,na.rm=T)) %>% 
    mutate(NPX= ifelse(abs(NPX)>4,4*sign(NPX),NPX)) %>% ##Replace any values beyond 4sd with 4
    mutate(NPX=ifelse(is.na(NPX)|!is.finite(NPX),0,NPX)) %>% ##if missing, replace with mean/0
    ungroup() %>% 
    pivot_wider(names_from="Assay_OlinkID",values_from="NPX") %>% 
    mutate(COVID=factor(COVID),
           age.group=factor(age.group))  
  
  tmp.wide <- tmp.wide %>% 
    mutate(covid.color=c("#FE1F04FF","#00C7E1FF")[as.numeric(tmp.wide$COVID)],
           age.color=c("#00559EFF","#FFC700FF")[as.numeric(tmp.wide$age.group)])
  
  rowside.cols <- cbind(COVID=tmp.wide$covid.color,
                        Age_Group = tmp.wide$age.color)
  
  x <- as.matrix(tmp.wide %>% select(-c(SampleID,COVID,covid.color,age.group,age.color)))
  rownames(x) <- tmp.wide$SampleID
  
  CairoPDF(paste0("sample heatmaps/samples_heatmap_top",useN,"_significant.pdf"),
           width=20,
           height=20)
  
  tmp.covid <- tmp.wide %>% select(COVID,covid.color) %>% unique() %>% arrange(COVID)
  tmp.age.group <- tmp.wide %>% select(age.group,age.color) %>% unique() %>% arrange(age.group)
  
  tmp.mod <- hclust(as.dist(1-cor(t(x),method="pearson")))
  cuth <- getCutHeight(tmp.mod,search.seq = cut.seq,n.clusters=pref.clusters)
  
  cluster.out[[paste0("cluster_top_",useN,"_sig")]] <- cutree(tmp.mod,k=2)
  
  min.col <- which.min(min(x) > seq(range(x)[1],range(x)[2],length.out = 1024))
  max.col <- which.min(max(x) > seq(range(x)[1],range(x)[2],length.out = 1024))
  
  heatmap3(t(x),
           scale="none",
           ColSideColors=(rowside.cols),
           method="complete",
           ColSideCut = cuth,
           distfunC = function(x) as.dist(1-cor(t(x),method="pearson")),
           distfunR = function(x) as.dist(1-cor(t(x),method="pearson")),
           col = heatmap.color.scale[min.col:max.col],
           legendfun =  function() {showLegend(legend = c("COVID",as.character(tmp.covid$COVID),
                                                          "Age Group",as.character(tmp.age.group$age.group)),
                                               col=c("white",tmp.covid$covid.color,
                                                     "white",tmp.age.group$age.color),
                                               cex=2,
                                               lwd=5)}
  )
  ##Color bar
  plotfunctions::gradientLegend(valRange=c(-4,4),
                                color = heatmap.color.scale,
                                pos=c(0.11,.85,0.13,1.02), # xleft, ybottom, xright, ytop
                                fit.margin = F,cex=1.75)
  
  dev.off()
  
  
}


olinkr::output_dfs_to_excel(list(cluster.out),
                            df_names = "COVID clusters",
                            filename="sample heatmaps/sample_clusters.xlsx")








# use Residuals -----------------------------------------------------------

top50.label <- covid.lm.results %>% filter(Threshold=="Significant") %>% arrange(-abs(estimate)) %>% 
  slice(1:50) %>% mutate(OlinkID = paste0(Assay,"_",OlinkID)) %>% pull(OlinkID)

### Top N significant proteins using residuals
single_lm_return_residuals <- function(formula,data){
  mod <- lm(as.formula(formula),data=data)
  data$lm.residual <- residuals(mod)
  return(data)
}

npx.info.resids <- npx.info.covid %>% 
  group_by(OlinkID) %>% 
  mutate(NPX=ifelse(is.na(NPX),median(NPX,na.rm=T),NPX)) %>% 
  group_modify(~single_lm_return_residuals(formula=paste0("NPX~",paste(lmm.covariates,collapse="+")),data=.x)) %>% 
  ungroup()

topn <- c(25,50,100,200,sum(covid.lm.results$Threshold=="Significant"),nrow(covid.lm.results))

for(useN in topn){
  
  sig.oids <- covid.lm.results %>% filter(Threshold=="Significant") %>% arrange(Adjusted_pval) %>% 
    slice(1:useN) %>% pull(OlinkID)
  
  tmp.wide <- npx.info.resids %>% mutate(NPX=lm.residual) %>% 
    filter(OlinkID %in% sig.oids ) %>% 
    mutate(Assay_OlinkID=paste0(Assay,"_",OlinkID)) %>% 
    select(Assay_OlinkID,SampleID,NPX,COVID,age.group) %>% 
    group_by(Assay_OlinkID) %>% 
    mutate(NPX=(NPX-mean(NPX,na.rm=T))/sd(NPX,na.rm=T)) %>% 
    mutate(NPX= ifelse(abs(NPX)>4,4*sign(NPX),NPX)) %>% ##Replace any values beyond 4sd with 4
    mutate(NPX=ifelse(is.na(NPX)|!is.finite(NPX),0,NPX)) %>% ##if missing, replace with mean/0
    ungroup() %>% 
    pivot_wider(names_from="Assay_OlinkID",values_from="NPX") %>% 
    mutate(COVID=factor(COVID),
           age.group=factor(age.group))  
  
  tmp.wide <- tmp.wide %>% 
    mutate(covid.color=c("#FE1F04FF","#00C7E1FF")[as.numeric(tmp.wide$COVID)],
           age.color=c("#00559EFF","#FFC700FF")[as.numeric(tmp.wide$age.group)])
  
  rowside.cols <- cbind(COVID=tmp.wide$covid.color,
                        Age_Group = tmp.wide$age.color)
  
  x <- as.matrix(tmp.wide %>% select(-c(SampleID,COVID,covid.color,age.group,age.color)))
  rownames(x) <- tmp.wide$SampleID
  
  CairoPDF(paste0("sample heatmaps/on_residuals/samples_heatmap_top",useN,"_useResids.pdf"),
           width=20,
           height=20)
  
  tmp.covid <- tmp.wide %>% select(COVID,covid.color) %>% unique() %>% arrange(COVID)
  tmp.age.group <- tmp.wide %>% select(age.group,age.color) %>% unique() %>% arrange(age.group)
  
  tmp.mod <- hclust(as.dist(1-cor(t(x),method="pearson")))
  cuth <- getCutHeight(tmp.mod,search.seq = cut.seq,n.clusters=2)
  
  cluster.out[[paste0("cluster_top_",useN,"_sig")]] <- cutree(tmp.mod,k=2)
  
  min.col <- which.min(min(x) > seq(range(x)[1],range(x)[2],length.out = 1024))
  max.col <- which.min(max(x) > seq(range(x)[1],range(x)[2],length.out = 1024))
  
  tmp.row.labels <- rep("",nrow(t(x)))
  tmp.row.labels[rownames(t(x)) %in% top50.label] <- rownames(t(x))[rownames(t(x)) %in% top50.label]
  
  heatmap3(t(x),
           scale="none",
           ColSideColors=(rowside.cols),
           method="complete",
           labRow = tmp.row.labels,
           ColSideCut = cuth,
           # cexRow = 1.5,
           margins = c(5,7.5),
           distfunC = function(x) as.dist(1-cor(t(x),method="pearson")),
           distfunR = function(x) as.dist(1-cor(t(x),method="pearson")),
           col = heatmap.color.scale[min.col:max.col],
           legendfun =  function() {showLegend(legend = c("COVID",as.character(tmp.covid$COVID),
                                                          "Age Group",as.character(tmp.age.group$age.group)),
                                               col=c("white",tmp.covid$covid.color,
                                                     "white",tmp.age.group$age.color),
                                               cex=2,
                                               lwd=5)}
  )
  ##Color bar
  plotfunctions::gradientLegend(valRange=c(-4,4),
                                color = heatmap.color.scale,
                                pos=c(0.11,.85,0.13,1.02), # xleft, ybottom, xright, ytop
                                fit.margin = F,cex=1.75)
  
  dev.off()
  
  
}



cut.seq <- seq(.5,1.8,.01)
pref.clusters <- 8

cluster.list <- list()


##all proteins
tmp <- npx.info.resids %>% mutate(NPX=lm.residual) %>% 
  mutate(Assay_OlinkID=paste0(Assay,"_",OlinkID)) %>% 
  select(Assay_OlinkID,SampleID,NPX) %>% 
  group_by(Assay_OlinkID) %>% 
  mutate(NPX=(NPX-mean(NPX,na.rm=T))/sd(NPX,na.rm=T)) %>% 
  mutate(NPX=ifelse(is.na(NPX)|!is.finite(NPX),0,NPX)) %>% ##18 values are missing. replace with mean (scaled data, mean=0)
  ungroup() %>% 
  pivot_wider(names_from="Assay_OlinkID",values_from="NPX") %>% 
  select(-SampleID) %>% 
  cor(method="pearson")
min.col <- which.min(min(tmp) > seq(-1,1,length.out = 1024))
max.col <- which.min(max(tmp) > seq(-1,1,length.out = 1024))


##### All proteins

CairoPDF("assay correlations/on_residuals/all_proteins_corr_heatmap_useResids.pdf",
         width=20,
         height=20)


tmp.mod <- hclust(as.dist(1-as.matrix(tmp)))
cuth <- getCutHeight(tmp.mod,search.seq = cut.seq,n.clusters=pref.clusters)
out.clusters <- cutree(tmp.mod,h=cuth)
cluster.list[["All Proteins"]] <- data.frame(Assay = names(out.clusters),
                                             Cluster=as.character(as.numeric(out.clusters))) %>% arrange(Cluster)

heatmap3(tmp,
         scale="none",
         method="complete",
         revC=T,
         ColSideCut = cuth,
         distfun = function(x) as.dist(1 - (x)),
         col = heatmap.color.scale[min.col:max.col],
         legendfun =  function() {showLegend(legend = "",col="white")})
plotfunctions::gradientLegend(valRange=c(-1,1),
                              color = heatmap.color.scale,
                              pos=c(0,.85,.025,.99), # xleft, ybottom, xright, ytop
                              fit.margin = F,cex=2)
dev.off()




##Significant proteins
tmp <- npx.info.resids %>% mutate(NPX=lm.residual) %>% 
  filter(OlinkID %in% covid.lm.results$OlinkID[covid.lm.results$Threshold=="Significant"]) %>% 
  mutate(Assay_OlinkID=paste0(Assay,"_",OlinkID)) %>% 
  select(Assay_OlinkID,SampleID,NPX) %>% 
  group_by(Assay_OlinkID) %>% 
  mutate(NPX=(NPX-mean(NPX,na.rm=T))/sd(NPX,na.rm=T)) %>% 
  mutate(NPX=ifelse(is.na(NPX)|!is.finite(NPX),0,NPX)) %>% ##18 values are missing. replace with mean (scaled data, mean=0)
  ungroup() %>% 
  pivot_wider(names_from="Assay_OlinkID",values_from="NPX") %>% 
  select(-SampleID) %>% 
  cor(method="pearson")
min.col <- which.min(min(tmp) > seq(-1,1,length.out = 1024))
max.col <- which.min(max(tmp) > seq(-1,1,length.out = 1024))


CairoPDF("assay correlations/on_residuals/significant_proteins_corr_heatmap_useResids.pdf",
         width=20,
         height=20)


tmp.mod <- hclust(as.dist(1-as.matrix(tmp)))
cuth <- getCutHeight(tmp.mod,search.seq = cut.seq,n.clusters=pref.clusters)
out.clusters <- cutree(tmp.mod,h=cuth)
cluster.list[["Significant Proteins"]] <- data.frame(Assay = names(out.clusters),
                                                     Cluster=as.character(as.numeric(out.clusters))) %>% arrange(Cluster)

heatmap3(tmp,
         scale="none",
         method="complete",
         revC=T,
         ColSideCut = cuth,
         distfun = function(x) as.dist(1 - (x)),
         col = heatmap.color.scale[min.col:max.col],
         legendfun =  function() {showLegend(legend = "",col="white")})
plotfunctions::gradientLegend(valRange=c(-1,1),
                              color = heatmap.color.scale,
                              pos=c(0,.85,.025,.99), # xleft, ybottom, xright, ytop
                              fit.margin = F,cex=2)
dev.off()




olinkr::output_dfs_to_excel(cluster.list,
                            filename="assay correlations/on_residuals/clusters_Proteins_useResids.xlsx",
                            df_names = names(cluster.list))



#### Upregulated in COVID
tmp <- npx.info.resids %>% mutate(NPX=lm.residual) %>% 
  filter(OlinkID %in% covid.lm.results$OlinkID[covid.lm.results$Threshold=="Significant" &
                                                 covid.lm.results$estimate>0]) %>%
  mutate(Assay_OlinkID=paste0(Assay,"_",OlinkID)) %>%
  select(Assay_OlinkID,SampleID,NPX) %>%
  group_by(Assay_OlinkID) %>%
  mutate(NPX=(NPX-mean(NPX,na.rm=T))/sd(NPX,na.rm=T)) %>%
  mutate(NPX=ifelse(is.na(NPX)|!is.finite(NPX),0,NPX)) %>% ##18 values are missing. replace with mean (scaled data, mean=0)
  ungroup() %>%
  pivot_wider(names_from="Assay_OlinkID",values_from="NPX") %>%
  select(-SampleID) %>%
  cor(method="pearson")
min.col <- which.min(min(tmp) > seq(-1,1,length.out = 1024))
max.col <- which.min(max(tmp) > seq(-1,1,length.out = 1024))
tmp.mod <- hclust(as.dist(1-as.matrix(tmp)))
cuth <- getCutHeight(tmp.mod,search.seq = cut.seq,n.clusters=8)

CairoPDF("assay correlations/on_residuals/upregulated_proteins_corr_heatmap_useResids.pdf",
         width=20,
         height=20)
heatmap3(tmp,
         scale="none",
         method="complete",
         revC=T,
         ColSideCut = cuth,
         distfun = function(x) as.dist(1 - (x)),
         col = heatmap.color.scale[min.col:max.col],
         legendfun =  function() {showLegend(legend = "",col="white")})
plotfunctions::gradientLegend(valRange=c(-1,1),
                              color = heatmap.color.scale,
                              pos=c(0,.85,.025,.99), # xleft, ybottom, xright, ytop
                              fit.margin = F,cex=2)
dev.off()
#### Down regulated in COVID
tmp <- npx.info.resids %>% mutate(NPX=lm.residual) %>% 
  filter(OlinkID %in% covid.lm.results$OlinkID[covid.lm.results$Threshold=="Significant" &
                                                 covid.lm.results$estimate<0]) %>%
  mutate(Assay_OlinkID=paste0(Assay,"_",OlinkID)) %>%
  select(Assay_OlinkID,SampleID,NPX) %>%
  group_by(Assay_OlinkID) %>%
  mutate(NPX=(NPX-mean(NPX,na.rm=T))/sd(NPX,na.rm=T)) %>%
  mutate(NPX=ifelse(is.na(NPX)|!is.finite(NPX),0,NPX)) %>% ##18 values are missing. replace with mean (scaled data, mean=0)
  ungroup() %>%
  pivot_wider(names_from="Assay_OlinkID",values_from="NPX") %>%
  select(-SampleID) %>%
  cor(method="pearson")
min.col <- which.min(min(tmp) > seq(-1,1,length.out = 1024))
max.col <- which.min(max(tmp) > seq(-1,1,length.out = 1024))
tmp.mod <- hclust(as.dist(1-as.matrix(tmp)))
cuth <- getCutHeight(tmp.mod,search.seq = cut.seq,n.clusters=8)

CairoPDF("assay correlations//on_residuals/downregulated_proteins_corr_heatmap_useResids.pdf",
         width=20,
         height=20)
heatmap3(tmp,
         scale="none",
         method="complete",
         revC=T,
         ColSideCut = cuth,
         distfun = function(x) as.dist(1 - (x)),
         col = heatmap.color.scale[min.col:max.col],
         legendfun =  function() {showLegend(legend = "",col="white")})
plotfunctions::gradientLegend(valRange=c(-1,1),
                              color = heatmap.color.scale,
                              pos=c(0,.85,.025,.99), # xleft, ybottom, xright, ytop
                              fit.margin = F,cex=2)
dev.off()
