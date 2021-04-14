library(OlinkAnalyze)
library(corrplot)
library(heatmap3)
library(Cairo)
# library(parallel)



#load image rather than running lines 12:292
rm(list=ls())
# load("20201023_Analysis_Full_Covariates.RData")




lmm.covariates <- c("heart.disease","Diabetes","HTN","HLD","Pulmonary","kidney.disease","Immuno.comprimised","sex","ethnicity","age.decile")

load("../../NPX data/20201023_clean_merged_data")
source("../heatmap_functions.R")

npx.info <- npx.info %>% 
  filter(COVID=="Positive") %>% 
  filter(Timepoint != "DE")


# 
# quickLmer <- function(formula,data){
#   mod <- anova(lmer(formula = formula,data=data))
#   out <- matrix(as.numeric(mod$`Pr(>F)`),nrow=1)
#   colnames(out) <- rownames(mod)
#   return(as.data.frame(out))
# }
# permuteWHOLMM <- function(formula,data){
#   
#   # orig.who.map$perm.who <- sample(orig.who.map$WHO.max)
#   data <- data %>% left_join(orig.who.map %>% 
#                                mutate(perm.who=sample(WHO.max)) %>% 
#                                select(Patient,perm.who),by="Patient")
#   
#   out <- data %>% 
#     group_by(OlinkID) %>% 
#     group_modify(~quickLmer(as.formula(formula),
#                             data=.x)) %>% 
#     ungroup()
#   
#   return(out)
#   
# }
# 
# 
# perm.who.results <- vector("list",100)
# orig.who.map <- npx.info %>% select(Patient,WHO.max) %>% unique()
# 
# for(i in 1:100){
#   perm.who.results[[i]] <- permuteWHOLMM(formula=paste("NPX~Timepoint*perm.who + ",paste(lmm.covariates,collapse="+"),"+(1|Patient)"),
#                                          data=npx.info)
#   cat(i," ", as.character(Sys.time()),"\n")
# }
# 
# 
# saveRDS(perm.who.results,
#         file="perumute_whomax.RDS")

orig.who.results <- read_excel("/Users/jamey.guess/Desktop/20201023_Broad_Hacohen - Sponsored Covid/R code/20201023_analysis/lmm/20201023_LMER_all_covariates_results.xlsx",
                               sheet=4)

perm.who.results <- readRDS("perumute_whomax.RDS")

names(perm.who.results) <- as.character(1:length(perm.who.results))
perm.who.results.df <- plyr::ldply(perm.who.results,data.frame,.id="Rep") %>% as_tibble()









  



orig.who.results %>% 
  mutate(sig=p.value<0.05) %>% 
  select(term,sig) %>% 
  table()


line.dat <- data.frame(grp=c("Max WHO","Time:WHO Interaction"),
                       v=c(1131,963))


pval.hists <- perm.who.results.df %>% 
  pivot_longer(-c(Rep,OlinkID),names_to="term",values_to="p.value") %>% 
  filter(term %in% c("perm.who","Timepoint.perm.who")) %>% 
  group_by(Rep,term) %>% 
  # mutate(BH.pvalue=p.adjust(p.value,method ="BH")) %>% 
  mutate(BH.pvalue=p.adjust(p.value,method="BH")) %>% 
  group_by(Rep,term) %>% 
  summarize(N.Sig = sum(p.value<0.05)) %>% 
  ungroup() %>% 
  mutate(grp = ifelse(term=="perm.who","Max WHO","Time:WHO Interaction")) %>% 
  ggplot(aes(N.Sig))+
  geom_histogram(color="white",bins=35)+
  facet_wrap(~grp,ncol=1)+
  geom_vline(data=line.dat,aes(xintercept=v),color="red")+
  theme_bw()+
  labs(x="Number of p-values < 0.05",y="Count")
  
ggsave(filename="lmm/permutation_pvalue.pdf",
       plot=pval.hists,
       height=6,
       width=8,
       units = "in",
       device = cairo_pdf)



tmp.pal <- c("#00C7E1FF", "#FE1F04FF")




hist.orig.data <- orig.who.results %>% 
  filter(term %in% c("WHO.max","Timepoint:WHO.max")) %>% 
  mutate(cuts = 0.025*(as.numeric(cut(p.value,breaks = seq(0,1,by=.025),ordered_result = T)))-.025/2) %>% 
  mutate(cuts.bh = cut(Adjusted_pval,breaks = seq(0,1,by=.025))) %>% 
  group_by(term,cuts) %>% 
  summarize(n=n()) %>% 
  ungroup %>% 
  group_by(term,cuts) %>% 
  ungroup() %>% 
  mutate(Data="Original")



perm.hist95 <- perm.who.results.df %>% 
  pivot_longer(-c(Rep,OlinkID),names_to="term",values_to="p.value") %>% 
  filter(term %in% c("perm.who","Timepoint.perm.who")) %>% 
  group_by(Rep) %>% 
  mutate(BH.pvalue=p.adjust(p.value,method ="BH")) %>% 
  ungroup() %>% 
  mutate(cuts = 0.025*(as.numeric(cut(p.value,breaks = seq(0,1,by=.025),ordered_result = T)))-.025/2) %>% 
  mutate(cuts.bh = cut(BH.pvalue,breaks = seq(0,1,by=.025))) %>% 
  group_by(Rep,term,cuts) %>% 
  summarize(n=n()) %>% 
  ungroup %>% 
  group_by(term,cuts) %>% 
  summarize(n=quantile(n,.95)) %>% 
  ungroup() %>% 
  mutate(Data="Permuted") %>% 
  rbind(hist.orig.data) %>% 
  mutate(grp = ifelse(term=="perm.who","Max WHO",ifelse(term=="WHO.max","Max WHO","Time:WHO Interaction"))) %>%
  ggplot(aes(x=cuts,y=n,fill=Data))+
  geom_bar(stat="identity",color=NA,alpha=0.4) +
  facet_wrap(~grp,scales="fixed",ncol=1)+
  theme_bw()+
  scale_fill_manual(values=tmp.pal)+
  labs(x="p-value",y="Count",fill="WHO")+
  theme(legend.position="bottom")


ggsave(filename="lmm/permutation_histograms.pdf",
       plot=perm.hist95,
       height=6,
       width=8,
       units = "in",
       device = cairo_pdf)



# perm.who.results.df %>% 
#   pivot_longer(-c(Rep,OlinkID),names_to="term",values_to="p.value") %>% 
#   filter(term %in% c("perm.who","Timepoint.perm.who")) %>% 
#   group_by(Rep) %>% 
#   mutate(BH.pvalue=p.adjust(p.value,method ="BH")) %>% 
#   ungroup() %>% 
#   mutate(cuts = 0.025*(as.numeric(cut(p.value,breaks = seq(0,1,by=.025),ordered_result = T)))-.025/2) %>% 
#   mutate(cuts.bh = cut(BH.pvalue,breaks = seq(0,1,by=.025))) %>% 
#   select(term,p.value) %>% 
#   mutate(Data="Permuted") %>% 
#   rbind(orig.who.results %>% 
#           filter(term %in% c("WHO.max","Timepoint:WHO.max")) %>% 
#           select(term,p.value) %>% mutate(Data="Original")) %>% 
#   mutate(grp = ifelse(term=="perm.who","Max WHO",ifelse(term=="WHO.max","Max WHO","Time:WHO Interaction"))) %>%
#   ggplot(aes(x=p.value,color=Data))+
#   geom_density() + 
#   facet_wrap(~grp,scales="fixed",ncol=1)+
#   theme_bw()+
#   scale_fill_manual(values=tmp.pal)+
#   labs(x="p-value",y="Count",fill="WHO")+
#   theme(legend.position="bottom")
#   


