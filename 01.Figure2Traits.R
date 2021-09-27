#############################################
#         Prj: Cotton-Heterosis
#         Assignment: Figure2 phenotype
#         Author: Shanw wang
#         Date: Sep 27,2021
#############################################
###### fig2
install.packages('tidyverse',)
library(tidyverse)
library(ggprism)
#BiocManager::install("ggpol")
library(ggpol)
library(ggsci)
library(ggsignif)
library(patchwork)

# yield -------------------------------------------------------------

yield <- read.table("~/03.project/杂优大论文/02.data/yieldRaw.txt",header = T,sep = "\t")
head(yield)
yield %>% 
  mutate(`LY/plant` = MBN*LY,
         BW = SCY/MBN) -> yield
## box_jitter
boxjitter = function(data,trait) {
  p = data %>% 
    filter(Sample == "B985S" | Sample == "Z41S") %>% 
    mutate(Sample = gsub("S"," × S",Sample)) %>% 
    ggplot(.,aes(x = Type,y = .[,trait]))+
    geom_boxjitter(mapping = aes(fill = Type),
                   errorbar.draw = T,outlier.shape = NA,
                   jitter.size = .8,jitter.alpha = 0.7,alpha = 0.6,lwd = 1)+
    theme_bw()+
    facet_grid(~Sample)+
    geom_signif(test="wilcox.test", 
                comparisons = list(
                  c("Hybrid","Male"),
                  c("Hybrid","Female"),
                  c("Female","Male")), 
                map_signif_level = T,
                size = .8,
                textsize = 4,vjust = 0.2)+
    ylab(trait)+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 12),
          axis.text = element_text(size = 10),
          line = element_line(size = 1),
          panel.border = element_rect(size = 1),
          panel.grid = element_blank(),
          strip.background = element_rect(size = 1),
          strip.text = element_text(size = 12,face = "bold")
    )+
    scale_fill_lancet()
  return(p)
}

p_mbn = boxjitter(data = yield,trait = "MBN")
p_LY.plant = boxjitter(data = yield,trait = "LY/plant")  
p_LP = boxjitter(data = yield,trait = "LP")
p_BW = boxjitter(data = yield,trait = "BW")

p_final <- p_BW+p_mbn+p_LY.plant+p_LP+
  plot_layout(nrow = 1,
              guides = 'collect') & theme(legend.position = "bottom")

ggsave(p_final,filename = "../04.result/01.Figure/02.Figure2/yield.pdf",width = 18,height = 5)

ggsave(plot = p_LY.plant,filename = "../04.result/01.Figure/02.Figure2/LY_yield.pdf",height = 5,width = 8)
ggsave(plot = p_BW,filename = "../04.result/01.Figure/02.Figure2/BW.pdf",height = 5,width = 8)

dev.off()
# quality -----------------------------------------------------------------

quality <- read.table("../02.data/qualityRaw.txt",header = T,sep = "\t")

head(quality)
p_FS = boxjitter(data = quality,trait = "FS")
p_FL = boxjitter(data = quality,trait = "FL")
p_FE = boxjitter(data = quality,trait = "FE")
p_MIC = boxjitter(data = quality,trait = "MIC")

ggsave(plot = p_FS,filename = "../04.result/01.Figure/02.Figure2/FS.pdf",height = 5,width = 8)
ggsave(plot = p_FL,filename = "../04.result/01.Figure/02.Figure2/FL.pdf",height = 5,width = 8)

p_final <- p_FS + p_FL + p_FE + p_FE + p_BW + p_mbn + p_LY.plant + p_LP + plot_layout(ncol = 2, guides = 'collect') & theme(legend.position = "bottom")
ggsave(p_final,filename = "../04.result/01.Figure/02.Figure2/quality.pdf",width = 200,height = 280,units = "mm")


# corr of traits ----------------------------------------------------------

cor.test(as.numeric(yield$MBN),as.numeric(quality$FL) )

yield %>% 
  select(-cultivar) %>% 
  pivot_longer(!c("Sample","Type"),names_to = "trait",values_to = "value") %>% 
  group_by(Sample,Type,trait) %>% 
  summarise(mean = mean(value)) %>% 
  arrange(trait) %>% 
  pivot_wider(names_from = trait,values_from = mean) %>% 
  as.data.frame() %>% 
  filter(Sample == "B985S" | Sample == "Z41S")-> yield_mean

quality%>% 
  pivot_longer(!c("Sample","Type"),names_to = "trait",values_to = "value") %>% 
  group_by(Sample,Type,trait) %>% 
  summarise(mean = mean(value)) %>% 
  arrange(trait) %>% 
  pivot_wider(names_from = trait,values_from = mean)%>% 
  as.data.frame()  %>% 
  filter(Sample == "B985S" | Sample == "Z41S")-> quality_mean

all = cbind(yield_mean,quality_mean[,-c(1,2)])

all_final = all %>% 
  mutate(
    `LY/yield` = LY*MBN,
    BW = SCY/MBN,
    cultivar = paste(Sample,Type,sep = "_")
  ) %>% 
  select(cultivar,LP,MBN,`LY/yield`,BW,FE,FL,FS,MIC) %>% 
  column_to_rownames("cultivar") %>% 
  as.matrix()

library(corrplot)
library(RColorBrewer)
trait_cor = cor(all_final)
trait_cor_all = cor.mtest(all_final,met)
col4 <- colorRampPalette(c("#7F0000","red","#FF7F00","yellow","#7FFF7F", 
                           "cyan", "#007FFF", "blue","#00007F"))  
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
                           "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))
pdf("../04.result/01.Figure/02.Figure2/corr_plot.pdf",width = 10,height = 10)
corrplot(trait_cor,type = "upper",method = 'shade',order = "FPC",addCoef.col="grey",
         col = col4(100))
dev.off()


# heterosis ---------------------------------------------------------------


sep = all %>% 
  mutate(
    `LY/yield` = LY*MBN,
    BW = SCY/MBN
  ) %>% 
  select(Sample,Type,LP,MBN,`LY/yield`,BW,FE,FL,FS,MIC) %>% 
  pivot_longer(!c("Sample","Type")) %>%  
  group_by(Sample,name) %>% 
  pivot_wider(names_from = Type,values_from = value) %>% 
  mutate(h_parent = case_when(
    Female > Male ~ Female,
    Female < Male ~ Male
  ),
  m_parent = (Female+Male)/2
  ) %>% 
  mutate(
    MPH = (Hybrid-m_parent)/m_parent,
    OPH = (Hybrid-h_parent)/h_parent
  ) %>% 
  mutate(
    trait_type = case_when(
      name == "FS" | name == "FL" | name == "FE" | name == "MIC" ~ "q",
      TRUE ~ "y"
    )
  )
qq = sep %>% 
  filter(trait_type == "q") %>% 
  select(Sample,name,MPH,OPH) %>% 
  pivot_longer(!c(Sample,name),names_to = "Heterosis",values_to = "value") %>% 
  ggplot(mapping = aes(x = name,y = value,fill = Heterosis))+
  geom_bar(stat = "identity", position = "dodge", width = 0.5)+
  facet_wrap("Sample")+
  theme_bw()+   
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(size = 10),
        line = element_line(size = 1),
        panel.border = element_rect(size = 1),
        panel.grid = element_blank(),
        strip.background = element_rect(size = 1),
        strip.text = element_text(size = 12,face = "bold")
  )+
  scale_fill_lancet(alpha = 0.7)

yy = sep %>% 
  filter(trait_type == "y") %>% 
  select(Sample,name,MPH,OPH) %>% 
  pivot_longer(!c(Sample,name),names_to = "Heterosis",values_to = "value") %>% 
  ggplot(mapping = aes(x = name,y = value,fill = Heterosis))+
  geom_bar(stat = "identity", position = "dodge", width = 0.5)+
  facet_wrap("Sample")+
  theme_bw()+   
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(size = 10),
        line = element_line(size = 1),
        panel.border = element_rect(size = 1),
        panel.grid = element_blank(),
        strip.background = element_rect(size = 1),
        strip.text = element_text(size = 12,face = "bold")
  )+
  scale_fill_lancet(alpha = 0.7)

ggsave(qq,filename = "../04.result/01.Figure/02.Figure2/heterosis.quality.pdf",width = 5,height = 3)
ggsave(yy,filename = "../04.result/01.Figure/02.Figure2/heterosis.yield.pdf",width = 5,height = 3)