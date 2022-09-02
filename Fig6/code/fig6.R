# data wrangling
library(tidyverse)
library(janitor)
library(dplyr)

# misc
library(lubridate)
library(conflicted)
library(Hmisc)

library(broom)
library(hrbrthemes)
library(here)

#PDF import
library(magick)
library(pdftools)
library(grImport2)

# other plotting
library(ggbeeswarm)
library(gridExtra)
library(ggpubr)
library(cowplot)
library(magrittr)
library(ggtext)
library(ZamanianLabThemes)


conflict_prefer("filter", "dplyr")

#set working directory-------------------------------- 
setwd<-here() 

#-------------------------------------read in data--------------------------------------
bodipy.data <- readRDS("../data/bodipy.rds")

#set theme 
aesthetic <- 
  theme_zlab() +
  theme(
    panel.grid.major.y = element_line(),
    axis.text.x = element_markdown(angle = 45, hjust = 1, face="plain",size=10),
    axis.text.y = element_text(face = "plain", size = 10),
    strip.text.x = element_markdown(size=10),
    strip.text.y = element_markdown(size=12),
    strip.text = ggplot2::element_text(face = "plain", size = 10),
    axis.title.y = element_text(size=12, face= "plain"),
    panel.border = element_rect(colour= "black", fill=NA))


#-----------------------------------optimization BODIPY+HB101 assay---------------------------------

fig6.opt<- bodipy.data%>%
  filter( Metadata_Date == "20210823", other != "GreenBodipy")

#rename and set order
fig6.opt$other<-factor(fig6.opt$other, levels = c("RedBodipy_1ulHB", "RedBodipy_1ulM9","RedBodipy_2ulHB", "RedBodipy_2ulM9"),
                       labels=c("BODIPY1x+HB101","BODIPY1x+M9" ,"BODIPY2x+HB101","BODIPY2x+M9"))

fig6.opt$treatment<-factor(fig6.opt$treatment, levels = c("CON", "5HT","5HT_NIC", "HB101","HB101_NIC"), 
                           labels= c("CON", "5HT","5HT+NIC", "HB101","HB101+NIC"))

#stats 
stats<-compare_means(Intensity_StdIntensity_binarymaskTx~treatment , data = fig6.opt, group.by= "other",
                     method = "t.test", ref.group = "CON")
stats$p.signif<- factor(stats$p.signif, levels=c("ns","*","**","***","****"),labels= c("","*","**","***","***"))

stats2<-compare_means(Intensity_MassDisplacement_binarymaskTx~treatment , data = fig6.opt, group.by= "other",
                      method = "t.test", ref.group = "5HT")

stats2$p.signif<- factor(stats2$p.signif, levels=c("ns","*","**","***","****"),labels= c("","*","**","***","***"))

#plot
(bodipy.data.plot <- ggplot(filter(fig6.opt), aes(treatment, y =Intensity_StdIntensity_binarymaskTx )) +
    geom_errorbar(aes(color = other),
                  position = position_dodge(width = 0.97),
                  stat = "summary", width = 0.5, alpha = 0.75 ) +
    geom_bar(aes(fill = treatment, color=treatment), position = position_dodge(width = 0.97), stat = "summary", fun = "mean", alpha=.8) +
    geom_quasirandom(aes(color = strains),
                     shape = 21, fill = "white", size = 1,
                     alpha = 0.75, dodge.width = 0.9, width = 0.1)+
    stat_pvalue_manual(
      stats,x="group2", y.position = .018,
      label = "p.signif",
      hide.ns = TRUE)+
    stat_pvalue_manual(
      stats2,x="group2", y.position = .020,
      label = "p.signif",color="brown"
    )+
    facet_grid(cols = vars(other)) +
    scale_y_continuous(expand = c(0, 0))+
    scale_color_manual(values = c("#000000", "#000000", "#000000","#000000","#000000", "#000000", "#000000","#000000","#000000","#000000"))+
    scale_fill_manual(values =  c("#000000","brown","#F0F0F0", "#F0F0F0", "#F0F0F0","#F0F0F0", "#F0F0F0", "#F0F0F0","#F0F0F0"))+
    labs(x = "", y = "FLU", title = " ") +
    aesthetic +
    theme(plot.margin = margin(0, 0, -5, 0, unit = "pt"))+
    NULL)


#--------------------------------------- BODIPY Cholinergic panel----------------------------
panel.data <- readRDS("../data/bodipy_panel.rds")

# filter GFP pos 
fig6<- panel.data%>%
  mutate(transgenic = ifelse(nonoverlappingworms_Intensity_StdIntensity_binarymaskGFP >= 0.0015, "T", "F"))%>%
  filter(transgenic == "T")#%>%

#get mean red of each treatment for each strain, new table, new column

control_summary <- filter(fig6, treatment == "CON") %>%
  group_by(strains,batch, treatment) %>%
  summarise(mean_red = mean(nonoverlappingworms_Intensity_StdIntensity_binarymaskTx))%>%
  select(-treatment, mean_CON_red = mean_red)

fig6 <- left_join(fig6, control_summary) %>%
  mutate(percent_red = (nonoverlappingworms_Intensity_StdIntensity_binarymaskTx / mean_CON_red))

#divide pumps by set treatment reference (control) pumps

#rename and set order  
fig6$strains<- factor(fig6$strains, levels= c("ZAM19","ZAM20","ZAM11"), labels= c("N2", "*gar-3(gk305)*", "*gar-3(gk305);<br>myo-2p::Bma-gar-3*"))

fig6$treatment <-factor(fig6$treatment, levels= c("CON","5HT",'ACH','CAR','OXO','ATR','NIC','LEV'),
                        labels= c("CON", "5HT",'ACH','CAR','OXO','ATR','NIC','LEV'))

#outlier trimming

fig6<- fig6%>%
  select(1:13,57,61,86,87,88)

trim<-fig6%>%
  select(10,11,15,16:18)

#############
Size_filter<- trim%>%
  filter(straightenedworms_AreaShape_Area >1000)%>%
  filter(straightenedworms_AreaShape_Area <7000)

z_scores <- trim%>%
  group_by(strains, treatment)%>%
  mutate(z_scores = (percent_red-mean(percent_red))/sd(percent_red))

no_outliers <- z_scores%>%
  filter(z_scores <3)%>%
  filter(z_scores >-3)

##########

# stats
stats3<-compare_means(percent_red~treatment , data = no_outliers, group.by=c("strains"),
                      method = "t.test", ref.group = "5HT")

stats3$group2 <-factor(stats3$group2, levels= c("CON", "5HT",'ACH','CAR','OXO','ATR','NIC','LEV'))
stats3$p.signif<- factor(stats3$p.signif, levels=c("ns","*","**","***","****"),labels= c("","*","**","***","***"))

stats4<-compare_means(percent_red~treatment , data = no_outliers, group.by=c("strains"),
                      method = "t.test", ref.group = "CON")

stats4$group2 <-factor(stats4$group2, levels= c("CON", "5HT",'ACH','CAR','OXO','ATR','NIC','LEV'))
stats4$p.signif<- factor(stats4$p.signif, levels=c("ns","*","**","***","****"),labels= c("","*","**","***","***"))



#plot 
plot<- ggplot(no_outliers, aes(x = treatment, y = percent_red)) +
  geom_errorbar(aes(color = strains),
                position = position_dodge(width = 0.97,),
                stat = "summary", width = 0.5, alpha = 0.75 ) +
  geom_bar(aes(fill = treatment, color= strains), position = position_dodge(width = 0.97), stat = "summary", fun = "mean", alpha=.8) +
  geom_quasirandom(aes(color = strains),
                   shape = 21, fill = "white", size = 2,
                   alpha = 0.75, dodge.width = 0.9, width = 0.1)+
  stat_pvalue_manual(
    stats3, x="group2", y.position = 3,
    label = "p.signif", color="brown")+
  stat_pvalue_manual(
    stats4, x="group2", y.position = 3.2,
    label = "p.signif", color="black")+
  facet_grid(cols = vars(strains)) +
  labs(title= "",x = "", y = "% CON FLU") +
  theme_minimal(base_size = 10, base_family = "Helvetica") +
  scale_fill_manual(values =  c("#000000","brown","#F0F0F0", "#F0F0F0", "#F0F0F0","#F0F0F0","#F0F0F0","#F0F0F0","#F0F0F0","#F0F0F0","#F0F0F0", "#F0F0F0", "#F0F0F0","#F0F0F0"))+
  scale_color_manual(values = c("black", "gray","#050298"))+
  scale_y_continuous(labels = scales::percent,limits = c(0,3.8,1), breaks = seq(0, 3.8, 1),expand= c(0,0))+
  geom_hline(yintercept = 1 )+
  aesthetic+
  theme(panel.border = element_rect(color="black", fill=NA),
        strip.text.y = element_blank())+
  
  NULL
plot

#---------------------------------------------Figure 6------------------------------------------ 
illust.hb <- image_read_pdf("/Users/kendra/Library/CloudStorage/Box-Box/ZamanianLab/LabMembers/Kendra/projects/Bma_GAR-ms/schematics/feed.pdf", density = 600)

fig6A <- ggdraw() +
  draw_image(illust.hb, scale = 1)

fig6<-plot_grid(bodipy.data.plot,plot,fig6A, labels = c('A','B','C'), ncol=1,align = "v", axis="rl", rel_widths = c(1,1,1,1), rel_heights = c(1,1,1))

save_plot("./plots/fig6.pdf", fig6, base_height = 8, base_width = 6)

#-------------------------------------supplemental for GFP ID----------------------------------------

worms<-readRDS("../data/bodipy_sup.rds")

worm.corr <- worms %>% 
  select(contains("GFP"),Transgenic)

worm.corr <- worm.corr %>% 
  mutate(temp = 1) %>% 
  select(-FileName_GFP, -PathName_GFP) %>% 
  pivot_longer(cols = !c(temp, Transgenic))%>%
  mutate_at("name", str_replace, "nonoverlappingworms_Intensity_", "")%>%
  mutate_at("name", str_replace, "_binarymaskGFP", "")


(worm_corr_plot <- worm.corr %>% 
    ggplot() +
    geom_quasirandom(aes(x = temp, y = value, color = Transgenic)) +
    facet_wrap(facets = vars(name), scales = 'free_y', ncol=3)+
    scale_color_manual(values = c("#f4c33a",'#8E1F1F'))+
    labs(title='',x = "", y = "") +
    aesthetic+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()))

worm.corr2<-worm.corr%>%
  filter (name== "StdIntensity")

(worm_corr_plot2 <- worm.corr2 %>% 
    ggplot() +
    geom_quasirandom(aes(x = temp, y = value, color = Transgenic)) +
    #facet_wrap(facets = vars(name), scales = 'free_y', ncol=3)+
    scale_color_manual(values = c("#f4c33a",'#8E1F1F'))+
    scale_y_continuous( limits = c(0,.01), breaks = seq(0, .01, .0015),expand = c(0, 0))+
    labs(title='StdIntensity',x = "", y = "") +
    geom_hline(yintercept = 0.0015, color= "#84d6cc")+
    aesthetic+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()))

sup_GFP<-plot_grid(worm_corr_plot,worm_corr_plot2, labels = c('A','B'), ncol=1,rel_widths = c(1.5,.5), rel_heights = c(2,.8))
sup_GFP

save_plot("../sup/worm_gfp_sup.pdf", sup_GFP, base_height = 8, base_width = 8)



