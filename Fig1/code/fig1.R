# data wrangling
library(tidyverse)
library(janitor)
library(dplyr)

# other plotting
library(cowplot)
library(ggtext)
library(ggbeeswarm)
library(ZamanianLabThemes)

# misc
library(conflicted)
library(DescTools)
library(magrittr)
library(lubridate)
library(broom)

library(ggpubr)
library(purrr)
library(ggrepel)
library(Hmisc)
library(here)

#PDF import
library(magick)
library(pdftools)
library(grImport2)

#for phylogenetics 
#remotes::install_github("YuLab-SMU/ggtree", force = TRUE)
#remotes::install_github("YuLab-SMU/tidytree", force = TRUE)
#BiocManager::install("treeio")
library("ape")
library("ggtree")
library("tidytree")
library("treeio")

conflict_prefer("filter", "dplyr")
conflict_prefer("summarize", "dplyr")

#------------------------------------set working directory-----------------------------
setwd<- here()
#-------------------------------------read in data--------------------------------------
phylo_data<- read_csv("../data/phylo/gar_phylo.csv")
qPCR_data<- readRDS("../data/qPCR.rds")
f.data <- readRDS("../data/brugia_fecundity.rds") 
m.data<-  readRDS("../data/brugia_motility.rds") 
mf.data <- readRDS("../data/brugia_flow.rds") 
ct.data <- readRDS("../data/brugia_celltox.rds") 

#---------------------------------Phylogenetics----------------------------------------- 
#read in tree
nwk<- read.newick(here( "Fig1/data/phylo/GAR.combined_final.aln.contree"))
nwk<- ape::root(nwk, node = 123)
nwk <- tree_subset(nwk, node = 123, levels_back = 0) 

#turn tree into df and seperate into columns 
tree_data <- as_tibble(nwk)%>%
  left_join(.,phylo_data)%>%
  mutate(gene_ortholog = coalesce(gene,ortholog)) %>%
  unite(tip.label, c("species","gene_ortholog"),sep= "_", remove = "False")%>%
  mutate(
    species = case_when(
      species  == 'Cel' ~ '*C. elegans*',
      species  == 'Bma' ~ '*B. malayi*',
      species  == 'Asu' ~ '*A. suum*',
      species  == 'Hco' ~ '*H. contortus*',
      species  == 'Hsa' ~ '*H. sapiens*',
      species  == 'Sra' ~ '*S. ratti*',
      species == 'Tmu' ~ '*T. muris*'))

#remove isoforms/duplicate nodes
remove.nodes <- c(1,2,3,4,5,6,7,8,9,10,12,13,14,15,16,20,21,22,23,24,25,26,28,30,32,33,43)

nwk.nr <- treeio::drop.tip(nwk, remove.nodes)

#rename species 

tree_data.nr <- as_tibble(nwk.nr) %>%
  left_join(.,phylo_data)%>%
  mutate(gene_ortholog = coalesce(gene,ortholog)) %>%
  unite(tip.label, c("species","gene_ortholog"),sep= "_", remove = "False")%>%
  mutate(
    species = case_when(
      species  == 'Cel' ~ '*C. elegans*',
      species  == 'Bma' ~ '*B. malayi*',
      species  == 'Asu' & gene_ortholog == 'TYRA-3' ~ 'Outgroup',
      species  == 'Asu' ~ '*A. suum*',
      species  == 'Hco' ~ '*H. contortus*',
      species  == 'Hsa' ~ '*H. sapiens*',
      species  == 'Sra' ~ '*S. ratti*',
      species == 'Tmu' ~ '*T. muris*'))

tree_data.nr$gene_ortholog <- str_replace(tree_data.nr$gene_ortholog, "TYRA-3", "Outgroup")

# plot
phylo_plot<- ggtree(nwk.nr) %<+% tree_data.nr  + 
  geom_tiplab(aes(label=gene_ortholog),  hjust=-.2,size= 2.8, color='brown') +
  geom_tippoint(aes(colour=species), size = 2.8) +
  geom_nodelab(aes(label=label),size=2.5, hjust=1.5, vjust=-.8)+
  xlim(0,3.2)+
  scale_color_manual(values = c("#C2DDC7","#dd0d2c","#f4c33a","#E2B8A3", "#749092", "#72502E" , "#747B7D","#310A15"))+
  theme(legend.position= "bottom",legend.title=element_blank(), 
        legend.text = element_markdown(size=10),  )
phylo_plot

#---------------------------------set theme for plots-----------------------------------
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
    axis.title.x = element_text(angle = 0,size=12, face= "plain"))
#-----------------------------------qPCR-------------------------------------------------

#plot
qPCR_plot <- ggplot(qPCR_data, aes(x = target_name, y = dct)) +
  geom_errorbar(aes(color = "black"),
                position = position_dodge(width = 0.97),
                stat = "summary", width = 0.5, alpha = 0.9 ) +
  geom_bar(aes(color = "black", fill= "white"), position = position_dodge(width = 0.97), stat = "summary", fun = "mean") +
  geom_quasirandom(aes(color = target_name),
                   shape = 21,fill ="#f4c33a" , size = 2,
                   alpha = 0.75, dodge.width = 0.9, width = 0.1)+
  labs(x = "", y = expression(dC[T]), title = " ") +
  facet_grid(cols = vars(sample_name), scales= "free_y") +
  theme_minimal(base_size = 16, base_family = "Helvetica") +
  scale_fill_manual(values = c("#FFFFFF","#FFFFFF","#FFFFFF"))+
  scale_color_manual(values = c("#000000","#000000","#000000","#000000"))+
  aesthetic+
  theme(panel.border = element_rect(color="black", fill=NA))+
  NULL
qPCR_plot

#---------------------------------Whole Organisms--------------------------------------

#adult motility 
plot_data <- m.data %>%
  mutate(
    species = factor(species,
                     levels = c("Bpa", "Bma"),
                     labels = c("*B. pahangi*", "*B. malayi*")),
    smooth_time_point = case_when(
      time_point == 0 ~ 0,
      time_point == 0.1 ~ 1,
      time_point == 1 ~ 2,
      time_point == 24 ~ 3,
      time_point == 48 ~ 4),
    conc = factor(conc,
                  levels = c("0uM", "10uM", "100uM"),
                  labels = c("0 uM", "10 µM", "100 µM"))) %>%
    mutate(treatment = factor(treatment, levels = c("Water", "Acetylcholine", "Atropine", "Arecoline", "Carbachol", "Oxotremorine M", "Nicotine","Levamisole"),labels = c("CON", "ACH","ATR","ARE","CAR","OXO","NIC","LEV")))

plot_data<- plot_data%>%
  filter(species != "*B. malayi*")

#statistics 
stat_layer <- plot_data %>%
  # filter(time_point %in% c(24, 48)) %>%
  group_nest(species, smooth_time_point, stages, treatment, conc) %>%
  mutate(t = map(data, ~ t.test(log2(.x$self_norm_raw_flow), mu = 0))) %>%
  mutate(tidy = map(t, broom::tidy)) %>%
  unnest(tidy) %>%
  select(-t, -data) %>%
  # mutate(p.adjust = p.adjust(p.value, 'hochberg')) %>%
  mutate(
    sig = case_when(
      p.value <= 0.0001 ~ "***",
      p.value <= 0.001 ~ "***",
      p.value <= 0.01 ~ "**",
      p.value <= 0.05 ~ "*",
      p.value > 0.05 ~ "",))


stat_layer1<-stat_layer%>%
  filter(conc =="10 µM")%>%
  mutate(y= 1.5)

stat_layer2<-stat_layer%>%
  filter(conc =="100 µM")%>%
  mutate(y= 2)

#plot 
(adult_plot <- plot_data %>%
    ggplot(aes(x = as.numeric(smooth_time_point), y = log2(self_norm_raw_flow), color = conc, fill = conc)) +
    geom_line(aes(group = interaction(exp_date, stages, well)), alpha = 0.2, size = 0.3) +
    geom_hline(yintercept = log2(1), color = "black", alpha = 0.5, linetype = "dashed") +
    geom_smooth(aes(group = interaction(species, stages, conc)), size = 0.5, se = TRUE) +
    geom_point(alpha = 0.35, size = 1) +
    geom_text(
      data = stat_layer1, aes(label = sig, y = y), position = position_dodge(0.6),
      size = 3, show.legend = FALSE ) +
    geom_text(
      data = stat_layer2, aes(label = sig, y = y), position = position_dodge(0.6),
      size = 3, show.legend = FALSE ) +
    scale_x_continuous(labels = c("\u00d8", "0", "1", "24", "48")) +
    scale_color_manual(values = c("#F4C33A","#8E1F1F","#84D6CC")) +
    scale_fill_manual(values = c("#F4C33A","#8E1F1F","#84D6CC")) +
    scale_y_continuous(limits = c(-2, 2.5), breaks = seq(-2, 2, 1))+
    facet_grid(rows = vars(stages), cols = vars(treatment)) +
    labs(
      x = "Hours after drug addition", y = expression("Log"[2]*"FC movement"),
      color = "Concentration", fill = "Concentration") +
    aesthetic +
   theme(plot.margin = margin(0, 0, 0, 10, unit = "pt")))

#Adult female fecundity plot

f.data <- f.data %>%
  mutate(
    species = factor(species,
                     levels = c("Bpa", "Bma"),
                     labels = c("*B. pahangi*", "*B. malayi*")))

f.data <- f.data %>%
  mutate(
  treatment = factor(treatment, 
                     levels = c("Water", "Acetylcholine", "Atropine", "Arecoline", "Carbachol", "Oxotremorine M", "Nicotine","Levamisole"),
                     labels = c("CON", "ACH","ATR","ARE","CAR","OXO","NIC","LEV")))

#normalization
mean<- f.data%>%
  filter(treatment == "CON")%>%
  group_by(species,other) %>%
  summarize(mean_pixel = mean(pixel_count))

corrected.f<-left_join (f.data, mean) 

data<- corrected.f%>%          
  mutate(percent = (pixel_count / mean_pixel))

data$conc <- factor(data$conc , levels = c("0uM","10uM", "100uM"), labels = c("0µM","10µM","100µM"))

#stats 
sapply(lapply(data, unique), length)

data<-data %>%
  filter(other != "0hr")%>%
  filter(species != "*B. malayi*")

statsf<-data%>%
  filter(conc %in% c("0µM","10µM"))
         
statsf<-compare_means(percent~treatment , data = statsf, group.by= "other",
                         method = "t.test", ref.group = "CON")
statsf2<-data%>%
  filter(conc %in% c("0µM","100µM"))

statsf2<-compare_means(percent~treatment , data = statsf2, group.by= "other",
                      method = "t.test", ref.group = "CON")

statsf$p.signif<- factor(statsf$p.signif, levels=c("ns","*","**","***","****"),labels= c("","*","**","***","***"))

#plot
fecundity.plot <- data %>%
  ggplot( aes(x = treatment, y = percent, fill = conc)) +
  facet_grid(cols = vars(), rows= vars(other),scales= "free") +
  geom_errorbar(aes( color = conc),
                position = position_dodge(width = .8),
                stat = "summary", width = 0.5, alpha = 0.75) +
  geom_bar(aes(fill= conc),position = "dodge", stat = "summary", fun = "mean") +
  geom_quasirandom(aes(color = conc),
                   shape = 21, fill = "white", size = 2,
                   alpha = 0.75, dodge.width = 0.9, width = 0.1)+
  scale_y_continuous(labels = scales::percent,limits = c(0, 4), breaks = seq(0, 4, 1),expand = c(0, 0))+
  labs(title= "",x = "group", y = "% Fecundity") +
  scale_color_manual(values = c("#F4C33A","#8E1F1F","#84D6CC")) +
  scale_fill_manual(values = c("#F4C33A","#8E1F1F","#84D6CC")) +
  aesthetic+ 
NULL
fecundity.plot

#plot add stats and legend 
f.plot<-fecundity.plot+
  geom_text(data = tibble(percent = 3, treatment="LEV", conc= "10µM", other="24hr", label = '***'),
          aes(label = label),
          size = 4, color="#8E1F1F")+
  aesthetic +
  theme(legend.position= "bottom",
        legend.text=element_text(size=10),
        legend.title=element_blank(), 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
  
f.plot

#-------------------------------------microfilariae flow--------------------------------

mf.data$conc <- factor(mf.data$conc, levels = c("NA", "10mM", "1mM","100uM","10uM","1uM", "100nM","10nM","1nM"))
mf.data$treatment <- factor(mf.data$treatment, levels = c("CON", "SER","ACH","ATR","CAR", "OXO","NIC","LEV","HK"))

#normalize 

mf.data<- mf.data%>%
  mutate(normalized= (optical_flow)/(log2(worm_area)))

c_model <- mf.data %>%
  mutate(smooth_row = as.numeric(factor(row))) %>%
  mutate(smooth_col = as.numeric(col)) %>%
  # remove the controls, adjust for your data
  filter(str_detect(treatment, "CON")) %>%
  lm(normalized ~ smooth_row + smooth_col, data = .) %>%
  tidy()

row_coef <- as.numeric(c_model[2, 2])
col_coef <- as.numeric(c_model[3, 2])

corrected <-mf.data %>%
  mutate(smooth_row = as.numeric(factor(row))) %>%
  mutate(smooth_col = as.numeric(col)) %>%
  mutate(corrected_flow = normalized - smooth_row * row_coef - smooth_col * col_coef)

mean<- corrected%>%
  filter(treatment == "CON")%>%
  group_by(other, treatment) %>%
  summarize(mean_flow = mean(corrected_flow))%>%
  select(- treatment, mean_CON_flow = mean_flow)

corrected<-left_join (corrected, mean) 
corrected<- corrected%>%          
  mutate(percent = (corrected_flow / mean_CON_flow))

mean_percent<- corrected%>%
  group_by(other, treatment, conc,plate) %>%
  summarize(mean_percent = mean(percent))
#select(- treatment, mean_CON_flow = mean_flow)

# set levels and groups
corrected<- mean_percent%>% 
  mutate(
    conc= case_when(             
      conc == NA ~ '0',
      conc == '10mM' ~ '10mM',
      conc == '1mM' ~ '1mM',
      conc== '100uM' ~ '100uM',
      conc == '10uM' ~ '10uM',
      conc== '1uM' ~ '1uM',
      conc== '100nM' ~ '100nM',
      conc== '10nM' ~ '10nM',
      conc== '1nM' ~ '1nM',
      TRUE ~ as.character(conc)))

corrected$conc <- factor(corrected$conc , levels = c("1nM","10nM","100nM","1uM","10uM","100uM","1mM","10mM"), labels=c('-9','-8','-7','-6','-5','-4','-3','-2'))

corrected<-corrected%>%
  filter(other!= "72hr")

corrected<- corrected%>%
  filter(!grepl('p03', plate))

corrected<-corrected%>%
  filter(other!= "72hr")%>%
  filter(treatment!= "SER")

set1<-corrected %>%
  filter(!treatment %in% c("CON", "HK"))
set2<-corrected %>%
  filter(!treatment %in% c("ACH","ATR","CAR", "OXO","NIC","LEV"))

#stats 
stat_layermfm <-corrected %>%
  group_nest( other, treatment, conc) %>%
  mutate(t = map(data, ~ t.test(log2(.x$mean_percent), mu = 0))) %>%
  mutate(tidy = map(t, broom::tidy)) %>%
  unnest(tidy) %>%
  select(-t, -data) %>%
  # mutate(p.adjust = p.adjust(p.value, 'hochberg')) %>%
  mutate(
    sig = case_when(
      p.value <= 0.0001 ~ "***",
      p.value <= 0.001 ~ "***",
      p.value <= 0.01 ~ "**",
      p.value <= 0.05 ~ "*",
      p.value > 0.05 ~ "",
      
    ),
    y = .5
  )

#plot
(motility.plot <- corrected %>% 
    ggplot(aes(x=conc, y =log2(mean_percent),color="black", fill='black', group= plate)) +
    geom_line(color="black")+
    geom_text(
     data = stat_layermfm, aes( label = sig,x= as.numeric(conc), y = y),inherit.aes = FALSE, position = position_dodge(0.6),
     size = 3, show.legend = FALSE) +
    geom_hline(yintercept = 0, color = "black", alpha = 0.5, linetype = "dashed") +
    geom_point(data=set1,aes(color = "black", fill="black"),alpha = 0.35, size = 1) +
    geom_point(data=set2,aes(color = "black",fill="black"),alpha = 0.8, size = 1) +
    geom_smooth(aes(group = interaction(other, treatment)), color= "black", fill="black", size = 0.5, se = TRUE) +
    facet_grid(cols=vars(treatment),rows=vars(other), scales="free")+
    scale_y_continuous(limits = c(-3, 1), breaks = seq(-2, 1, 2))+
    labs(title= ,x = expression("Log"[10]* "[Concentration]"), y = expression("Log"[2]*"FC movement"))+
    aesthetic)


#celltox

ct.data$conc <- factor(ct.data$conc , levels = c("1nM","10nM","100nM","1uM","10uM","100uM","1mM","10mM"), labels=c('-9','-8','-7','-6','-5','-4','-3','-2'))
ct.data$treatment <- factor(ct.data$treatment, levels = c("CON", "SER","ACH","ATR","CAR", "OXO","NIC","LEV","HK"))

#normalization
ct.data<- ct.data%>%
  mutate(normalized= (1)/(AreaOccupied_AreaOccupied_GreenWorms))%>%
  filter(treatment !="SER")
  
ct.log<- ct.data%>%
  mutate(log = (log2(AreaOccupied_AreaOccupied_GreenWorms)))

ct.mean<- ct.log%>%
  filter(treatment == "CON")%>%
  group_by(other, treatment) %>%
  summarize(mean_flu = mean(log))%>%
  select(- treatment, mean_CON_flu = mean_flu)

ct.norm<-left_join (ct.data,ct.mean)%>%
  left_join(., ct.log)

ct.norm<- ct.norm%>%
  mutate(percent = (log / mean_CON_flu))

#stats 

statsmf2<-compare_means(AreaOccupied_AreaOccupied_GreenWorms~treatment, data = ct.norm,
                       method = "t.test", ref.group = "CON")

statsmf2$p.signif<- factor(statsmf2$p.signif, levels=c("ns","*","**","***","****"),labels= c("","*","**","***","***"))

#plot 
(ct.data.plot <- ct.norm %>% 
    ggplot(aes(conc, y = percent))+
    geom_errorbar(aes(color = "black"),
                  position = position_dodge(width = 0.97),
                  stat = "summary", width = 0.5, alpha = 0.75 ) +
    geom_bar(aes(fill = "grey"), position = position_dodge(width = 0.5), stat = "summary", fun = "mean") +
    geom_quasirandom(aes(color = "#f4c33a"),
                     shape = 21, fill = "white", size = 2,
                     alpha = 0.75, dodge.width = 0.9, width = 0.1)+
    facet_grid(cols=vars(treatment),rows=vars(other), scales="free")+
    scale_color_manual(values = c("#F06E73", "grey"))+
    scale_fill_manual(values = c("grey", "#F06E73"))+ 
   scale_y_continuous(labels = scales::percent,limits = c(0, 1.5), breaks = seq(0, 2, .25),expand = c(0, 0))+
    labs(title= ,x = expression("Log"[10]* "[Concentration]"), y = "% FLU") +
    aesthetic+
    theme(strip.background = element_blank(),
                   strip.text.x = element_blank()))


#-----------------------------------------figure plot-----------------------------------

fig1_row1<- plot_grid( phylo_plot, qPCR_plot, labels = c('A','B'), label_size = 16, ncol=2, rel_widths= c(1.2,1))
fig1_row1

fig1_adult <- plot_grid(adult_plot, f.plot, labels = c('C','D'), label_size = 16, ncol=1, rel_widths= c(1,1), rel_heights= c(1,1.2), align = "vh", axis="rl")
fig1_adult

illust.hb2 <- image_read_pdf("../images/fig1_mf_flow.pdf", density = 100)

Figill <- ggdraw() +
  draw_image(illust.hb2, scale = 1)

fig1_mf<-plot_grid(motility.plot,Figill,ct.data.plot,labels = c('E','','F'), label_size = 16, ncol=1, rel_widths= c(1,1,1),rel_heights= c(1.8,1,1.8), align = "vh", axis="rl")

fig1<-plot_grid(fig1_row1, fig1_adult, fig1_mf, labels = c('','',''), nrow=3, rel_heights= c(0.9,1.5,2),align = "vh", axis="rl")

save_plot("../plots/fig1_boot.pdf", fig1, base_height = 16, base_width = 12)

