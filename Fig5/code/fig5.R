# other plotting
library(ggplot2)
library(tidyr)
library(shades)
library(ggbeeswarm)
library(cowplot)

# data wrangling
library(janitor)
library(tidyverse)
library(dplyr)

#PDF import
library(magick)
library(pdftools)
library(grImport2)

# misc
library(magrittr)
library(ggtext)
library(ggdendro)
library(readr)
library(ggpubr)
library(stringr)
library(lubridate)
library(conflicted)
library(viridis)
library(here)
library(ZamanianLabThemes)


conflict_prefer("filter", "dplyr")
#------------------------------------set working directory-----------------------------
setwd<- here()

#-------------------------------------load in EPG data--------------------------------------
#Update RDS file with new data as needed

EPG.data <- readRDS("../data/EPG.rds")

data<- EPG.data%>%
  select(1:11,14:17, 20, 21, 10)
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
    panel.border = element_rect(colour= "black", fill=NA))

#-------------------------------------EPG validation 

figB.data<- EPG.data %>%
  filter(date < "2020-07-30") 

figB.data$treatment <- factor(figB.data$treatment, levels = c("Untreated","Serotonin","Serotonin_Arecoline","Arecoline","Serotonin_Atropine","Serotonin_Nicotine","Serotonin_Levamisole","Serotonin_CaEPGchol","Serotonin_Acetylcholine","Serotonin_Octopamine","Serotonin_Oxotremorine"))
figB.data$strain <- factor(figB.data$strain, levels = c("N2","*gar-3(gk305)*","*gar-3(gk305);<br>myo-3p::Bma-gar-3*","*gar-3(gk305);<br>myo-2p::Bma-gar-3*"))

#stats 

stats<-compare_means(Number.of.Pumps~treatment , data = figB.data, group.by= "strain",
                      method = "t.test", ref.group = "Untreated")
stats$p.signif<- factor(stats$p.signif, levels=c("ns","*","**","***","****"),labels= c("","*","**","***","***"))

stats2<-compare_means(Number.of.Pumps~treatment , data = figB.data, group.by= "strain",
                     method = "t.test", ref.group = "Serotonin")
stats2$p.signif<- factor(stats2$p.signif, levels=c("ns","*","**","***","****"),labels= c("","*","**","***","***"))

#plot
fig5B.plot <- ggplot(filter(figB.data, date > "2019-08-29"), aes( treatment, y = Number.of.Pumps / Recording.Duration * 60, color=strain)) +
  facet_grid(cols= vars(strain)) +
  geom_errorbar(aes(color = strain),
                position = position_dodge(width = 0.97),
                stat = "summary", width = 0.5, alpha = 0.75 ) +
  geom_boxplot(aes(color=strain),width=0.8, alpha=2) +
  stat_pvalue_manual(
  stats, x="group2", y.position = 315,
  label = "p.signif",
 )+
  stat_pvalue_manual(
    stats2, x="group2", y.position = 340,
    label = "p.signif", color= "brown",
    )+
  scale_x_discrete(labels= c("Untreated"="Untreated","Serotonin"= "5HT","Arecoline"="ARE","Serotonin_Arecoline"="5HT+ARE"))+
  scale_color_manual(values = c("black", "gray", "#005900","#050298"))+
  scale_fill_manual(values =  c("#F0F0F0", "#F0F0F0", "#F0F0F0","#F0F0F0"))+
  scale_y_continuous( limits = c(0, 400), breaks = seq(0, 350, 50),expand = c(0, 0))+
  labs(x = "", y = "Pumps / Minute") +
  aesthetic+
  geom_beeswarm(aes(color=strain, fill= "grey"),stroke = .2) 
  

fig5B.plot

#----------------------------------cholinergic panel------------------------------------

fig5c<-filter(data, date > "2020-07-30")

#normalization
fig5c_oxo<- group_by(fig5c, date, strain, treatment) %>%
  mutate(mean_pumps = Mean.Frequency)
#get mean pump of each treatment for each strain, new table, new column 

fig5c_control_summary <- filter(fig5c_oxo, treatment == "Untreated") %>%
  group_by(date, strain, treatment) %>% 
  summarize(mean_pumps = mean(Mean.Frequency)) %>%
  select(-treatment, mean_control_pumps = mean_pumps)

fig5c_5HT_summary <- filter(fig5c_oxo, treatment == "Serotonin") %>%
  group_by(date, strain, treatment) %>% 
  summarize(mean_pumps = mean(Mean.Frequency)) %>%
  select(-treatment, mean_5HT_pumps = mean_pumps)
# sets mean control pumps for strain= mean of treatment 

fig5c <- left_join(fig5c, fig5c_control_summary) %>%
  left_join(., fig5c_5HT_summary) %>%
  mutate(normalized_pumps = (Mean.Frequency - mean_control_pumps)/(mean_5HT_pumps - mean_control_pumps))

fig5c<- fig5c%>%          
  mutate(percent = (Mean.Frequency / mean_5HT_pumps))
# divide pumps by set treatment refrence (control) pumps

#set order and names
fig5c$treatment <- factor(fig5c$treatment, levels = c("Untreated","Serotonin","Serotonin_Acetylcholine","Serotonin_Atropine","Serotonin_Carbachol","Serotonin_Oxotremorine","Serotonin_Octopamine","Serotonin_Nicotine","Serotonin_Levamisole"),
                          labels = c("Untreated","5HT","5HT+ACH","5HT+ATR","5HT+CAR","5HT+OXO","5HT+OCT","5HT+NIC","5HT+LEV"))

fig5c<-fig5c%>%
  filter(treatment!= "5HT+OCT")

#stats
stats3<-compare_means(percent~treatment , data = fig5c, group.by= "strain",
                     method = "t.test", ref.group = "5HT")
stats3$p.signif<- factor(stats3$p.signif, levels=c("ns","*","**","***","****"),labels= c("","*","**","***","***"))

#plot

fig5C.plot <- ggplot(filter (fig5c), aes(x = treatment, y = percent)) +
  geom_errorbar(aes(color = strain),
                position = position_dodge(width = 0.97),
                stat = "summary", width = 0.5, alpha = 0.75 ) +
  geom_bar(aes(fill = strain, color= strain), position = position_dodge(width = 0.5), stat = "summary", fun = "mean") +
  geom_quasirandom(aes(color = strain),
                   shape = 21, fill = "white", size = 2,
                   alpha = 0.75, dodge.width = 0.9, width = 0.1)+
  #geom_text(aes(label = format( (..y..*100), digits= 2)),stat = "summary",vjust= -1.5, colour = "Black")+
  facet_grid(cols = vars(strain)) +
  labs(title= "",x = "", y = "% 5HT pumps") +
  theme_minimal(base_size = 10, base_family = "Helvetica") +
  stat_pvalue_manual(
    stats3, x="group2", y.position = 2.25,
    label = "p.signif", size=3,
    position = position_dodge(0.8))+
  scale_fill_manual(values =  c("#F0F0F0", "#F0F0F0", "#F0F0F0","#F0F0F0"))+
  scale_color_manual(values = c("black", "gray", "#005900","#050298"))+
  scale_y_continuous(labels = scales::percent,limits = c(0, 3), breaks = seq(0, 2.50, .5),expand = c(0, 0))+
  aesthetic

fig5C.plot


#--------------------------------Heat map---------------------------------------------- 
#data refresh
EPG.data <- readRDS("../data/EPG.rds")

recent.data <- EPG.data %>%
  filter(date > "2020-07-30") 

long.data<- recent.data %>%
  pivot_longer(cols=c( 13:25), names_to = "metric", values_to = "value")%>%
  mutate(value = ifelse(is.infinite(value), 0,value))

#normalize
summary <-long.data%>%
  filter(value != "na")%>%
  group_by(strain, treatment, metric)%>%
  summarise(mean_phenotype = mean(value))%>%
  mutate(mean_phenotype= ifelse(is.infinite(mean_phenotype), 0,mean_phenotype))

N2.summary <-summary%>%
  #filter(metric %in% c('Mean.Amplitude','Mean.Frequency','Mean.IPI.Duration','Mean.Pump.Duration'))%>%
  filter( treatment== 'Serotonin', strain == 'N2')%>%
  group_by( metric)%>%
  mutate(norm = mean_phenotype)%>%
  select(-treatment,-mean_phenotype, -strain, norm_5HT = norm)

data<-left_join(N2.summary, summary)%>%
  mutate(normalized = (mean_phenotype / norm_5HT))

# scale normalization (prepare matrix, scale / log normalize)
data.m <- data %>% select(strain, treatment,metric,mean_phenotype) %>%
  dplyr::mutate(strain_treatment = paste0(strain,"%",treatment)) %>%
  select(strain_treatment,metric,mean_phenotype) %>%
  pivot_wider(names_from = metric, values_from = mean_phenotype) %>%
  column_to_rownames(var = "strain_treatment")
data.m <- data.matrix(data.m, rownames.force = TRUE)
data.m <- t(data.m)
data.m <- t(scale(t(log(data.m+1)),center=TRUE,scale=TRUE))

# calculate phenotype distances and cluster
pDists <- dist(data.m, method = "euclidean") 
pclust.dist <- hclust(pDists, method="ward.D2")
ord <- pclust.dist$order #if you want to arrange phenotypes by distance

#convert count matrix to df, re-order based on clustering, and tidy (long form + metadata) for plotting
data <- as.data.frame(data.m) %>%
  rownames_to_column(var = "metric")
data$metric <- factor(data$metric, levels = c(data$metric[ord]))
data <- data %>%
  pivot_longer(2:ncol(data), names_to = "strain_treatment", values_to = "mean_phenotype") %>%
  separate(strain_treatment, c("strain","treatment"), sep = "%", extra = "warn")

#set  order and names 
data$strain <- factor(data$strain, levels = c("N2","*gar-3(gk305)*","*gar-3(gk305);<br>myo-3p::Bma-gar-3*","*gar-3(gk305);<br>myo-2p::Bma-gar-3*"))
data$treatment <- factor(data$treatment, levels = c("Untreated","Serotonin","Serotonin_Acetylcholine","Serotonin_Atropine","Serotonin_Carbachol","Serotonin_Oxotremorine","Serotonin_Octopamine","Serotonin_Nicotine","Serotonin_Levamisole"),
                         labels = c("Untreated","5HT","5HT+ACH","5HT+ATR","5HT+CAR","5HT+OXO","5HT+OCT","5HT+NIC","5HT+LEV"))

#plot
heatmap <- ggplot(filter(data, treatment != "5HT+OCT"),scale="column",aes( x= treatment, y= metric, fill=mean_phenotype))+
  geom_tile() + 
  facet_grid(cols = vars(strain)) +
  scale_fill_gradient2(low= "#4D6291", mid= "white", high = "#FDA50F", limits=c(-3.5,3.5))+
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) + 
  #scale_fill_viridis(discrete=FALSE)+
  labs(title= "",x = "", y = "", fill = "Z-score") +
  aesthetic+
  geom_hline(yintercept = c(2.5,4.5,6.5,8.5,11.5), color = "#3E3D53", alpha = 1) +
  theme(
    panel.grid.major.y = element_line(),
    axis.text.x = element_markdown(angle = 45, hjust = 1, face="plain",size=10),
    axis.text.y = element_markdown(face="plain",size=9),
    legend.position="bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10))
heatmap


# generate dendrogram
dend <- ggdendrogram(pclust.dist, rotate = 90, size = 3) + 
  theme(axis.text.x = element_blank(),axis.text.y = element_blank())
dend

#------------------------------------Amplitude phenotype------------------------------------

amplitude<- EPG.data%>%
  select(1:11,14:17, 20, 21, 10)%>%
  mutate(Mean.Amplitude = ifelse(is.na(Mean.Amplitude), 0,Mean.Amplitude))%>%
  filter(date > "2020-07-30")

#normalize
fig5E5ht_summary <- filter(amplitude, treatment == "Serotonin") %>%
  group_by(date, strain, treatment) %>% 
  summarize(mean_amp = mean (Mean.Amplitude)) %>%
  select(-treatment, mean_5ht_amp = mean_amp)

fig5E <- left_join(amplitude, fig5E5ht_summary) %>%
  mutate(norm = ((Mean.Amplitude) / (mean_5ht_amp)))

#set order and names
fig5E$treatment <- factor(fig5E$treatment, levels = c("Untreated","Serotonin","Serotonin_Acetylcholine","Serotonin_Atropine","Serotonin_Carbachol","Serotonin_Oxotremorine","Serotonin_Octopamine","Serotonin_Nicotine","Serotonin_Levamisole"),
                          labels = c("Untreated","5HT","5HT+ACH","5HT+ATR","5HT+CAR","5HT+OXO","5HT+OCT","5HT+NIC","5HT+LEV"))

fig5E<-fig5E%>%
  filter(treatment!= "5HT+OCT")

#stats
stats4<-compare_means(norm~treatment , data = fig5E, group.by= "strain",
                      method = "t.test", p.adjust.method = "none", ref.group = "5HT")
stats4$p.signif<- factor(stats4$p.signif, levels=c("ns","*","**","***","****"),labels= c("","*","**","***","***"))

#plot
fig5E.plot <- ggplot(filter(fig5E), aes(treatment, y = norm, color=strain)) +
  facet_grid(cols=vars(strain), scales="free") +
  geom_errorbar(aes(color = strain),
                position = position_dodge(width = 0.97),
                stat = "summary", width = 0.5, alpha = 0.75 ) +
  geom_bar(aes(fill = strain, color= strain), position = position_dodge(width = 0.5), stat = "summary", fun = "mean") +
  labs(x = "", y = "% 5HT amplitude") +
  geom_hline(yintercept = 1)+
  geom_quasirandom(aes(color = strain),
                   shape = 21, fill = "white", size = 2,
                   alpha = 0.75, dodge.width = 0.9, width = 0.1)+
  stat_pvalue_manual(
    stats4, x="group2", y.position = 2.75,
    label = "p.signif", size=3,
    position = position_dodge(0.8))+
  scale_x_discrete(labels=c("untreated"= "Untreated", "Serotonin"= "5HT", "Arecoline"="ARE","Serotonin_Arecoline"="5HT+ARE","Serotonin_Atropine"="5HT+ATR","Serotonin_Nicotine"="5HT+NIC","Serotonin_Levamisole"="5HT+LEV","Serotonin_Carbachol"="5HT+CAR","Serotonin_Acetylcholine"="5HT+ACH","Serotonin_Octopamine"="5HT+OCT","Serotonin_Oxotremorine"="5HT+OXO"))+
  theme_minimal(base_size = 10, base_family = "Helvetica") +
  scale_fill_manual(values =  c("#F0F0F0", "#F0F0F0", "#F0F0F0","#F0F0F0"))+
  scale_color_manual(values = c("black", "gray", "#005900","#050298"))+
  scale_y_continuous(labels = scales::percent,limits = c(0, 3), breaks = seq(0, 2.50, .5),expand = c(0, 0))+
  aesthetic+
  theme()

fig5E.plot

#----------------------------------------------figure 5 ---------------------------------

illust.hb <- image_read_pdf("../images/EPG_Schematic.pdf", density = 600)
fig5A <- ggdraw() +
  draw_image(illust.hb, scale = 1)

fig5BC <- plot_grid(fig5B.plot, fig5C.plot, labels = c('B','C'), label_size = 16, nrow=2, align = "vh", axis="lr", rel_heights= c(1,1))

fig5D<- plot_grid(heatmap + theme(plot.margin = unit(c(0, 0, 0, -0.25), "cm")), dend + theme(plot.margin = unit(c(1.75, 0.2, 3.35, 0), "cm")),labels = c('D', ''), label_size = 16, nrow=1, rel_widths= c(1,.075))

#fig5E <- plot_grid(fig5E.plot + theme(plot.margin=margin(0,0,0,0,unit="cm")), fig5Enic.plot + theme(plot.margin=margin(0,0,0,0,unit="cm")), labels = c('E',''), label_size = 16, nrow=1, align = "vh", rel_widths= c(1,.6))

fig5BCDE <- plot_grid(fig5BC,fig5D,fig5E.plot, labels= c('','','E'), nrow=3, rel_heights= c(1.2,0.8,0.55))

fig5<-plot_grid(fig5A, fig5BCDE, nrow=2, labels= c('A',''), label_size = 16, rel_heights= c(0.7,3))

save_plot("../figures/plots/fig5.pdf", fig5, base_height = 16, base_width = 12)


#------------------------------------------Supplementals------------------------------------ 

#-------------------------------EPG sup of all metrics 

long.plot.data<-long.data%>%
  filter(date > "2020-07-30")%>%
  filter(metric %in% c("Mean.Pump.Duration","Mean.IPI.Duration","Pump.Duration.SD","Mean.Amplitude","Median.R.to.E.Ratio", "Mean.Frequency"))%>%
  mutate(value = ifelse(is.na(value), 0,value))
  
long.plot.data$treatment <- factor(long.plot.data$treatment, levels = c("Untreated","Serotonin","Serotonin_Acetylcholine","Serotonin_Atropine","Serotonin_Carbachol","Serotonin_Oxotremorine","Serotonin_Nicotine","Serotonin_Levamisole"),
                         labels = c("Untreated","5HT","5HT+ACH","5HT+ATR","5HT+CAR","5HT+OXO","5HT+NIC","5HT+LEV"))

#plot
fig5S.plot <- ggplot(filter(long.plot.data, treatment != "Serotonin_Octopamine"), aes(treatment, y = value, color=strain)) +
  facet_grid(cols=vars(strain),rows=vars(metric), scales="free") +
  geom_errorbar(aes(color = strain),
                position = position_dodge(width = 0.97),
                stat = "summary", width = 0.5, alpha = 0.75 ) +
   geom_boxplot(aes(color=strain),width=0.8, alpha=2) +
  labs(x = "", y = "") +
  geom_quasirandom(aes(color = strain),
                   shape = 21, fill = "white", size = 2,
                   alpha = 0.75, dodge.width = 0.9, width = 0.1)+
  theme_minimal(base_size = 10, base_family = "Helvetica") +
  scale_fill_manual(values =  c("#F0F0F0", "#F0F0F0", "#F0F0F0","#F0F0F0"))+
  scale_color_manual(values = c("black", "gray", "#005900","#050298"))+
  scale_y_continuous(expand = c(0, 0))+
  aesthetic
# amp uV, Duration ms, IPI ms, freq hz 

fig5S.plot

save_plot("../sup/fig5sup.pdf", fig5S.plot, base_height = 16, base_width = 12)
