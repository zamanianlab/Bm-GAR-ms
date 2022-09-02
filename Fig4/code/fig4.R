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

#set working directory 
setwd<-here()

#-------------------------------load in data---------------------------------------------

MC.data <- readRDS("../data/MC.rds") 
#---------------------------------set theme----------------------------------------------
aesthetic <- 
  theme_zlab() +
  theme(
    panel.grid.major.y = element_line(),
    axis.text.x = element_markdown(angle = 45, hjust = 1, face="plain",size=10),
    axis.text.y = element_text(face = "plain", size = 10),
    strip.text.x = element_markdown(size=10),
    strip.text.y = element_markdown(size=12),
    strip.text = ggplot2::element_text(face = "plain", size = 10),
    axis.title.y = element_text(size=12, face= "plain"))

#------------------------------------Rename strains---------------------------------------

comparisons <-list(c( "N2", "*gar-3(gk305)*"),
                   c("N2", "*gar-3(gk305);<br>myo-3p::Bma-gar-3*"), 
                   c("N2","*gar-3(gk305);<br>myo-2p::Bma-gar-3*"),
                   c("*gar-3(gk305)*","*gar-3(gk305);<br>myo-3p::Bma-gar-3*"), 
                   c("*gar-3(gk305)*","*gar-3(gk305);<br>myo-2p::Bma-gar-3*"),
                   c("*gar-3(gk305);<br>myo-2p::Bma-gar-3*","*gar-3(gk305);<br>myo-3p::Bma-gar-3*"))

date_fig4A <- list(c("2019-07-25"),c("2019-08-06"),c("2019-08-13"))
comparisons <- list(c("Untreated", "Serotonin"),  c("Untreated", "5HT + Arecoline"),c("Untreated", "Arecoline"), c("Serotonin", "5HT + Arecoline"), c("Serotonin", "Arecoline"), c("Arecoline", "5HT + Arecoline"))

new_strip <- c("2019-07-25" = "L1 +OP50", "2019-08-06" = "Adult +OP50", "2019-08-13" = "Adult -OP50")


#-------------------------------plot 4 A
#visual count optimization

#stats
stats.fil<-filter(MC.data, date %in% date_fig4A)

stats<-compare_means(pps~treatment , data = stats.fil, group.by= c("strain","date"),
                      method = "t.test", ref.group = "Untreated")

stats$p.signif<- factor(stats$p.signif, levels=c("ns","*","**","***","****"),labels= c("","*","**","***","***"))

stats2<-compare_means(pps~treatment , data = stats.fil, group.by= c("strain","date"),
                     method = "t.test", ref.group = "Serotonin")

stats2$p.signif<- factor(stats2$p.signif, levels=c("ns","*","**","***","****"),labels= c("","*","**","***","***"))

#plot 
fig4A.plot <- ggplot(filter(MC.data, date %in% date_fig4A), aes(x = treatment, y = pps*60)) +
  geom_boxplot(aes(color=strain),width=0.8, alpha=2) +
  geom_beeswarm(aes(color=strain, fill= "grey"),stroke = .2) +
  stat_pvalue_manual(
    stats, x="group2", y.position = 420,
    label = "p.signif",
    position = position_dodge(0.8))+
  stat_pvalue_manual(
    stats2, x="group2", y.position = 480,
    label = "p.signif",color="brown",
    position = position_dodge(0.8))+
  labs(title= "",x = "", y = "Pumps / Minute") +
  scale_x_discrete(labels=c("untreated"= "Untreated", "Serotonin"= "5HT", "Arecoline"="ARE","5HT + Arecoline"="5HT+ARE" ))+
  theme_minimal(base_size = 16, base_family = "Helvetica") +
  scale_color_manual(values = c("black", "gray", "#005900","#050298"))+
  aesthetic +
  theme(panel.border = element_rect(color="black", fill=NA))+
  facet_grid(cols= vars(strain), rows = vars(date), labeller = labeller(date = new_strip), shrink = TRUE)+
  scale_y_continuous(limits = c(-5, 600), breaks = seq(0, 450, 100))


fig4A.plot

#-------------------------------fig4E 
#L1 cholinergic panel 

#normalization
fig4E_norm<-filter(MC.data, date %in% c("2020-07-15","2020-07-17"))

fig4E_summary <- group_by(MC.data, date, strain, treatment) %>%
  summarise(mean_pumps = mean(pumps))
#get mean pump of each treatment for each strain, new table, new column 

fig4E_control_summary <- filter(fig4E_summary , treatment == "Untreated") %>%
  select(-treatment, mean_control_pumps = mean_pumps)

fig4E_5HT_summary <- filter(fig4E_summary , treatment == "Serotonin") %>%
  select(-treatment, mean_5HT_pumps = mean_pumps)
# sets mean control pumps for strain= mean of treatment 

fig4E<- left_join(fig4E_norm,fig4E_control_summary) %>%
  left_join(., fig4E_5HT_summary) %>%
  mutate(normalized_pumps = (pumps - mean_control_pumps)/(mean_5HT_pumps - mean_control_pumps))
# divide pumps by set treatment refrence (control) pumps

fig4E$treatment <- factor(fig4E$treatment , levels = c("Untreated","Serotonin","Acetylcholine","Atropine","Carbachol","Levamisole"))

fig4E<- fig4E%>%          
  mutate(percent = (pumps / mean_control_pumps))

#stats
stats3<-compare_means(percent~treatment , data = fig4E, group.by= "strain",
                      method = "t.test", ref.group = "Untreated")

stats3$p.signif<- factor(stats3$p.signif, levels=c("ns","*","**","***","****"),labels= c("","*","**","***","***"))

#plot
fig4E.plot <- ggplot(filter(fig4E, treatment != "Nicotine"), aes(x = treatment, y = percent)) +
  geom_errorbar(aes(color = strain),position = position_dodge(width = 0.97),stat = "summary", width = 0.5, alpha = 0.75 ) +
  geom_bar(aes(color= strain, fill = strain), stat = "summary", fun = "mean") +
  geom_quasirandom(aes(color = strain),
                   shape = 21, fill = "white", size = 2,
                   alpha = 0.75, dodge.width = 0.9, width = 0.1)+
  facet_grid(cols = vars(strain)) +
  scale_x_discrete(labels= c("Serotonin"= "5HT","Atropine"="ATR","Nicotine"="NIC","Levamisole"="LEV","Carbachol"="CAR","Acetylcholine"="ACH","Octopamine"="OCT"))+
  labs(title= "L1",x = "", y = "% Control pumps") +
  stat_pvalue_manual(
    stats3, x="group2", y.position = 1.8,
    label = "p.signif")+
  theme_minimal(base_size = 16, base_family = "Helvetica") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 2), breaks = seq(0, 1.50, .5),expand = c(0, 0))+
  scale_fill_manual(values =  c("#F0F0F0", "#F0F0F0", "#F0F0F0","#F0F0F0"))+
  scale_color_manual(values = c("black", "gray", "#005900","#050298"))+
  aesthetic+
  theme(panel.border = element_rect(color="black", fill=NA))+
  NULL
fig4E.plot


#------------------------------------fig4F
#Adult cholinergic panel 

#filter dates
fig4F_norm<- filter(MC.data, date %in% c("2019-10-23", "2019-11-04"))

#normalize
fig4F_summary <- group_by(fig4F_norm, date, strain, treatment) %>%
  summarise(mean_pumps = mean(pumps))

#get mean pump of each treatment for each strain, new table, new column 

fig4F_5HT_summary <- filter(fig4F_summary, treatment == "Serotonin") %>%
  select(-treatment, mean_5HT_pumps = mean_pumps)
# sets mean control pumps for strain= mean of treatment 

fig4F <- left_join(fig4F_norm, fig4F_5HT_summary) %>%
  mutate(normalized_pumps = pumps / mean_5HT_pumps)
# divide pumps by set treatment refrence (control) pumps

fig4F$treatment <- factor (fig4F$treatment, levels = c("Untreated","Serotonin","Acetylcholine","Carbachol","Atropine","Levamisole"))

fig4F<- fig4F%>%          
  mutate(percent = (pumps / mean_5HT_pumps))

#stats
stats4<-compare_means(percent~treatment , data = fig4F, group.by= "strain",
                      method = "t.test", ref.group = "Serotonin")

stats4$p.signif<- factor(stats4$p.signif, levels=c("ns","*","**","***","****"),labels= c("","*","**","***","***"))

#plot
fig4F.plot <- ggplot(filter(fig4F,treatment != "Nicotine"), aes(x = treatment, y = percent)) +
  geom_errorbar(aes(color = strain),
                position = position_dodge(width = 0.97),
                stat = "summary", width = 0.5, alpha = 0.75 ) +
  geom_bar(aes(fill = strain, color= strain), position = position_dodge(width = 0.97), stat = "summary", fun = "mean") +
  geom_quasirandom(aes(color = strain),
                   shape = 21, fill = "white", size = 2,
                   alpha = 0.75, dodge.width = 0.9, width = 0.1)+
  facet_grid(cols = vars(strain)) +
  scale_x_discrete(labels= c("Untreated"="Untreated","Serotonin"= "5HT","Atropine"="5HT+ATR","Nicotine"="5HT+NIC","Levamisole"="5HT+LEV","Carbachol"="5HT+CAR","Acetylcholine"="5HT+ACH","Octopamine"="5HT+OCT"))+
  labs(title= "Adult",x = "", y = "% 5HT pumps") +
  stat_pvalue_manual(
    stats4, x="group2", y.position = 1.8,
    label = "p.signif")+
  scale_y_continuous(labels = scales::percent,limits = c(0, 2), breaks = seq(0, 1.50, .5),expand = c(0, 0))+
  theme_minimal(base_size = 16, base_family = "Helvetica") +
  scale_fill_manual(values =  c("#F0F0F0", "#F0F0F0", "#F0F0F0","#F0F0F0"))+
  scale_color_manual(values = c("black", "gray", "#005900","#050298"))+
  aesthetic+
  theme(panel.border = element_rect(color="black", fill=NA))+
  NULL

fig4F.plot


#------------------------------------- fig4c
#L1 OXO
#filter dates
fig4C_norm<-filter(MC.data, date %in% c("2020-02-26","2020-06-30","2020-07-01"))

#normalize
fig4C_summary <- group_by(fig4C_norm, date, strain, treatment) %>%
  summarise(mean_pumps = mean(pumps))
#get mean pump of each treatment for each strain, new table, new column 

fig4C_control_summary <- filter(fig4C_summary, treatment == "Untreated") %>%
  select(-treatment, mean_control_pumps = mean_pumps)

fig4C_5HT_summary <- filter(fig4C_summary, treatment == "Serotonin") %>%
  select(-treatment, mean_5HT_pumps = mean_pumps)
# sets mean control pumps for strain= mean of treatment 

fig4C <- left_join(fig4C_norm, fig4C_control_summary) %>%
  left_join(., fig4C_5HT_summary ) %>%
  mutate(normalized_pumps = (pumps - mean_control_pumps)/(mean_5HT_pumps - mean_control_pumps))

# divide pumps by set treatment refrence (control) pumps

fig4C<- fig4C%>%          
  mutate(percent = (pumps / mean_control_pumps))

fig4C$treatment <- factor(fig4C$treatment , levels = c("Untreated","Serotonin","Oxotremorine","Nicotine"))

#stats
stats5<-compare_means(percent~treatment , data = fig4C, group.by= "strain",
                      method = "t.test", ref.group = "Untreated")

stats5$p.signif<- factor(stats5$p.signif, levels=c("ns","*","**","***","****"),labels= c("","*","**","***","***"))

#plot
fig4C.plot <- ggplot(filter(fig4C), aes(x = treatment, y = percent)) +
  geom_errorbar(aes(color = strain),position = position_dodge(width = 0.97),
           stat = "summary", width = 0.5, alpha = 0.75 ) +
  geom_bar(aes(fill = strain, color= strain), position = position_dodge(width = 0.97), 
           stat = "summary", fun = "mean") +
  geom_quasirandom(aes(color = strain),
                   shape = 21, fill = "white", size = 2,
                   alpha = 0.75, dodge.width = 0.9, width = 0.1)+
  facet_grid(cols = vars(strain)) +
  scale_x_discrete(labels= c("Untreated"="Untreated","Serotonin"= "5HT","Oxotremorine"="OXO","Nicotine"="NIC"))+
  labs(title= "L1",x = "", y = "% Control pumps") +
  stat_pvalue_manual(
    stats5, x="group2", y.position = 1.8,
    label = "p.signif")+
  theme_minimal(base_size = 10, base_family = "Helvetica") +
  scale_fill_manual(values =  c("#F0F0F0", "#F0F0F0", "#F0F0F0","#F0F0F0"))+
  scale_color_manual(values = c("black", "gray", "#005900","#050298"))+
  scale_y_continuous(labels = scales::percent, limits = c(0, 2), breaks = seq(0, 1.50, .5), expand  = c(0,0))+
  aesthetic+
  theme(panel.border = element_rect(color="black", fill=NA))+
  
  NULL
fig4C.plot 

#----------------------------------fig4d
#adult OXO

#filter dates
fig4d_norm<-filter(MC.data, date %in% c("2019-11-26", "2019-12-06", "2019-12-27"))

#normalize
fig4d_summary <- group_by(fig4d_norm, date, strain, treatment) %>%
  summarise(mean_pumps = mean(pumps))
#get mean pump of each treatment for each strain, new table, new column 

fig4d_control_summary <- filter(fig4d_summary, treatment == "Untreated") %>%
  select(-treatment, mean_control_pumps = mean_pumps)

fig4d_5HT_summary <- filter(fig4d_summary, treatment == "Serotonin") %>%
  select(-treatment, mean_5HT_pumps = mean_pumps)
# sets mean control pumps for strain= mean of treatment 

fig4D <- left_join(fig4d_norm, fig4d_control_summary) %>%
  left_join(., fig4d_5HT_summary ) %>%
  mutate(normalized_pumps = (pumps - mean_control_pumps)/(mean_5HT_pumps - mean_control_pumps))
# divide pumps by set treatment refrence (control) pumps

fig4D<- fig4D%>%          
  mutate(percent = (pumps / mean_5HT_pumps))

fig4D$treatment <- factor(fig4D$treatment , levels = c("Untreated","Serotonin","Oxotremorine","Nicotine"))

#stats
stats6<-compare_means(percent~treatment , data = fig4D, group.by= "strain",
                      method = "t.test", ref.group = "Serotonin")

stats6$p.signif<- factor(stats6$p.signif, levels=c("ns","*","**","***","****"),labels= c("","*","**","***","***"))

#plot
fig4D.plot<- ggplot(fig4D, aes(x = treatment, y = percent)) +
  geom_errorbar(aes(color = strain),
                position = position_dodge(width = 0.97),
                stat = "summary", width = 0.5, alpha = 0.75 ) +
  geom_bar(aes(fill = strain, color= strain), position = position_dodge(width = 0.97), stat = "summary", fun = "mean") +
  geom_quasirandom(aes(color = strain),
                   shape = 21, fill = "white", size = 2,
                   alpha = 0.75, dodge.width = 0.9, width = 0.1)+
  facet_grid(cols = vars(strain)) +
 scale_x_discrete(labels= c("Untreated"="Untreated","Serotonin"= "5HT","Nicotine"="5HT+NIC","Oxotremorine"="5HT+OXO"))+
 labs(title= "Adult",x = "", y = "% 5HT pumps") +
  stat_pvalue_manual(
    stats6, x="group2", y.position = 1.8,
    label = "p.signif")+
  theme_minimal(base_size = 10, base_family = "Helvetica") +
  scale_fill_manual(values =  c("#F0F0F0", "#F0F0F0", "#F0F0F0","#F0F0F0"))+
  scale_color_manual(values = c("black", "gray", "#005900","#050298"))+
  scale_y_continuous(labels = scales::percent, limits = c(0, 2), breaks = seq(0, 1.50, .5), expand= c(0,0))+
 aesthetic+
  theme(panel.border = element_rect(color="black", fill=NA))+
  
  NULL
fig4D.plot

#------------------------------------Figure 4--------------------------------------------

illust.hb <- image_read_pdf("../images/MC_Schematic.pdf", density = 600)
Fig4b <- ggdraw() +draw_image(illust.hb, scale = 1)

AB<- plot_grid(fig4A.plot,Fig4b, labels=c('A','B'), label_size=16, ncol=2, rel_heights = c(1,.8))
AB

fig3row1<- plot_grid(fig4C.plot,fig4D.plot, labels = c('C','D'), label_size = 16, nrow=1, rel_widths=c(1,1), rel_heights= 1, align = "hv", axis="t" )
fig3row1

fig3row2<- plot_grid( fig4E.plot, fig4F.plot, labels = c('E','F'), label_size=16, nrow=1, rel_widths=c(1,1), rel_heights=1,align = "v", axis="t")
fig3row2

fig4 <- plot_grid(AB,fig3row1,fig3row2, nrow = 3, align = "v", axis = 'lr', rel_heights=c(1,.8,.8))
fig4

save_plot("../figures/plots/fig4.pdf", fig4, base_height = 12, base_width = 14)



