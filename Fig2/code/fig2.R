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
library(magrittr)
library(broom)

library(purrr)
library(here)

#PDF import
library(magick)
library(pdftools)
library(grImport2)

conflict_prefer("filter", "dplyr")
conflict_prefer("summarize", "dplyr")

#------------------------------------set working directory-----------------------------
setwd<- here()

#-------------------------------------read in data--------------------------------------
scope_data<- read_csv("../data/rnascope/RNAscope.csv")

#set theme 
aesthetic <- 
  theme_zlab() +
  theme(
    panel.grid.major.y = element_line(),
    axis.text.x = element_markdown(angle = 0, hjust = .5, face="plain",size=10),
    axis.text.y = element_text(face = "plain", size = 10),
    strip.text.x = element_markdown(size=10),
    strip.text.y = element_markdown(size=12),
    strip.text = ggplot2::element_text(face = "plain", size = 10),
    axis.title.y = element_text(size=12, face= "plain"),    
    axis.title.x = element_text(angle = 0,size=12, face= "plain"))
#-------------------------------------RNAscope-------------------------------------------

#data long form 
scope_data<- scope_data%>%
  mutate(Bad = ifelse(is.na(Bad), F, Bad))%>%
  filter(Bad == FALSE)%>%
  pivot_longer(cols= c(3:12,14:19), names_to = "Tissue", values_to = "Count")%>%
  select(1,2,4,5,6)


#rename
scope_data$Tissue <- factor(scope_data$Tissue , 
  levels = c("Total_KJG","Total_PMA","Body_wall","Cuticle","Excretory_nervous","Other","Pharynx","Intestine",
              "Reproductive","total_lateralcord", "Pharynx_PMA","Reproductive_PMA","Body_wall_PMA","Cuticle_PMA",
              "Intestine_PMA","total_lateralcord_PMA"), 
  labels=c('Total_kjg','Total_pma','Body wall','Dermal layers','es','Other','Pharynx','Intestine','Reproductive',
            'Lateral cords','Pharynx_pma','Reproductive_pma','Body wall_pma','Dermal layers_pma',
             'Intestine_pma','Lateral cords_pma'))
#tidy
scope_data<- scope_data%>%
  filter(!Tissue %in% c( "es", "Other"))

scope_data<- scope_data%>%
separate( col = Tissue, into = c("Tissue", "person"), sep = "_")

scope_data<- scope_data%>%
  mutate(person = ifelse(is.na(person), 'kjg', person))%>%
mutate(Count = ifelse(is.na(Count), 0, Count))

scope_mean<- scope_data%>%
  group_by(Distance, Tissue) %>% 
  summarize(Count = mean(Count))%>%
  mutate(person= NA)

#plot----------------------------------------------------------------------
scope_plot <- ggplot(scope_data, aes(x = Distance, y = Count, color= person, fill=person)) +
  geom_line(data= scope_data, orientation = "x" , size=1)+
  geom_point(size = .2)+
  facet_grid(rows = vars(Tissue), scales= "free_x") +
  theme_minimal(base_size = 16, base_family = "Helvetica") +
 # scale_y_reverse()+
  geom_line(data= scope_mean, orientation = "x" , size=.5, color= "black", alpha=.8)+
  aesthetic+
  ylab("Puncta") +
  scale_color_manual(values = c("#84d6cc","#FFCFA4","#84D6CC")) +
  scale_fill_manual(values = c("#84d6cc","#FFCFA4","#84D6CC")) +
  xlab(expression(paste("Distance from anterior (µm)")))+
  theme(panel.border = element_rect(color="black", fill=NA),
        axis.title.y = element_text(size=12, face= "plain"),
        panel.grid.major.x = element_line(color = "grey",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.minor.x = element_line(color = "grey",
                                          size = 0.5,
                                          linetype = 1))
scope_plot


#----------------------------percent staining 
#refresh data
scope_data<- read_csv("../data/rnascope/RNAscope.csv")

#calculate percent 
df<-as.data.frame(t(colSums(!is.na(scope_data))))%>%
  pivot_longer(cols= c(1:19), names_to = "Tissue", values_to = "Count")%>%
  mutate( 
    percent = case_when(
      Tissue == "Pharynx" ~ Count/123, 
      Tissue == "Pharynx_PMA" ~ Count/123,
      Tissue == "total_lateralcord" ~ Count/202,
      Tissue == "total_lateralcord_PMA" ~ Count/202,
      Tissue == "Total_KJG" ~ Count/214,
      Tissue == "Total_PMA" ~ Count/214,
      Tissue == "Body_wall" ~ Count/214,
      Tissue == "Body_wall_PMA" ~ Count/214,
      Tissue == "Reproductive" ~ Count/130,
      Tissue == "Reproductive_PMA" ~ Count/130,
      Tissue == "Intestine" ~ Count/91,
      Tissue == "Intestine_PMA" ~ Count/91,
      Tissue == "Cuticle" ~ Count/214,
      Tissue == "Cuticle_PMA" ~ Count/214,
    ))

#rename 
df$Tissue <- factor(df$Tissue , 
                            levels = c("Distance","Slice","lateral_cord","Total_KJG","Total_PMA","Body_wall","Cuticle","Excretory_nervous","Other","Pharynx","Intestine",
                                       "Reproductive","total_lateralcord", "Pharynx_PMA","Reproductive_PMA","Body_wall_PMA","Cuticle_PMA",
                                       "Intestine_PMA","total_lateralcord_PMA"), 
                            labels=c("Distance","Slice","lateral_cord",'Total_kjg','Total_pma','Body wall','Dermal layers','es','Other','Pharynx','Intestine','Reproductive',
                                     'Lateral cords','Pharynx_pma','Reproductive_pma','Body wall_pma','Dermal layers_pma',
                                     'Intestine_pma','Lateral cords_pma'))
df<- df%>%
  filter(!Tissue %in% c( "es", "Other","Distance","Slice","lateral_cord"))

df<- df%>%
  separate( col = Tissue, into = c("Tissue", "person"), sep = "_")

df<- df%>%
  mutate(person = ifelse(is.na(person), 'kjg', person))%>%
  mutate(Count = ifelse(is.na(Count), 0, Count))
  

#refresh theme 
aesthetic <- 
  theme_zlab() +
  theme(
    panel.grid.major.y = element_line(),
    axis.text.x = element_markdown(angle = 0, hjust = .5, face="plain",size=10),
    axis.text.y = element_text(face = "plain", size = 10),
    strip.text.x = element_markdown(size=10),
    strip.text.y = element_markdown(size=12),
    strip.text = ggplot2::element_text(face = "plain", size = 10),
    axis.title.y = element_text(size=12, face= "plain"),    
    axis.title.x = element_text(angle = 0,size=12, face= "plain"))

#plot
box_plot <- ggplot(df, aes(y =Tissue , x = percent, color= person, fill= person)) +
  geom_point(size=3)+
  theme_minimal(base_size = 16, base_family = "Helvetica") +
  stat_summary(
    geom = "point",
    fun.y = "mean",
    col = "black",
    size = 2,fill = "black" )+
aesthetic+
  xlab("Sections with Staining") +
  scale_color_manual(values = c("#84d6cc","#FFCFA4","#84D6CC")) +
  scale_fill_manual(values = c("#84d6cc","#FFCFA4","#84D6CC")) +
  scale_x_continuous(labels = scales::percent, limits = c(0, 1.1), breaks = seq(0, 1.50, .2),expand = c(0, 0))+
  scale_color_manual(values = c("#84d6cc","#FFCFA4","#84D6CC")) +
  scale_fill_manual(values = c("#84d6cc","#FFCFA4","#84D6CC")) +
  theme(panel.border = element_rect(color="black", fill=NA),
        axis.title.y = element_text(size=12, face= "plain"),
        panel.grid.major.x = element_line(color = "grey",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.minor.x = element_line(color = "grey",
                                          size = 0.5,
                                          linetype = 1))
box_plot


#------------------------------------------FIGURE 2------------------------------------- 

bm_adult <- image_read_pdf("../images/bm_adult.pdf", density = 600)
bm_sections <- image_read_pdf("../images/bm_sections_inset_red.pdf", density = 550)


Fig2_top <- ggdraw() + draw_image(bm_adult, scale = 1)
bm_sections <- ggdraw() + draw_image(bm_sections, scale = 1)

Fig2_bottom_left <- scope_plot
Fig2_bottom_right <- plot_grid(bm_sections, box_plot, labels = c('C', 'D'), rel_heights = c(2.5,0.7), ncol = 1)

Fig2_bottom <- plot_grid(Fig2_bottom_left, Fig2_bottom_right,labels = c('B', ''), ncol =2)

Fig2 <- plot_grid(Fig2_top,Fig2_bottom,labels = c('A', ''), rel_heights = c(0.5,2), ncol =1, scale = 0.99)

save_plot("../figures/plots/fig2.pdf", Fig2, base_height = 11.5, base_width = 13)




