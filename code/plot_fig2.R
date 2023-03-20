#example: https://r-charts.com/distribution/box-plot-jitter-ggplot2/
# install.packages("ggplot2")
library(ggplot2)
# ========================== example plot (using random data) =========================
# Generate Random Data
set.seed(8)
y <- rnorm(200)
group <- sample(LETTERS[1:3], size = 200, replace = TRUE)
df <- data.frame(y, group)          # View(df)

# Box plot by group with jitter
ggplot(df, aes(x = group, y = y,
               colour = group,
               shape = group)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter()

#========================= boxplot with jitter on data =============================

# theme_set(theme_pubr())
setwd("/Volumes/Elements_1/PROJECTS/SONJA-LEBER/") #getwd()

plot_function <- function(gene_name){     #   gene_name = "Art3"
  
  # unique(data_fig2$condition)
  data_fig2$condition <- factor(data_fig2$condition , 
                                levels= c("DEEP_Sh",  "DEEP_Ct", "Snupe_Sh", "Snupe_Ct",  
                                          "RRBS_Y",   "RRBS_O",   "miseq_Y",  "miseq_O"))
  # data_fig2sel = data_fig2[((data_fig2$genes == gene_name) && (data_fig2$panel == "Sh vs Ct")) ,]
  data_fig2sel <- subset(data_fig2, genes==gene_name & panel == "Sh vs Ct")
  data_fig2sel$condition <- factor(data_fig2sel$condition , 
                                levels= c("DEEP_Sh",  "DEEP_Ct", "Snupe_Sh", "Snupe_Ct"))
  s1 <- ggplot(data_fig2sel, aes(x = factor(condition), y = CpG, colour = factor(clr))) + 
    geom_boxplot(outlier.shape = NA, colour = c("gray", "red","gray", "red")) +
    geom_jitter(colour = factor(data_fig2sel$clr)) +
    labs(x = gene_name, y = "CpG Methylation", colour = "Groups")  +
    theme(axis.text.x=element_text(angle=45, size = 10, hjust=1)) + 
    facet_wrap(~ panel, scales = "free_x") + ylim(10, 100) + labs(x = NULL) +
    theme(
      legend.position = "none",
      panel.spacing.x = unit(-8, "lines")) 
      
  s1
  
  if(gene_name == "Lrig1" || gene_name == "Inpp5a" || gene_name == "Smarca2"){
    boxcolour = c("dark green", "blue")
  }else{
    boxcolour = c("dark green", "blue","dark green", "blue")
  }
  data_fig2sel <- subset(data_fig2, genes==gene_name & panel == "Y vs O")
  data_fig2sel$condition <- factor(data_fig2sel$condition , 
                                levels= c("RRBS_Y",   "RRBS_O",   "miseq_Y",  "miseq_O"))
  s2 <- ggplot(data_fig2sel, aes(x = factor(condition), y = CpG, colour = factor(clr))) + 
    geom_boxplot(outlier.shape = NA, colour = boxcolour) +
    geom_jitter(colour = factor(data_fig2sel$clr)) +
    labs(x = gene_name, colour = "Groups")  +
    theme(axis.text.x=element_text(angle=45, size = 10, hjust=1)) + 
    facet_wrap(~ panel) + labs(y = "") + ylim(10, 100) + labs(x = NULL) +
    theme(legend.position = "none", axis.text.y=element_blank(), panel.spacing.x = unit(-8, "lines")) 
  s2

  # library(gridtext)
  
  bottom <- textGrob(gene_name, gp = gpar(fontsize = 10, col="black"))
  p <- grid.arrange(s1, s2, ncol=2, bottom = bottom)
  # p <- plot_grid(s1, s2, labels = "AUTO")
  return(p)
}

data_fig2 <- read.csv("data_fig2.csv", sep="\t") #View(data_fig2) glimpse(data_fig2)

Art3 = plot_function(gene_name = "Art3")
Tns1 = plot_function(gene_name = "Tns1")
Smarca2 = plot_function(gene_name = "Smarca2")
Ndrg2 = plot_function(gene_name = "Ndrg2")
Fgfr3=plot_function(gene_name = "Fgfr3")
Cpn2=plot_function(gene_name = "Cpn2")
Ahdc1=plot_function(gene_name = "Ahdc1")
Lrig1=plot_function(gene_name = "Lrig1")
Inpp5a = plot_function(gene_name = "Inpp5a")
dev.off()
figure <- ggarrange(Art3, Tns1, Smarca2, Ndrg2,Fgfr3,Cpn2,Ahdc1,Lrig1,Inpp5a,
                    labels = c("A", "B", "C", "D","E","F", "G","H", "I"),
                    ncol = 3, nrow = 3)
getwd()
ggsave("fig2_new.pdf", plot = figure, dpi=300)
plot_function(gene_name = "Prdm16") #not run due to missing values



# ::::::::::::::::::::::::::::::::: figure 2 Version 2 :::::::::::::::::::::::::::::::::::#
# load data
location <- "/Volumes/Elements_1/PROJECTS/SONJA-LEBER/"
data_fig2 <- read.csv(paste0(location,"data_fig2.csv"), sep="\t") #View(data_fig2) glimpse(data_fig2)

View(data_fig2)

# load function
source(paste0(location, "function.R"))

Art3 = plot_function_v2(data_fig2, gene_name = "Art3", ypos = 97, ybpos = 97); Art3
Tns1 = plot_function_v2(data_fig2, gene_name = "Tns1", ypos = 97, ybpos = 97)
Smarca2 = plot_function_v2(data_fig2, gene_name = "Smarca2", ypos = 95, ybpos = 92)
Ndrg2 = plot_function_v2(data_fig2, gene_name = "Ndrg2", ypos = 90, ybpos = 87)
Fgfr3 = plot_function_v2(data_fig2, gene_name = "Fgfr3", ypos = 87, ybpos = 84)
Cpn2 = plot_function_v2(data_fig2, gene_name = "Cpn2", ypos = 97, ybpos = 97)
# Ahdc1 = plot_function_v2(data_fig2, gene_name = "Ahdc1", ypos = 20, ybpos = 17)
# Lrig1 = plot_function_v2(data_fig2, gene_name = "Lrig1", ypos = 25, ybpos = 23)
Inpp5a = plot_function_v2(data_fig2, gene_name = "Inpp5a", ypos = 30, ybpos = 27)


figure_cd <- ggarrange(Art3, Tns1, #Smarca2, Ndrg2, Fgfr3, Cpn2, Inpp5a, #Lrig1,Ahdc1,
                    labels = c("C", "D") ,
                    font.label = list(size = 10),
                    ncol = 2 , nrow = 1)
ggsave(plot = figure_cd, paste0(location, "fig2_v3_cd.png"), limitsize = FALSE,
       width = 6.81, height = 3.72, units = "in") #cat(location)

figure_ef <- ggarrange(Smarca2, Ndrg2, #Fgfr3, Cpn2, Inpp5a, #Lrig1,Ahdc1,
                    labels = c("E", "F") ,
                    font.label = list(size = 10),
                    ncol = 2 , nrow = 1)
ggsave(plot = figure_ef, paste0(location, "fig2_v3_ef.png"), limitsize = FALSE,
       width = 6.81, height = 3.72, units = "in") #cat(location)

figure_gh <- ggarrange(Fgfr3, Cpn2, # Inpp5a, #Lrig1,Ahdc1,
                    labels = c("G", "H") ,
                    font.label = list(size = 10),
                    ncol = 2 , nrow = 1)
ggsave(plot = figure_gh, paste0(location, "fig2_v3_gh.png"), limitsize = FALSE,
       width = 6.81, height = 3.72, units = "in") #cat(location)

figure_i <- ggarrange(Inpp5a, #Lrig1,Ahdc1,
                       labels = c("I") ,
                       font.label = list(size = 10),
                       ncol = 1 , nrow = 1)
ggsave(plot = figure_i, paste0(location, "fig2_v3_i.png"), limitsize = FALSE,
       width = 3.4, height = 3.72, units = "in") #cat(location)



figure_all_new<-ggarrange(Tars, Scx, Art3, Tns1, Smarca2, Ndrg2, Fgfr3, Cpn2, Inpp5a,
          labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I") ,
          font.label = list(size = 10),
          ncol = 3 , nrow = 3)
ggsave(plot = figure_all_new, paste0(location, "fig2_v3_all_new.png"), 
       limitsize = FALSE,
       width = 6.81, height = 8.72, units = "in")



# -------------------------------- new loci ---------------------------------
location <- "/Volumes/Elements_1/PROJECTS/SONJA-LEBER/"
data_fig2 <- read.csv(paste0(location,"data_fig2.csv"), sep="\t") #View(data_fig2) glimpse(data_fig2)
source(paste0(location, "function.R"))

# Arhgef19 <- plot_function_new_loci_v2(data_fig2, "Arhgef19", ypos = 0.2, ybpos = 0.17)
# Tent5a <- plot_function_new_loci_v2(data_fig2, "Tent5a", ypos = 1.00, ybpos = 0.97)
Tars <- plot_function_new_loci_v2(data_fig2, "Tars", ypos = 98.0, ybpos = 98); Tars
Scx <- plot_function_new_loci_v2(data_fig2, "Scx", ypos = 90, ybpos = 85); Scx

figure <- ggarrange(Tars, Scx, #Arhgef19, Tent5a, 
                    labels = c("A", "B"), #, "L", "M"),
                    font.label = list(size = 8),
                    ncol = 2, nrow = 1)
figure
ggsave(plot = figure, paste0(location, "fig2_newloci_v3_ab.png"), limitsize = FALSE,
       width = 6.81, height = 3.72, units = "in") #cat(location)

Arhgef19 <- plot_function_new_loci_v2(data_fig2, "Arhgef19", ypos = 0.2, ybpos = 0.17)
Tent5a <- plot_function_new_loci_v2(data_fig2, "Tent5a", ypos = 1.00, ybpos = 0.97)
# Tars <- plot_function_new_loci_v2(data_fig2, "Tars", ypos = 1.0, ybpos = 0.98)
# Scx <- plot_function_new_loci_v2(data_fig2, "Scx", ypos = 0.9, ybpos = 0.85)

figure <- ggarrange(Arhgef19, Tent5a, #Arhgef19, Tent5a, 
                    labels = c("A", "B"), #, "L", "M"),
                    font.label = list(size = 8),
                    ncol = 2, nrow = 1)
figure
ggsave(plot = figure, paste0(location, "fig2_newloci_v3_supl.png"), limitsize = FALSE,
       width = 6.81, height = 3.72, units = "in") #cat(location)

## ================= generating for the Prdm16 control ============
# load data
gene_name = "Prdm16"

data_mod <- data_fig2 %>%
  filter(genes %in% gene_name) %>% #c("Prdm16")) %>%
  mutate(condition = factor(condition, 
                            levels= c("DEEP_Ct", "DEEP_Sh",  "Snupe_Ct", "Snupe_Sh", 
                                      "RRBS_Y",   "RRBS_O",   "miseq_Y",  "miseq_O")))

prdm16 <- ggplot(data_mod, aes(x = condition, y = CpG, colour = factor(clr))) + #, color = factor(clr)
  geom_boxplot(outlier.shape = NA , lwd = 0.4, fatten = 0.4   #color = c("black", "red", "black", "red",
               #"gray", "brown", "gray", "brown"
               # )
               ) + 
               geom_jitter(size = 1.8, alpha = 0.5) +
                 facet_wrap(vars(genes), scales = "free_y", ncol = 3) +
                 #to rename the x axis labels
                 scale_x_discrete(labels=c("DEEP_Ct" = "Co RRBS",
                                           "DEEP_Sh" = "LDC RRBS",  
                                           "Snupe_Ct" = "Co SNuPE",
                                           "Snupe_Sh"= "LDC SNuPE", 
                                           "RRBS_Y" = "Y RRBS",
                                           "RRBS_O" = "O RRBS",
                                           "miseq_Y" = "Y MiSeq",  
                                           "miseq_O" = "O MiSeq"))+
  scale_color_manual(values = c("1" = "red", "2" = "black", "3" = "grey50", "4" = "brown")) +
  ylab("CpG Methylation") +
  xlab("") +
  theme_minimal() +
  theme(
    plot.margin = margin(0,0,0,0, "cm"),
    legend.position = "none",
    axis.text.y = element_text(size = 8),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(angle = 45, size = 8, hjust=1),
    strip.background = element_blank(), 
    strip.placement = "outside", 
    strip.text = element_text(face="bold", size=9))
figurep <- ggarrange(prdm16, #Arhgef19, Tent5a, 
                    labels = c("C"), #, "L", "M"),
                    font.label = list(size = 8),
                    ncol = 1, nrow = 1)
ggsave(plot = figurep, paste0(location, "fig2_newloci_prdm16.png"), limitsize = FALSE,
       width = 6.81, height = 3.72, units = "in") #cat(location)
