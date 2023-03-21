
library(ggplot2)


# ::::::::::::::::::::::::::::::::: figure 2 :::::::::::::::::::::::::::::::::::#

location <- ""                              # give path location where you place the data file
data_fig2 <- read.csv(paste0(location,"data_fig2.csv"), sep="\t") 


source(paste0(location, "function.R"))      # load function; location same as that of the data file

Tars <- plot_function_new_loci_v2(data_fig2, "Tars", ypos = 98.0, ybpos = 98); Tars
Scx <- plot_function_new_loci_v2(data_fig2, "Scx", ypos = 90, ybpos = 85); Scx

Art3 = plot_function_v2(data_fig2, gene_name = "Art3", ypos = 97, ybpos = 97); Art3
Tns1 = plot_function_v2(data_fig2, gene_name = "Tns1", ypos = 97, ybpos = 97)
Smarca2 = plot_function_v2(data_fig2, gene_name = "Smarca2", ypos = 95, ybpos = 92)
Ndrg2 = plot_function_v2(data_fig2, gene_name = "Ndrg2", ypos = 90, ybpos = 87)
Fgfr3 = plot_function_v2(data_fig2, gene_name = "Fgfr3", ypos = 87, ybpos = 84)
Cpn2 = plot_function_v2(data_fig2, gene_name = "Cpn2", ypos = 97, ybpos = 97)
Inpp5a = plot_function_v2(data_fig2, gene_name = "Inpp5a", ypos = 30, ybpos = 27)


figure_all_new<-ggarrange(Tars, Scx, Art3, Tns1, Smarca2, Ndrg2, Fgfr3, Cpn2, Inpp5a,
          labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I") ,
          font.label = list(size = 10),
          ncol = 3 , nrow = 3)
ggsave(plot = figure_all_new, paste0(location, "fig2_v3_all_new.png"), 
       limitsize = FALSE,
       width = 6.81, height = 8.72, units = "in")



## ================= generating for Arhgef19, Tent5a in the Supp.   ======================

Arhgef19 <- plot_function_new_loci_v2(data_fig2, "Arhgef19", ypos = 0.2, ybpos = 0.17)
Tent5a <- plot_function_new_loci_v2(data_fig2, "Tent5a", ypos = 1.00, ybpos = 0.97)


figure <- ggarrange(Arhgef19, Tent5a, 
                    labels = c("A", "B"), 
                    font.label = list(size = 8),
                    ncol = 2, nrow = 1)

ggsave(plot = figure, paste0(location, "fig2_newloci_v3_supl.png"), limitsize = FALSE,
       width = 6.81, height = 3.72, units = "in") #cat(location)

## ========================= generating for the Prdm16 control ===========================
# load data
gene_name = "Prdm16"

data_mod <- data_fig2 %>%
  filter(genes %in% gene_name) %>% #c("Prdm16")) %>%
  mutate(condition = factor(condition, 
                            levels= c("DEEP_Ct", "DEEP_Sh",  "Snupe_Ct", "Snupe_Sh", 
                                      "RRBS_Y",   "RRBS_O",   "miseq_Y",  "miseq_O")))

prdm16 <- ggplot(data_mod, aes(x = condition, y = CpG, colour = factor(clr))) + 
  geom_boxplot(outlier.shape = NA , lwd = 0.4, fatten = 0.4   
               ) + 
               geom_jitter(size = 1.8, alpha = 0.5) +
                 facet_wrap(vars(genes), scales = "free_y", ncol = 3) +
                 # rename the x axis labels
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
figurep <- ggarrange(prdm16, 
                    labels = c("C"), 
                    font.label = list(size = 8),
                    ncol = 1, nrow = 1)
ggsave(plot = figurep, paste0(location, "fig2_newloci_prdm16.png"), limitsize = FALSE,
       width = 6.81, height = 3.72, units = "in") 
