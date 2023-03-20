library(ggplot2)
library(ggpubr)
library(cowplot)
library(gridExtra)
library(grid)
library(dplyr)

plot_function_v2 <- function(
  data ,
  gene_name, ypos, ybpos
  ){ 
  
  
  cat("Within the plot function", "\n") 
  
  
  data_mod <- data %>% #mutate(CpG = CpG * 100) %>%
    filter(genes %in% gene_name) %>% #c("Prdm16")) %>%
    mutate(condition = factor(condition, 
                              levels= c("DEEP_Ct", "DEEP_Sh",  "Snupe_Ct", "Snupe_Sh", 
                              "RRBS_Y",   "RRBS_O",   "miseq_Y",  "miseq_O")))
  
  # checking the wilcox test results
  #View(compare_means(CpG ~ condition, data = data_mod))
  
  Deep <- compare_means(CpG ~ condition, data = data_mod, method="t.test", p.adjust.method="BH"
) %>%
    filter(group1=="DEEP_Ct"	& group2=="DEEP_Sh")
  snupe <- compare_means(CpG ~ condition, data = data_mod, method="t.test" , p.adjust.method="BH") %>%
    filter(group1=="Snupe_Ct"	& group2=="Snupe_Sh")
  if(! gene_name %in% c("Smarca2","Inpp5a", "Lrig1")){
    rrbs <- compare_means(CpG ~ condition, data = data_mod, method="t.test", p.adjust.method="BH") %>%
    filter(group1=="RRBS_Y"	& group2=="RRBS_O")
    if(rrbs$p.signif == "ns"){rrbs$p.signif = ""}
    if(rrbs$p.signif == "**"){rrbs$p.signif = "*"}
  }
  miseq <- compare_means(CpG ~ condition, data = data_mod, method="t.test", p.adjust.method="BH") %>%
    filter(group1=="miseq_Y"	& group2=="miseq_O")
  
  print(Deep$p.signif)
  print(snupe$p.signif)
  print(miseq$p.signif)
  
  if(miseq$p.signif == "ns"){miseq$p.signif = ""}
  if(Deep$p.signif == "ns"){Deep$p.signif = ""}
  if(snupe$p.signif == "ns"){snupe$p.signif = ""}
  if(miseq$p.signif == "**"){miseq$p.signif = "*"}
  if(Deep$p.signif == "**"){Deep$p.signif = "*"}
  if(snupe$p.signif == "**"){snupe$p.signif = "*"}
  
  
  p <- ggplot(data_mod, aes(x = condition, y = CpG, colour = factor(clr))) + #, color = factor(clr)
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
                              "RRBS_O" = "A RRBS",
                              "miseq_Y" = "Y MiSeq",  
                              "miseq_O" = "A MiSeq")) +
    scale_color_manual(values = c("1" = "red", "2" = "black", "3" = "grey50", "4" = "brown")) +
    ylab("CpG Methylation") +
    xlab("") +
    theme_minimal() +
    theme(
      plot.margin = margin(0,0,0,0, "cm"),
      legend.position = "none",
      axis.text.y = element_text(size = 8),
      axis.title.y = element_text(size = 10),
      axis.text.x = element_text(angle = 45, size = 10, hjust=1),
      strip.background = element_blank(), 
      strip.placement = "outside", 
      strip.text = element_text(face="bold", size=9)) + #+ stat_compare_means()
    geom_text(x=1.5, y=ypos, label=Deep$p.signif, size = 8, colour = "black") +
    geom_text(x=3.5, y=ypos, label=snupe$p.signif, size = 8, colour = "black")
    
  
    if(! gene_name %in% c("Smarca2","Inpp5a", "Lrig1")){
      p <- p + geom_text(x=5.5, y=ypos, label=rrbs$p.signif, size = 8, colour = "black") +
        geom_text(x=7.5, y=ypos, label=miseq$p.signif, size = 8, colour = "black") #+
      #   geom_line(data = tibble(x = c(1,2), y = c(ybpos,ybpos)), aes(x=x, y=y), colour = "black" ) +
      #   geom_line(data = tibble(x = c(3,4), y = c(ybpos,ybpos)), aes(x=x, y=y), colour = "black" ) +
      #   geom_line(data = tibble(x = c(5,6), y = c(ybpos,ybpos)), aes(x=x, y=y), colour = "black" ) +
      #   geom_line(data = tibble(x = c(7,8), y = c(ybpos,ybpos)), aes(x=x, y=y), colour = "black" )
    }else{
      p <- p + geom_text(x=5.5, y=ypos, label=miseq$p.signif, size = 8, colour = "black") #+
      #   geom_line(data = tibble(x = c(1,2), y = c(ybpos,ybpos)), aes(x=x, y=y), colour = "black" ) +
      #   geom_line(data = tibble(x = c(3,4), y = c(ybpos,ybpos)), aes(x=x, y=y), colour = "black" ) +
      #   geom_line(data = tibble(x = c(5,6), y = c(ybpos,ybpos)), aes(x=x, y=y), colour = "black" )
    }
  
  p
  # cat("Figure saved in ", location, "\n")
  # ggsave(paste0(location, "fig2_v2.pdf"), dpi=300) #cat(location)
  
}

plot_function_new_loci_v2 <- function(
  data ,
  gene_name,
  ypos, ybpos
){ 
  
  
  cat("Within the plot function of novel loci", "\n") 
  
  data_mod <- data %>% mutate(CpG = CpG * 100) %>%
    filter(genes %in% gene_name) %>% #c("Prdm16")) %>%
    mutate(condition = factor(condition, 
                              levels= c("miseq_DEEP_Ct", "miseq_DEEP_Sh",  "Snupe_DEEP_Ct", "Snupe_DEEP_Sh", 
                                        "Snupe_Y",   "Snupe_O",   "miseq_Y",  "miseq_O")))
  
  
  # generate labels for significance
  # View(compare_means(CpG ~ condition, data = data_mod))
  
  miseq_DEEP <- compare_means(CpG ~ condition, data = data_mod, method="t.test",p.adjust.method = "BH") %>%
    filter(group1=="miseq_DEEP_Ct"	& group2=="miseq_DEEP_Sh")
  
  Snupe_DEEP <- compare_means(CpG ~ condition, data = data_mod, method="t.test",p.adjust.method = "BH") %>%
     filter(group1=="Snupe_DEEP_Ct"	& group2=="Snupe_DEEP_Sh")
  
  Snupe <- compare_means(CpG ~ condition, data = data_mod, method="t.test",p.adjust.method = "BH") %>%
       filter(group1=="Snupe_Y"	& group2=="Snupe_O")
  
  miseq <- compare_means(CpG ~ condition, data = data_mod, method="t.test",p.adjust.method = "BH") %>%
     filter(group1=="miseq_Y"	& group2=="miseq_O")
  
  print(miseq_DEEP$p.signif)
  print(Snupe_DEEP$p.signif)
  print(Snupe$p.signif)
  print(miseq$p.signif)

  if(miseq_DEEP$p.signif == "ns"){miseq_DEEP$p.signif = ""}
  if(Snupe_DEEP$p.signif == "ns"){Snupe_DEEP$p.signif = ""}
  if(Snupe$p.signif == "ns"){Snupe$p.signif = ""}
  if(miseq$p.signif == "ns"){miseq$p.signif = ""}
  if(Snupe$p.signif == "***"){Snupe$p.signif = "*"}
  if(miseq$p.signif == "***"){miseq$p.signif = "*"}

  ggplot(data_mod, aes(x = condition, y = CpG, colour = factor(clr))) + #, color = factor(clr)
    geom_boxplot(outlier.shape = NA , lwd = 0.4, fatten = 0.4   #color = c("black", "red", "black", "red",
                 #"gray", "brown", "gray", "brown"
                 # )
    ) + 
    geom_jitter(size = 1.8, alpha = 0.5) +
    facet_wrap(vars(genes), scales = "free_y", ncol = 3) +
    #to rename the x axis labels
    scale_x_discrete(labels=c("miseq_DEEP_Ct" = "Co MiSeq",      #\n(DEEP)",
                              "miseq_DEEP_Sh" = "LDC MiSeq",     #\n(DEEP)",  
                              "Snupe_DEEP_Ct" = "Co SNuPE",      #\n(DEEP)",
                              "Snupe_DEEP_Sh"= "LDC SNuPE",      #\n(DEEP)", 
                              "Snupe_Y" = "Y SNuPE",
                              "Snupe_O" = "A SNuPE",
                              "miseq_Y" = "Y MiSeq",  
                              "miseq_O" = "A MiSeq")) +
    scale_color_manual(values = c("1" = "red", "2" = "black", "3" = "grey50", "4" = "brown")) +
    ylab("CpG Methylation") +
    xlab("") +
    theme_minimal() +
    theme(
      plot.margin = margin(0,0,0,0, "cm"),
      legend.position = "none",
      axis.text.y = element_text(size = 8),
      axis.title.y = element_text(size = 10),
      axis.text.x = element_text(angle = 45, size = 10, hjust=1),
      strip.background = element_blank(), 
      strip.placement = "outside", 
      strip.text = element_text(face="bold", size=9)) + #+ stat_compare_means()
    geom_text(x=1.5, y=ypos, label=miseq_DEEP$p.signif, size = 8, colour = "black" ) + #miseq_DEEP$p.signif
    geom_text(x=3.5, y=ypos, label=Snupe_DEEP$p.signif, size = 8, colour = "black" ) + #Snupe_DEEP$p.signif
    geom_text(x=5.5, y=ypos, label=Snupe$p.signif, size = 8, colour = "black" ) +  #Snupe$p.signif
    geom_text(x=7.5, y=ypos, label=miseq$p.signif, size = 8, colour = "black" ) #+ #miseq$p.signif
#     geom_line(data = tibble(x = c(1,2), y = c(ybpos,ybpos)), aes(x=x, y=y), colour = "black" ) +
#     geom_line(data = tibble(x = c(3,4), y = c(ybpos,ybpos)), aes(x=x, y=y), colour = "black" ) +
#     geom_line(data = tibble(x = c(5,6), y = c(ybpos,ybpos)), aes(x=x, y=y), colour = "black" ) +
#     geom_line(data = tibble(x = c(7,8), y = c(ybpos,ybpos)), aes(x=x, y=y), colour = "black" )
  
  # cat("Figure saved in ", location, "\n")
  # ggsave(paste0(location, "fig2_v2.pdf"), dpi=300) #cat(location)
  
}