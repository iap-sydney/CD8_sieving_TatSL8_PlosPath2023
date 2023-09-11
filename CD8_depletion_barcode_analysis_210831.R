#Steffen Docken
#31-08-21
#This code will examine the barcodes detected before and after CD8 depletion

rm(list=ls())
graphics.off()


library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(readxl)
library(cowplot)
library(latex2exp)

span_val = 0.5

#color palettes:
plotPalette <- c(brewer.pal(n = 12,name = "Paired"), "#000000", "#999999", "#FFFFFF") 
#black, grey, and white are last two colors

names(plotPalette) <- factor(c("unused", "unused", "unused", "unused",
                               "Pre", "ATI-1",
                               "Allen et al. 2000.Primary", 
                               "Immonen et al. 2020.Primary",
                               "unused", "Post","unused", "unused",
                               "FTY.Pre-ART", "unused", "Never"))
  
Barcode_timingScale_plotPalette_fill <- 
  scale_colour_manual(name = "Barcode Timing",values = plotPalette,
                      breaks = c("Pre", "Post"),
                      labels = c("Pre-CD8 depletion", "Post-CD8 depletion"),
                      aesthetics = "fill", limits = force)

Barcode_timingScale_plotPalette_color <- 
  scale_colour_manual(name = "Barcode Timing",values = plotPalette,
                      breaks = c("Pre", "Post"),
                      labels = c("Pre-CD8 depletion","Post-CD8 depletion"),
                      aesthetics =  "color", limits = force)

Barcode_timingNdetScale_plotPalette_fill <- 
  scale_colour_manual(name = "Barcode Timing",values = plotPalette,
                      breaks = c("Pre", "Post", "Never"),
                      aesthetics =  "fill", limits = force)



ART_annotation = plotPalette[14]

save_paper_fig_flag = 1; #indicate if save figures

paper_symbol_size = 1.5
paper_text_size = 10

paper_fig_location = "Figures/"

symbol_size = paper_symbol_size

perc_linear_ticks = seq(0, 1, .25)
perc_linear_labels = c(TeX(r"($0\%$)"), TeX(r"($25\%$)"), TeX(r"($50\%$)"),
                          TeX(r"($75\%$)"), TeX(r"($100\%$)"))

Barcode_data_from_MATLAB <- read.csv("CD8_analysis/PrePost_CD8depl_barcodes_data_for_R_230120.csv")


CD8depl_bcode_data_df <- data.frame(Animal = Barcode_data_from_MATLAB$Animal,
                                    ATI2_timing = Barcode_data_from_MATLAB$ATI2.timing,
                                    Barcode_ID = Barcode_data_from_MATLAB$Barcode.ID,
                                    Pre_Post = Barcode_data_from_MATLAB$Pre.vs.Post,
                                    Last_det = Barcode_data_from_MATLAB$Last.Detection,
                                    Frac_WT_last_det = Barcode_data_from_MATLAB$WT.frac.at.Last.Det)

## Viral Load Data
VL_data_from_MATLAB <- read.csv("VL_timecourse_data_for_R.csv")

VL_df = data.frame(VL = VL_data_from_MATLAB$Viral.Load..copies.ml.,
                   VL_det = VL_data_from_MATLAB$Viral.Load.detection,
                   Animal = VL_data_from_MATLAB$Animal,
                   Collection_Group = VL_data_from_MATLAB$Collection.Group,
                   FTY_treat = VL_data_from_MATLAB$FTY.Treatment,
                   dpi = VL_data_from_MATLAB$dpi,
                   CD8depl_dpi = rep(0, length(VL_data_from_MATLAB$dpi)))

Animal_treat_df = unique(VL_df[c("Animal","FTY_treat")])

t_rescale_pars= list(PreARTend = 45, gap = 40, ATI1start = 210, ATI1end = 290,
                     ATI2start = 365, CD8_contract = 1/3, CD8depl = 470)
ART_disp = data.frame(ART_start = c(t_rescale_pars$PreARTend, 
                                    (t_rescale_pars$ATI1end- t_rescale_pars$ATI1start)+
                                      t_rescale_pars$PreARTend + t_rescale_pars$gap))
ART_disp = mutate(ART_disp, ART_end = ART_start + t_rescale_pars$gap)

CD8depl_vis_df = data.frame(t = rep(t_rescale_pars$CD8depl, 2),
                            log10VL = c(0, 2))

last_det_timepoints = c('Never', 'Primary', 'Early ATI-1', 'Late ATI-1')
last_det_t_plot = c(-t_rescale_pars$PreARTend/2, t_rescale_pars$PreARTend/2, 
                    t_rescale_pars$PreARTend + t_rescale_pars$gap + 
                      (t_rescale_pars$ATI1end-t_rescale_pars$ATI1start)/4,
                    t_rescale_pars$PreARTend + t_rescale_pars$gap + 
                      3*(t_rescale_pars$ATI1end-t_rescale_pars$ATI1start)/4)

num_last_det_timepoints = length(last_det_timepoints)

Last_det_time_count_df = data.frame(bcode_count = rep(0, 2*num_last_det_timepoints))

CD8depl_bcode_data_df[CD8depl_bcode_data_df$Last_det == "Never", ]$Frac_WT_last_det =
  1#setting fraction WT at last detection to be 1, since last detection was 
#the inoculum

CD8depl_bcode_data_df = merge(x = CD8depl_bcode_data_df, Animal_treat_df,
                              by = "Animal", all.x = TRUE)
 
for (kk in 1:4) {
  
  print(switch(kk, "All Animals", "Not RNd15 or RYm15", "FTY", "Control"))

  CD8depl_bcode_data_df_it = switch(kk, CD8depl_bcode_data_df, 
                                    filter(CD8depl_bcode_data_df, Animal != 'RNd15' &Animal != 'RYm15'),
                                    filter(CD8depl_bcode_data_df, FTY_treat == 'FTY'),
                                    filter(CD8depl_bcode_data_df, FTY_treat == 'Control'))
  print(length(CD8depl_bcode_data_df_it$Frac_WT_last_det))
  
  for (ii in 1:2) {
    Pre_Post = switch(ii, 'Pre', 'Post')
    for (jj in 1:num_last_det_timepoints) {
      last_det_time = last_det_timepoints[jj]
      
      Last_det_time_count_df$Pre_Post[(ii-1)*num_last_det_timepoints + jj] = Pre_Post
      Last_det_time_count_df$Last_det[(ii-1)*num_last_det_timepoints + jj] = last_det_time
      
      Last_det_time_count_df$bcode_count[(ii-1)*num_last_det_timepoints + jj] =
        length(which((CD8depl_bcode_data_df_it$Pre_Post == Pre_Post)&(CD8depl_bcode_data_df_it$Last_det == last_det_time)))
    }
  }
  
  Last_det_time_count_df = Last_det_time_count_df %>% group_by(Pre_Post) %>%
    mutate(Last_det_plot_t = case_when(Last_det == last_det_timepoints[1] ~ last_det_t_plot[1],
                                       Last_det == last_det_timepoints[2] ~ last_det_t_plot[2],
                                       Last_det == last_det_timepoints[3] ~ last_det_t_plot[3],
                                       Last_det == last_det_timepoints[4] ~ last_det_t_plot[4]),
           bcode_frac = bcode_count/sum(bcode_count),
           labels = scales::percent(bcode_frac))
  
  Last_det_time_bcode_panel = ggplot(Last_det_time_count_df,
                                 aes(x = Last_det_plot_t, 
                                     y = bcode_count, fill = Pre_Post)) +
                                 geom_bar(stat = "identity", 
                                          position = position_dodge2(reverse = TRUE)) + 
    Barcode_timingScale_plotPalette_fill +
    scale_x_continuous(breaks = last_det_t_plot, 
                       labels = c("Never", "Primary", "Early ATI 1", "Late ATI 1"))  +
    scale_y_continuous(limits = c(0, 65)) +
    geom_rect(data = ART_disp, aes(xmin = ART_start, xmax = ART_end, 
                                ymin=0, ymax=65), 
              fill = ART_annotation, alpha = .2, inherit.aes = FALSE)+
    geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
    annotate("text", x = c(t_rescale_pars$PreARTend+.5*t_rescale_pars$gap,
                           t_rescale_pars$PreARTend+1.5*t_rescale_pars$gap + 
                             (t_rescale_pars$ATI1end - t_rescale_pars$ATI1start)), 
             y = 20, label = "ART", angle = 90) +
    theme_minimal()+
    labs(x = "Last Detection", y = c(expression(atop("Number of", "Barcodes"))))+
    theme(title = element_text(size = paper_text_size), 
          axis.text = element_text(size = paper_text_size),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, 
                                     margin = margin(t = 0, r = 0, b = 0, l = 0)),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())
  
  print(Last_det_time_bcode_panel)
  
  #pie charts of fraction detected each time
  pie_chart_Pre_init = ggplot(subset(Last_det_time_count_df, Pre_Post == "Pre"), 
                          aes(x = "", y = bcode_frac, 
                              fill = factor(Last_det,
                                            levels = c("Never", "Primary", 
                                                       "Early ATI-1", 
                                                       "Late ATI-1")))) + 
    geom_bar(width = 1, stat = "identity", color = "black") + coord_polar(theta = "y", start = 0,
                                                         direction = -1)+
    theme_void()+theme(legend.position = "none")# + geom_text(aes(label = labels),
                              #position = position_stack(vjust = 0.5))
  
  pie_chart_Post_init = ggplot(subset(Last_det_time_count_df, Pre_Post == "Post"), 
                               aes(x = "", y = bcode_frac, 
                                   fill = factor(Last_det,
                                                 levels = c("Never", "Primary", 
                                                            "Early ATI-1", 
                                                            "Late ATI-1")))) + 
    geom_bar(width = 1, stat = "identity", color = "black") + coord_polar(theta = "y", start = 0,
                                                         direction = -1)+
    theme_void()+ theme(legend.position = "none")# + geom_text(aes(label = labels),
                               #position = position_stack(vjust = 0.5)) 
  
  Pre_color = plotPalette
  Post_color = plotPalette
  names(Pre_color)[5] = "Never"
  names(Post_color)[10] = "Never" #will change the name of the color needed 
  #for each version of the plot
  
  pie_chart_Pre_Never = pie_chart_Pre_init + 
    scale_fill_manual(values = Pre_color, breaks = "Never", na.value = "white")
  
  pie_chart_Post_Never = pie_chart_Post_init + 
    scale_fill_manual(values = Post_color, breaks = "Never", na.value = "white")
  
  names(Pre_color)[5] = "Primary"
  pie_chart_Pre_Primary = pie_chart_Pre_init + 
    scale_fill_manual(values = Pre_color, breaks = "Primary", na.value = "white")
  
  names(Post_color)[10] = "Primary"
  pie_chart_Post_Primary = pie_chart_Post_init + 
    scale_fill_manual(values = Post_color, breaks = "Primary", na.value = "white")
  
  names(Pre_color)[5] = "Early ATI-1"
  pie_chart_Pre_Early_ATI1 = pie_chart_Pre_init + 
    scale_fill_manual(values = Pre_color, breaks = "Early ATI-1", na.value = "white")
  
  names(Post_color)[10] = "Early ATI-1"
  pie_chart_Post_Early_ATI1 = pie_chart_Post_init + 
    scale_fill_manual(values = Post_color, breaks = "Early ATI-1", na.value = "white")
  
  names(Pre_color)[5] = "Late ATI-1"
  pie_chart_Pre_Late_ATI1 = pie_chart_Pre_init + 
    scale_fill_manual(values = Pre_color, breaks = "Late ATI-1", na.value = "white")

  names(Post_color)[10] = "Late ATI-1"
  pie_chart_Post_Late_ATI1 = pie_chart_Post_init + 
    scale_fill_manual(values = Post_color, breaks = "Late ATI-1", na.value = "white")
  
  #Full Time of last detection Figure
  
  legend_last_det <- get_legend(Last_det_time_bcode_panel + 
                                  theme(legend.box.margin = margin(0, 0, 0, 0),
                                        title = element_text(size = paper_text_size), 
                                        text = element_text(size = paper_text_size)))
  
  if (save_paper_fig_flag == 1) {
    
    #ggsave(filename = switch(kk, '../../../HIV_Research/paper_manuscripts/SL8_Evolution_FTY_Study_2020/Figures/CD8_figure/Last_det_Old_New_barcodes_all_animals.eps',
    #              '../../../HIV_Research/paper_manuscripts/SL8_Evolution_FTY_Study_2020/Figures/CD8_figure/Last_det_Old_New_barcodes_notRNd15_or_RYm15.eps'), 
    #       plot = Last_det_time_plot + theme(legend.position = 'none'), width = 8, height = 5, units = "cm", device = "eps")#previous version of figure
    
    ggsave(filename = switch(kk, paste0(paper_fig_location, 'Last_det_pre_postCD8depl_barcodes_all_animals_220816.pdf'),
                             paste0(paper_fig_location, 'Extra_CD8_depletion_figs/Last_det_pre_postCD8depl_barcodes_notRNd15_or_RYm15_220816.pdf'),
                             paste0(paper_fig_location, 'Extra_CD8_depletion_figs/Last_det_pre_postCD8depl_barcodes_FTY_animals_220816.pdf'),
                             paste0(paper_fig_location, 'Extra_CD8_depletion_figs/Last_det_pre_postCD8depl_barcodes_Control_animals_220816.pdf')), 
           plot = Last_det_time_bcode_panel + theme(legend.position = 'none'), width = 7, height = 8, units = "cm", device = "pdf")
    
    if(kk == 1){
      ggsave(filename = paste0(paper_fig_location, 'pie_chart_Pre_Never_all_animals_221116.pdf'), 
             plot = pie_chart_Pre_Never, width = 8, height = 8, units = "cm", device = "pdf") 
      ggsave(filename = paste0(paper_fig_location, 'pie_chart_Pre_Primary_all_animals_221116.pdf'), 
             plot = pie_chart_Pre_Primary, width = 8, height = 8, units = "cm", device = "pdf") 
      ggsave(filename = paste0(paper_fig_location, 'pie_chart_Pre_Early_ATI1_all_animals_221116.pdf'), 
             plot = pie_chart_Pre_Early_ATI1, width = 8, height = 8, units = "cm", device = "pdf") 
      ggsave(filename = paste0(paper_fig_location, 'pie_chart_Pre_Late_ATI1_all_animals_221116.pdf'), 
             plot = pie_chart_Pre_Late_ATI1, width = 8, height = 8, units = "cm", device = "pdf") 
      
      ggsave(filename = paste0(paper_fig_location, 'pie_chart_Post_Never_all_animals_221116.pdf'), 
             plot = pie_chart_Post_Never, width = 8, height = 8, units = "cm", device = "pdf")
      ggsave(filename = paste0(paper_fig_location, 'pie_chart_Post_Primary_all_animals_221116.pdf'), 
             plot = pie_chart_Post_Primary, width = 8, height = 8, units = "cm", device = "pdf") 
      ggsave(filename = paste0(paper_fig_location, 'pie_chart_Post_Early_ATI1_all_animals_221116.pdf'), 
             plot = pie_chart_Post_Early_ATI1, width = 8, height = 8, units = "cm", device = "pdf") 
      ggsave(filename = paste0(paper_fig_location, 'pie_chart_Post_Late_ATI1_all_animals_221116.pdf'), 
             plot = pie_chart_Post_Late_ATI1, width = 8, height = 8, units = "cm", device = "pdf") 
      
    }
    
  }
  
  M <- as.table(cbind(Last_det_time_count_df$bcode_count[Last_det_time_count_df$Pre_Post == "Pre"],
                      Last_det_time_count_df$bcode_count[Last_det_time_count_df$Pre_Post == "Post"]))
  dimnames(M) <- list(Last_Det = c("Never", "Primary", "Early ATI-1", "Late ATI-1"),
                      Pre_Post = c("Pre", "Post"))
  
  Last_det_test = chisq.test(M)
  
  
  print(M)
  print(Last_det_test$p.value)
  
  # Fraction WT Last detection
  CD8depl_bcode_data_df_it = CD8depl_bcode_data_df_it %>% 
    mutate(Pre_Post_Never = if_else(Last_det == "Never", Last_det, Pre_Post))
  #creating a dummy varialbe that contains "Never" for if the barcode was never
  #previously detected and when detected during ATI-2 if previously detected
  
  WT_last_det_plot = ggplot(data = CD8depl_bcode_data_df_it,
                            aes(x = Pre_Post,y = Frac_WT_last_det, 
                                color = Pre_Post, fill = Pre_Post_Never)) +
    geom_dotplot(binaxis = 'y', stackdir = 'center', binwidth = 0.02)+ 
    Barcode_timingScale_plotPalette_color + 
    Barcode_timingNdetScale_plotPalette_fill +
    scale_y_continuous(breaks = perc_linear_ticks, labels = perc_linear_labels) +
    scale_x_discrete(limits = c("Pre", "Post"), breaks = c("Pre", "Post"), 
                     labels = c(expression(atop("Pre-CD8", "Depletion")), 
                                expression(atop("Post-CD8", "Depletion")))) + 
    theme_classic() + labs(y = "Tat-SL8 Wild-Type (%)") + 
    theme(legend.position = 'none', axis.title.x = element_blank(), 
          axis.text = element_text(size = paper_text_size),
          axis.title = element_text(size = paper_text_size))
  
  
  print(WT_last_det_plot)
  
  if (save_paper_fig_flag == 1) {
    ggsave(filename = switch(kk, paste0(paper_fig_location, 'Frac_WT_Pre_Post_barcodes_all_animals.eps'),
                             paste0(paper_fig_location, 'Extra_CD8_depletion_figs/Frac_WT_Pre_Post_barcodes_notRNd15_or_RYm15.eps'),
                             paste0(paper_fig_location, 'Extra_CD8_depletion_figs/Frac_WT_Pre_Post_barcodes_FTY_animals.eps'),
                             paste0(paper_fig_location, 'Extra_CD8_depletion_figs/Frac_WT_Pre_Post_barcodes_Control_animals.eps')), 
           plot = WT_last_det_plot, width = 6, height = 5, units = "cm", device = "eps")
  }
  
  Frac_WT_test = wilcox.test(x = CD8depl_bcode_data_df_it$Frac_WT_last_det[CD8depl_bcode_data_df_it$Pre_Post == "Pre"],
              y = CD8depl_bcode_data_df_it$Frac_WT_last_det[CD8depl_bcode_data_df_it$Pre_Post == "Post"],
              paired = FALSE)
  
  print(Frac_WT_test$p.value)
  
  # Fraction WT Last detection cdf
  WT_last_det_cdf_plot = ggplot(CD8depl_bcode_data_df_it, aes(group = Pre_Post,
                                                          x = Frac_WT_last_det, 
                                                          color = Pre_Post,
                                                          fill = Pre_Post)) +
    stat_ecdf(geom = "step") +
    Barcode_timingScale_plotPalette_color + 
    labs(x = "Tat-SL8 Wild Type (%)", title = 
                             switch(kk, "All Animals", "Not RNd15 or RYm15")) + 
    theme(axis.text = element_text(size = paper_text_size),
          axis.title = element_text(size = paper_text_size))
  
  Post_depl_WT_stat_val = 0.5
  print(paste0(length(filter(CD8depl_bcode_data_df_it, 
                             ((Pre_Post == "Post")&
                                (Frac_WT_last_det >= Post_depl_WT_stat_val)))[,1]),
        " of ", length(filter(CD8depl_bcode_data_df_it, Pre_Post == "Post")[,1]),
        " barcodes first detected Post-CD8 depletion were greater than or equal to ",
        Post_depl_WT_stat_val*100, "% WT"))
  
  print(paste0(length(filter(CD8depl_bcode_data_df_it, 
                             ((Pre_Post == "Pre")&
                                (Frac_WT_last_det >= Post_depl_WT_stat_val)))[,1]),
        " of ", length(filter(CD8depl_bcode_data_df_it, Pre_Post == "Pre")[,1]),
        " barcodes first detected Pre-CD8 depletion were greater than or equal to ",
        Post_depl_WT_stat_val*100, "% WT"))
  
  print(WT_last_det_cdf_plot)
  
  #browser()
}



