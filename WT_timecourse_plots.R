#Steffen Docken
#10-05-21
#This code will visualize the time course data of WT SL8, VL, and total TL8 
#specific CD8s

rm(list=ls())
graphics.off()


library(R.matlab)
library(ggplot2)
library(ggridges)
library(RColorBrewer)
library(readxl)
library(cowplot)
library(latex2exp)
library(dplyr)

#color palettes:
plotPalette <- c(brewer.pal(n = 12,name = "Paired"), "#000000", "#999999") 
#black and grey are last two colors

plotPalette[1:2] = c("#bbffff","#00d8d8")

names(plotPalette) <- factor(c("unused", "unused", "Control", "FTY"))
TreatScale_plotPalette <- scale_colour_manual(name = "Treatment-Product",values = plotPalette,
                                                  breaks = c("Control", "FTY"),
                                                  labels = c("Control", "FTY720"),
                                                  aesthetics = c("color", "fill"),
                                                  limits = force)

plotPalette_w_DNA = c(plotPalette[1:2], plotPalette[1:4])
names(plotPalette_w_DNA) <- factor(c("Control.PBMC", "FTY.PBMC", "Control.LNMC", 
                               "FTY.LNMC", "Control.Plasma", "FTY.Plasma"))
TreatProdScale_plotPalette <- scale_colour_manual(name = "Treatment-Sample",values = plotPalette_w_DNA,
                                                breaks = c("Control.Plasma", "FTY.Plasma",
                                                           "Control.PBMC", "FTY.PBMC",
                                                           "Control.LNMC", "FTY.LNMC"),
                                                labels = c("Control-Plasma RNA", "FTY720-Plasma RNA",
                                                           "Control-PBMC DNA", "FTY720-PBMC DNA",
                                                           "Control-LNMC DNA", "FTY720-LNMC DNA"),
                                                aesthetics = c("color", "fill"),
                                                limits = force)

Sample_shape_def <- c("Control.Plasma" = 21, "FTY.Plasma" = 21,
                      "Control.PBMC" = 24, "FTY.PBMC" = 24,
                      "Control.LNMC" =22, "FTY.LNMC" =22)
Treat_shape_def <- c("FTY" = 21, "Control" = 24)
WT_det_fill_def <- c("1" = "black", "0" = "NA")

TreatProdScale_shape <- scale_shape_manual(name = "Treatment-Sample", 
                                           breaks = c("Control.Plasma", "FTY.Plasma",
                                                      "Control.PBMC", "FTY.PBMC",
                                                      "Control.LNMC", "FTY.LNMC"),
                                           labels = c("Control-Plasma RNA", "FTY720-Plasma RNA",
                                                      "Control-PBMC DNA", "FTY720-PBMC DNA",
                                                      "Control-LNMC DNA", "FTY720-LNMC DNA"), 
                                           values = Sample_shape_def)

names(plotPalette) <- factor(c("unused", "unused", "unused", "unused",
                              "FTY.ATI-2", "FTY.ATI-1",
                              "Allen et al. 2000.Primary", 
                              "Immonen et al. 2020.Primary",
                              "unused", "unused","unused", "unused",
                              "FTY.Pre-ART"))
StudyPhaseScale_plotPalette <- scale_colour_manual(name = "Study/Exp. Phase",values = plotPalette,
                                                  breaks = c("FTY.Pre-ART", "FTY.ATI-1", "FTY.ATI-2", "Immonen et al. 2020.Primary", "Allen et al. 2000.Primary"),
                                                  labels = c("Primary", "ATI 1", "ATI 2",
                                                             "Immonen et al. 2020","Allen et al. 2000"),
                                                  aesthetics = c("color", "fill"),
                                                  limits = force)
StudyPhaseScale_plotPalette_colorOnly <- scale_colour_manual(values = plotPalette,
                                                  breaks = c("FTY.Pre-ART", "FTY.ATI-1", "FTY.ATI-2", "Immonen et al. 2020.Primary", "Allen et al. 2000.Primary"),
                                                  labels = c("Primary", "ATI 1", "ATI 2",
                                                             "Immonen et al. 2020","Allen et al. 2000"),
                                                  aesthetics = "color",
                                                  limits = force)


ART_annotation = plotPalette[14]
ART_annotation_days = data.frame(xmin = c(14, 282), xmax = c(219, 374))
ATI1_annotation = plotPalette[6]
ATI1_annotation_days = data.frame(xmin = 219, xmax = 282)
ATI2_annotation = plotPalette[5]
ATI2_annotation_days = data.frame(xmin = 374, xmax = 465)
CD8depl_annotation = plotPalette[10]
CD8depl_annotation_days = data.frame(xmin = 465, xmax = 600)

background_alpha = 0.4

save_manuscript_fig_flag = 1; #0 = don't save manuscript figs. 1 = save
#manuscript figs

symbol_size = 1.75
text_size = 10

WT_data_from_MATLAB <- read.csv("WT_timecourse_data_for_R.csv")

VL_data_from_MATLAB <- read.csv("VL_timecourse_data_for_R.csv")

ART_days_from_MATLAB <- read.csv("ART_days_data_for_R.csv")

CD8depl_days_from_MATLAB <- read.csv("CD8depl_days_data_for_R.csv")

#Importing data on percent of total CD8s
percentTot_data_raw <- read.csv("Data/3-9-22 TL8-PBMC.csv")


WT_df = data.frame(WT_frac = WT_data_from_MATLAB$WT.Fraction.above.Background,
                   WT_det = WT_data_from_MATLAB$WT.detection,
                   WT_count = WT_data_from_MATLAB$WT.count,
                   WT_LOD = WT_data_from_MATLAB$WT.LOD,
                   Animal = WT_data_from_MATLAB$Animal,
                   Collection_Group = WT_data_from_MATLAB$Collection.Group,
                   FTY_treat = WT_data_from_MATLAB$FTY.Treatment,
                   Sample = WT_data_from_MATLAB$Sample,
                   dpi = WT_data_from_MATLAB$dpi,
                   phase = WT_data_from_MATLAB$phase,
                   timing = WT_data_from_MATLAB$timing,
                   AUC_VL = WT_data_from_MATLAB$AUC.VL.post.detection,
                   AUC_log10VL = WT_data_from_MATLAB$AUC.log10VL.post.detection,
                   CD8depl_dpi = rep(0,length(WT_data_from_MATLAB$dpi))) %>%
  mutate(Treat_prod = interaction(FTY_treat, Sample))

VL_df = data.frame(VL = VL_data_from_MATLAB$Viral.Load..copies.ml.,
                   VL_det = VL_data_from_MATLAB$Viral.Load.detection,
                   Animal = VL_data_from_MATLAB$Animal,
                   Collection_Group = VL_data_from_MATLAB$Collection.Group,
                   FTY_treat = VL_data_from_MATLAB$FTY.Treatment,
                   dpi = VL_data_from_MATLAB$dpi,
                   CD8depl_dpi = rep(0, length(VL_data_from_MATLAB$dpi)))

CD8depl_df = data.frame(Animal = CD8depl_days_from_MATLAB$Animal,
                        CD8depl_dpi = CD8depl_days_from_MATLAB$CD8depl_dpi)

for (ii in 1:length(WT_data_from_MATLAB$dpi)) {
  CD8depl_it = CD8depl_df$CD8depl_dpi[CD8depl_df$Animal == WT_df$Animal[ii]]
  
  if (CD8depl_it == 0) {
    WT_df$CD8depl_dpi[ii] = 1000 #setting CD8depl_dpi to a very late time because
  #animal was not CD8 depleted (listed as 0 in MATLAB data)
  } else {
    WT_df$CD8depl_dpi[ii] = CD8depl_it
  }
}

for (ii in 1:length(VL_df$dpi)) {
  CD8depl_it = CD8depl_df$CD8depl_dpi[CD8depl_df$Animal == VL_df$Animal[ii]]
  
  if (CD8depl_it == 0) {
    VL_df$CD8depl_dpi[ii] = 1000 #setting CD8depl_dpi to a very late time because
    #animal was not CD8 depleted (listed as 0 in MATLAB data)
  } else {
    VL_df$CD8depl_dpi[ii] = CD8depl_it
  }
}

## CD8 phenotype dynamics. Cells assayed for Tat-TL8, but this gives an 
#equivalent response to Tat-SL8
totCD8_percent_TL8_PBMC_df <- data.frame(Animal = percentTot_data_raw$Monkey.ID,
                                         FTY_treat = percentTot_data_raw$Treatment,
                                         Sample_type = rep("PBMC", length(percentTot_data_raw$Monkey.ID)),
                                         dpi = percentTot_data_raw$Days.P.i.,
                                         TL8_percentTotal = percentTot_data_raw$X.TL8.Total.CD8)

totCD8_percent_TL8_PBMC_df$FTY_treat[totCD8_percent_TL8_PBMC_df$FTY_treat == "FTY720"] = 
  "FTY"

dpi_ticks = seq(0, 500, 50)
dpi_labels = c("0", "", "100", "", "200", "", "300", "", "400", "", "500")
perc_WT_linear_ticks = seq(0, 1, .25)
perc_WT_linear_labels = c(TeX(r"($0\%$)"), TeX(r"($25\%$)"), TeX(r"($50\%$)"),
                          TeX(r"($75\%$)"), TeX(r"($100\%$)"))
perc_WT_log_ticks = c(0, .001, .01, .1, 1)
perc_WT_log_labels = c(TeX(r"(LOD)"), TeX(r"($0.1\%$)"), TeX(r"($1\%$)"),
                       TeX(r"($10\%$)"), TeX(r"($100\%$)"))
timing_ticks = seq(0, 80, 20)
timing_labels = timing_ticks
rep_ticks = seq(0, 60, 20)
rep_labels = rep_ticks


#WT figure generation
WT_protocol_plot = ggplot(WT_df, 
                          aes(x = dpi, y = WT_frac, color = Treat_prod,
                              shape = Treat_prod)) +
  annotate("rect", xmin = ART_annotation_days$xmin, xmax = ART_annotation_days$xmax, ymin = 0, ymax = Inf, 
           fill = ART_annotation, color = NA, alpha = background_alpha)+
  annotate("rect", xmin = ATI1_annotation_days$xmin, xmax = ATI1_annotation_days$xmax, ymin = 0, ymax = Inf, 
           fill = ATI1_annotation, color = NA, alpha = background_alpha)+
  annotate("rect", xmin = ATI2_annotation_days$xmin, xmax = ATI2_annotation_days$xmax, ymin = 0, ymax = Inf, 
           fill = ATI2_annotation, color = NA, alpha = background_alpha)+
  annotate("rect", xmin = CD8depl_annotation_days$xmin, xmax = CD8depl_annotation_days$xmax, ymin = 0, ymax = Inf, 
           fill = CD8depl_annotation, color = NA, alpha = background_alpha)+
  geom_point(data = WT_df[WT_df$WT_det == 1,], size = symbol_size, 
             color = "black", aes(fill = Treat_prod))+
  geom_point(data = WT_df[WT_df$WT_det == 0,], size = symbol_size) +
  TreatProdScale_plotPalette + 
  TreatProdScale_shape +
  scale_y_log10(breaks = perc_WT_log_ticks, labels = perc_WT_log_labels) +
  scale_x_continuous(breaks = dpi_ticks, labels = dpi_labels) + 
  labs(y = "Percent Wild-Type", x = "Days Post Infection")+ 
  theme_classic() +
  theme(title = element_text(size = text_size), 
        axis.text = element_text(size = text_size))

#dummy figures for legend generation

WT_protocol_plot_dummy2 = ggplot(WT_df, 
                                 aes(x = dpi, y = WT_frac,
                                     fill = factor(WT_det))) + 
  geom_point()+
  scale_fill_manual(values = WT_det_fill_def, name = "Detection", breaks  = c("1", "0"), labels = c("Detected", "Below LOD")) +
  scale_y_continuous(breaks = perc_WT_linear_ticks, labels = perc_WT_linear_labels) +
  scale_x_continuous(breaks = dpi_ticks, labels = dpi_labels) + 
  labs(y = "Percent Wild-Type", x = "Days Post Infection")+ theme_classic() +
  guides(fill=guide_legend(override.aes=list(shape=21)))

WT_protocol_plot_dummy3 = ggplot(WT_df, 
                                 aes(x = dpi, y = WT_frac, shape = Treat_prod,
                                     fill = Treat_prod)) +
  geom_point(size = symbol_size, color = "black")+
  TreatProdScale_plotPalette + TreatProdScale_shape + 
  scale_y_continuous(breaks = perc_WT_linear_ticks, labels = perc_WT_linear_labels) +
  scale_x_continuous(breaks = dpi_ticks, labels = dpi_labels) + 
  labs(y = "Percent Wild-Type", x = "Days Post Infection")+ theme_classic() +
  guides(col = guide_legend("Treatment-Sample"))


print(WT_protocol_plot)


VL_ticks = c(1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8)
VL_labels =  c(TeX(r"($10^{1}$)"), TeX(r"($10^{2}$)"), TeX(r"($10^{3}$)"),
               TeX(r"($10^{4}$)"), TeX(r"($10^{5}$)"), TeX(r"($10^{6}$)"),
               TeX(r"($10^{7}$)"), TeX(r"($10^{8}$)"))

#dummy plot for VL dynamics legend
VL_protocol_plot_dummy = ggplot(VL_df, aes(x = dpi, y = VL, group = Animal, color = FTY_treat))+ 
  annotate("rect", xmin = ART_annotation_days$xmin, xmax = ART_annotation_days$xmax, ymin = 0, ymax = Inf, 
           fill = ART_annotation, color = NA, alpha = background_alpha) + 
  annotate("rect", xmin = ATI1_annotation_days$xmin, xmax = ATI1_annotation_days$xmax, ymin = 0, ymax = Inf,
           fill = ATI1_annotation, color = NA, alpha = background_alpha)+
  annotate("rect", xmin = ATI2_annotation_days$xmin, xmax = ATI2_annotation_days$xmax, ymin = 0, ymax = Inf, 
           fill = ATI2_annotation, color = NA, alpha = background_alpha)+
  annotate("rect", xmin = CD8depl_annotation_days$xmin, xmax = CD8depl_annotation_days$xmax, ymin = 0, ymax = Inf, 
           fill = CD8depl_annotation, color = NA, alpha = background_alpha)+
  geom_line() +
  scale_y_log10(breaks = VL_ticks, labels = VL_labels) + 
  scale_x_continuous(breaks = dpi_ticks, labels = dpi_labels) + 
  labs(y = "Viral Load (copies/ml)", x = "Days Post Infection") + theme_classic()  +TreatScale_plotPalette

#VL dynamics figure
VL_protocol_plot = VL_protocol_plot_dummy +
  geom_point(data = VL_df[VL_df$VL_det == 1,], shape = 21, size = symbol_size, 
             color = "black", aes(fill = FTY_treat))+
  geom_point(data = VL_df[VL_df$VL_det == 0,], shape = 21, size = symbol_size)+
  theme(title = element_text(size = text_size), 
        axis.text = element_text(size = text_size))


print(VL_protocol_plot)

#Total Tat-SL8 specific CD8s figure
TotCD8_percent_TL8_PBMC_dynamics = ggplot(data = totCD8_percent_TL8_PBMC_df,
                                          aes(x = dpi, y = TL8_percentTotal, 
                                              group = Animal, color = FTY_treat, 
                                              fill = FTY_treat)) + #geom_boxplot(aes(x = dpi, y = TL8_percentTotal), color = "grey") +
  annotate("rect", xmin = ART_annotation_days$xmin, xmax = ART_annotation_days$xmax, ymin = 0, ymax = Inf, 
           fill = ART_annotation, color = NA, alpha = background_alpha) +
  annotate("rect", xmin = ATI1_annotation_days$xmin, xmax = ATI1_annotation_days$xmax, ymin = 0, ymax = Inf,
           fill = ATI1_annotation, color = NA, alpha = background_alpha)+
  annotate("rect", xmin = ATI2_annotation_days$xmin, xmax = ATI2_annotation_days$xmax, ymin = 0, ymax = Inf, 
           fill = ATI2_annotation, color = NA, alpha = background_alpha)+
  annotate("rect", xmin = CD8depl_annotation_days$xmin, xmax = CD8depl_annotation_days$xmax, ymin = 0, ymax = Inf, 
           fill = CD8depl_annotation, color = NA, alpha = background_alpha)+
  geom_line() + geom_point(shape = 21, color = "black", size = symbol_size) + TreatScale_plotPalette+
  scale_x_continuous(breaks = dpi_ticks, labels = dpi_labels) + 
  labs(y = "SL8-specific CD8\nT cells (%)", 
       x = "Days Post Infection") + theme_classic() +
  theme(title = element_text(size = text_size), 
        axis.text = element_text(size = text_size))

print(TotCD8_percent_TL8_PBMC_dynamics)

#Generation of full Figure 1
fig1 <- plot_grid(VL_protocol_plot + theme(legend.position = 'none'), 
          WT_protocol_plot + theme(legend.position = 'none'), 
          TotCD8_percent_TL8_PBMC_dynamics + theme(legend.position = 'none'),
          align = 'v', ncol = 1)

legend1 <- get_legend(WT_protocol_plot_dummy3 + theme(legend.box.margin = margin(0, 0, 0, 0),
                    title = element_text(size = text_size)))
legend2 <- get_legend(WT_protocol_plot_dummy2 + theme(legend.box.margin = margin(0, 0, 0, 0),
                    title = element_text(size = text_size)))
legends<-plot_grid(legend1, legend2, ncol = 2,align='h',axis='l',
                   rel_widths = c(2, 1))


print(legends)

if (save_manuscript_fig_flag == 1) {
  ggsave('Figures/VL_and_WT_vs_dpi.pdf', 
         fig1, width = 8.5, height = 15.5, units = "cm", dpi = 300)
  
  ggsave('Figures/VL_and_WT_vs_dpi_legend.eps', 
         legends, width = 8, height = 5, units = "cm", dpi = 300)
}

##VL around CD8 depletion
VL_ATI2_df = subset(VL_df, dpi > 370)
VL_CD8depl_df = subset(VL_ATI2_df, CD8depl_dpi != 1000) #not including animals that
#were not CD8 depleted

dpCD8depl_ticks = seq(-100, 100, 50)
dpCD8depl_labels = c("-100", "-50", "0", "50", "100")

#Dummy figure for VL at CD8 depletion figure
VL_CD8depl_plot_dummy = ggplot(VL_CD8depl_df, aes(x = dpi - CD8depl_dpi, y = VL, group = Animal)) + geom_line() +
  scale_y_log10(breaks = VL_ticks, labels = VL_labels) + 
  scale_x_continuous(breaks = dpCD8depl_ticks, labels = dpCD8depl_labels) + 
  labs(y = "Viral Load (copies/ml)", x = "Days Post-CD8 Depletion") + theme_classic()

#generating VL at CD8 depletion figure
VL_CD8depl_plot = VL_CD8depl_plot_dummy + 
  geom_point(data = VL_CD8depl_df[VL_CD8depl_df$VL_det == 1,], shape = 21, size = symbol_size, fill = "black")+
  geom_point(data = VL_CD8depl_df[VL_CD8depl_df$VL_det == 0,], shape = 21, size = symbol_size) +
  theme(title = element_text(size = text_size), 
        axis.text = element_text(size = text_size))


print(VL_CD8depl_plot)

VL_CD8depl_legend <- get_legend(WT_protocol_plot_dummy2 + 
                                  theme(legend.box.margin = margin(0, 0, 0, 0),
                                        title = element_blank(),
                                        legend.text = element_text(size = text_size)))

if (save_manuscript_fig_flag == 1) {
  ggsave('Figures/VL_at_CD8depl.eps', 
         VL_CD8depl_plot, width = 7, height = 4.5, units = "cm")
  
  ggsave('Figures/VL_at_CD8depl_legend.eps', 
         VL_CD8depl_legend, width = 2.7, height = 1.5, units = "cm")
}


## WT vs. days post reactivation and days of replication before CD8 depletion (Figure 2)
WT_plasma_FTY_df <- subset(subset(WT_df, Sample == "Plasma"),dpi <= CD8depl_dpi)
#FTY study animals
WT_plasma_FTY_df$Study = rep("FTY", length(WT_plasma_FTY_df$WT_frac)) #adding a column 
#for study name
WT_plasma_FTY_df$num_days_rep = rep(0, length(WT_plasma_FTY_df$WT_frac)) #adding a column
#for number of days of replication

### This code estimates the "total days of replication" by stitching ATIs onto
#start of treatment/end of previous ATI
num_WT_FTY_plasma_data = length(WT_plasma_FTY_df$Animal)

for (ii in 1:num_WT_FTY_plasma_data) {
  if (WT_plasma_FTY_df$phase[ii] == "Pre-ART") {
    WT_plasma_FTY_df$num_days_rep[ii] = WT_plasma_FTY_df$timing[ii]
  } else if (WT_plasma_FTY_df$phase[ii] == "ATI-1") {
    WT_plasma_FTY_df$num_days_rep[ii] = WT_plasma_FTY_df$timing[ii] + 14
  } else {
    WT_plasma_FTY_df$num_days_rep[ii] = NaN
  }
}

# compare WT dynamics in our experiment to Taina's data
WT_seq = "STPESANL" #sequence for WT to make sure pull in correct data from 
#Taina's results

TI_A13V019 = read_excel("Data/Immonen_etal_2020_SL8_evolution/Epitopes.xlsx", 
                        sheet = "A13V019", col_names = FALSE)
TI_A13V056 = read_excel("Data/Immonen_etal_2020_SL8_evolution/Epitopes.xlsx", 
                        sheet = "A13V056", col_names = FALSE)
TI_A13V136 = read_excel("Data/Immonen_etal_2020_SL8_evolution/Epitopes.xlsx", 
                        sheet = "A13V136", col_names = FALSE)
TI_A13V168 = read_excel("Data/Immonen_etal_2020_SL8_evolution/Epitopes.xlsx", 
                        sheet = "A13V168", col_names = FALSE)
TI_7811 = read_excel("Data/Immonen_etal_2020_SL8_evolution/Epitopes.xlsx", 
                        sheet = "7811", col_names = FALSE)
TI_25611 = read_excel("Data/Immonen_etal_2020_SL8_evolution/Epitopes.xlsx", 
                        sheet = "25611", col_names = FALSE)


TI_Animal_list = character()
TI_dpi_list = integer()
TI_WT_frac_list = double()

for (ii in 1:6) {
  Animal_name_it = switch(ii, "A13V019", "A13V056", "A13V136", "A13V168",
                          "7811", "25611")
  Data_it = switch(ii, TI_A13V019, TI_A13V056, TI_A13V136, TI_A13V168,
                   TI_7811, TI_25611)
  
  if (Data_it[2,1] != WT_seq) {
    print("Data not for WT sequenc")
  }
  
  for (jj in 2:length(Data_it[1,])) {
    TI_Animal_list = rbind(TI_Animal_list, Animal_name_it)
    TI_dpi_list = rbind(TI_dpi_list, as.numeric(Data_it[1,jj]))
    TI_WT_frac_list = rbind(TI_WT_frac_list, as.numeric(Data_it[2,jj]))
    #appending name of current animal, dpi and fraction WT for current time point
    
  }
}

Immonen_2020_WT = data.frame(Study = "Immonen et al. 2020",
                             Animal = TI_Animal_list, dpi = TI_dpi_list, 
                             WT_frac = TI_WT_frac_list, timing = TI_dpi_list,
                             num_days_rep = TI_dpi_list,
                             Sample = rep("Plasma", length(TI_dpi_list)),
                             phase = rep("Primary", length(TI_dpi_list)),
                             WT_det = rep(1, length(TI_dpi_list)),
                             FTY_treat = rep("Control", length(TI_dpi_list)),
                             CD8_depl = rep("not CD8 depl", length(TI_dpi_list)),
                             WT_count = rep(NA, length(TI_dpi_list)), #not including count of WT sequences 
                             WT_LOD = rep(NA, length(TI_dpi_list)), #not including LOD of WT 
                             AUC_VL = rep(NA, length(TI_dpi_list)), #not including AUC of VL
                             AUC_log10VL = rep(NA, length(TI_dpi_list)), #not including AUC of log10(VL)
                             Collection_Group = rep("Immonen et al. 2020", length(TI_dpi_list))) 
#creating dataframe of Immonen et al. 2020 data with columns to match those of 
#FTY study

# compare WT dynamics in our experiment to Allen et al. 2000's data (Fig 2C)
Allen_etal_2000_data = read.csv("Data/Allen_etal_2000/Allen_etal_2000_SL8_decay_data.csv")

Allen_etal_2000_WT = data.frame(Study = "Allen et al. 2000",
                             Animal = Allen_etal_2000_data$Animal,
                             dpi = Allen_etal_2000_data$Days,
                             WT_frac = pmax(Allen_etal_2000_data$WT.sequences/
                               Allen_etal_2000_data$Template.Number, 
                               1/Allen_etal_2000_data$Template.Number), #setting WT fraction to LOD when WT not detected. 
                             timing = Allen_etal_2000_data$Days,
                             num_days_rep = Allen_etal_2000_data$Days,
                             Sample = rep("Plasma", length(Allen_etal_2000_data$Days)),
                             phase = rep("Primary", length(Allen_etal_2000_data$Days)),
                             WT_det = as.numeric(Allen_etal_2000_data$WT.sequences>0),
                             FTY_treat = rep("Control", length(Allen_etal_2000_data$WT.sequences)),
                             CD8_depl = rep("not CD8 depl", length(Allen_etal_2000_data$WT.sequences)),
                             WT_count = Allen_etal_2000_data$WT.sequences,
                             WT_LOD = 1/Allen_etal_2000_data$Template.Number,
                             AUC_VL = rep(NA, length(Allen_etal_2000_data$Days)), #not including AUC of VL (would need to look if available)
                             AUC_log10VL = rep(NA, length(Allen_etal_2000_data$Days)), #not including AUC of log10(VL) (would need to look if available)
                             Collection_Group = rep("Allen et al. 2000", length(Allen_etal_2000_data$Days))) 
#creating dataframe of Allen et al. 2000 data with columns to match those of 
#FTY study

WT_plasma_df = rbind(WT_plasma_FTY_df[,c("Study", "Animal", "dpi", "phase", "timing",
                                     "num_days_rep", "WT_frac", "WT_det", "WT_count",
                                     "WT_LOD", "Collection_Group", "FTY_treat",
                                     "Sample", "AUC_VL", "AUC_log10VL")], 
                     Immonen_2020_WT[,c("Study", "Animal", "dpi", "phase", "timing",
                                        "num_days_rep", "WT_frac", "WT_det", "WT_count",
                                        "WT_LOD", "Collection_Group", "FTY_treat",
                                        "Sample", "AUC_VL", "AUC_log10VL")],
                     Allen_etal_2000_WT[,c("Study", "Animal", "dpi", "phase", "timing",
                                           "num_days_rep", "WT_frac", "WT_det", "WT_count",
                                           "WT_LOD", "Collection_Group", "FTY_treat",
                                           "Sample", "AUC_VL", "AUC_log10VL")]) #adding Immonen and 
#Allen data to WT data frame

# Generating plot of fraction WT vs. days post reactivation
WT_vs_dpReact = ggplot(WT_plasma_df[WT_plasma_df$Study == "FTY",], 
                       aes(x = timing, y = WT_frac, group = interaction(Animal, phase),
                           color = factor(interaction(Study, phase)), shape = FTY_treat)) + geom_line() +
  geom_point(data = WT_plasma_df[WT_plasma_df$WT_det == 1&WT_plasma_df$Study == "FTY",], 
             aes(fill = factor(interaction(Study, phase))), size = symbol_size) + 
  geom_point(data = WT_plasma_df[WT_plasma_df$WT_det == 0&WT_plasma_df$Study == "FTY",], size = symbol_size) + 
  scale_shape_manual(values = Treat_shape_def, name = "Treatment", 
                     labels =  c("FTY720", "Control"))+
  StudyPhaseScale_plotPalette +
  scale_y_log10(breaks = perc_WT_log_ticks, labels = perc_WT_log_labels) +
  scale_x_continuous(breaks = timing_ticks, labels = timing_labels) +
  labs(y = "Percent Wild-Type", x = "Day Post Infection/Reactivation") + theme_classic() +
  theme(title = element_text(size = text_size), 
        axis.text = element_text(size = text_size)) + xlim(0,80)

print(WT_vs_dpReact)

# Generating plot of fraction WT vs. total days of replication
WT_vs_rep_days = ggplot(WT_plasma_df[WT_plasma_df$phase != "ATI-2",], 
                        aes(x = num_days_rep, y = WT_frac, group = Animal,
                            color = factor(interaction(Study, phase)), shape = FTY_treat)) + geom_line() +
  geom_point(data = WT_plasma_df[WT_plasma_df$WT_det == 1&WT_plasma_df$phase != "ATI-2",], 
             aes(fill = factor(interaction(Study, phase))), size = symbol_size) + 
  geom_point(data = WT_plasma_df[WT_plasma_df$WT_det == 0&WT_plasma_df$phase != "ATI-2",], size = symbol_size) + 
  scale_shape_manual(values = Treat_shape_def, name = "Treatment")+
  StudyPhaseScale_plotPalette+
  scale_y_log10(breaks = perc_WT_log_ticks[3:5], labels = perc_WT_log_labels[3:5]) +
  scale_x_continuous(breaks = rep_ticks, labels = rep_labels) +
  labs(y = "Percent Wild-Type", x = "Total Days of Replication") + theme_classic() +
  theme(title = element_text(size = text_size), 
        axis.text = element_text(size = text_size))  + xlim(0,80)

print(WT_vs_rep_days)

#dummy figure for figure 2
legend1_Fig2_dummy = ggplot(data = filter(WT_plasma_df, Study == "FTY"), 
                        aes(x = timing, y = WT_frac, group = Animal,
                            color = factor(interaction(Study, phase)))) + geom_line() +
  geom_point(size = symbol_size) + StudyPhaseScale_plotPalette_colorOnly+ 
  theme_classic() + labs(color = "Exp. Phase") +
  theme(title = element_text(size = text_size), 
        axis.text = element_text(size = text_size))

legend2_Fig2_dummy = ggplot(WT_plasma_df[WT_plasma_df$Study == "FTY",], 
                            aes(x = timing, y = WT_frac, group = interaction(Animal, phase),
                                shape = FTY_treat)) + geom_line() +
  geom_point() + 
  scale_shape_manual(values = Treat_shape_def, name = "Treatment", 
                     labels = c("FTY720", "Control")) + 
  theme_classic()

legend3_Fig2_dummy = ggplot(data = filter(WT_plasma_df, Study != "FTY"), 
                            aes(x = timing, y = WT_frac, group = Animal,
                                color = factor(interaction(Study, phase)))) + geom_line() +
  geom_point(size = symbol_size) + StudyPhaseScale_plotPalette_colorOnly+ 
  theme_classic() + labs(color = "Previous Studies") +
  theme(title = element_text(size = text_size), 
        axis.text = element_text(size = text_size))

fig2 <- plot_grid(WT_vs_dpReact + theme(legend.position = 'none'), 
                  WT_vs_rep_days + theme(legend.position = 'none'), 
                  align = 'v', ncol = 1)

legend1_fig2 <- get_legend(legend1_Fig2_dummy + theme(legend.box.margin = margin(0, 0, 0, 0),
                                                              title = element_text(size = 8)))
legend2_fig2 <- get_legend(legend2_Fig2_dummy + theme(legend.box.margin = margin(0, 0, 0, 0),
                                                              title = element_text(size = 8)))
legend3_fig2 <- get_legend(legend3_Fig2_dummy + theme(legend.box.margin = margin(0, 0, 0, 0),
                                                      title = element_text(size = 8)))

legends_fig2<-plot_grid(legend1_fig2, legend2_fig2, legend3_fig2, legend2,
                        ncol = 1, align='v',axis='t', rel_heights = c(1,1,1,1))
fig2_full<-plot_grid(fig2,legends_fig2,
                     align = 'vh',
                     rel_widths = c(5, 2),
                     ncol=2)

print(fig2_full)


if (save_manuscript_fig_flag == 1) {
  ggsave('Figures/WT_vs_timing_and_Rep_days.eps', 
         fig2_full, width = 14, height = 10, units = "cm")
}
