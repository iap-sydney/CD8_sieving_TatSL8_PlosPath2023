#Steffen Docken
#16-12-21
#This code will examine the levels of TL8-specific CD8 phenotypes in PBMC before
#ATI-1 and ATI-2

rm(list=ls())
graphics.off()


library(R.matlab)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(readxl)
library(stats)
library(cowplot)
library(latex2exp)
library(reshape2)
library(dplyr)

save_results = 1

#Importing data on percent of memory CD8s
d190_percentMem_data_raw <- read_excel("Data/d190-d373_SSD_names_corrected.xlsx", 
                                       sheet = "d190" , skip = 3, 
                                       col_types = c("text", "text","text","text",
                                                     "numeric", "numeric", "numeric", 
                                                     "numeric", "numeric", "numeric", 
                                                     "numeric", "numeric", "numeric", 
                                                     "numeric"),
                                       .name_repair =  "universal")
d373_percentMem_data_raw <- read_excel("Data/d190-d373_SSD_names_corrected.xlsx", 
                                       sheet = "d373", skip = 3,  
                                       col_types = c("text", "text","numeric","text",
                                                     "text", "numeric", "numeric", 
                                                     "numeric", "numeric", "numeric", 
                                                     "numeric", "numeric", "numeric", 
                                                     "numeric"),
                                       .name_repair =  "universal")
#read.csv

#Importing data on percent of total CD8s
percentTot_data_raw <- read_excel("Data/3-9-22 TL8-PBMC-LNMC_SSD_names_corrected.xlsx", 
                                  .name_repair =  "universal")


## creating data frames from raw data for samples from around day 190 and around 
# day 373
ART1_percent_df <- data.frame(File = d190_percentMem_data_raw$FCS.FILES,
                              Animal = d190_percentMem_data_raw$Monkey.ID,
                              Treatment = d190_percentMem_data_raw$Treatment,
                              Sample_type = d190_percentMem_data_raw$Sample,
                              dpi = d190_percentMem_data_raw$day.p.i,
                              ART = rep("ART_1", length(d190_percentMem_data_raw$FCS.FILES)),
                              TL8_percentTotal = rep(NA, length(d190_percentMem_data_raw$FCS.FILES)),
                              TL8_percentMem = d190_percentMem_data_raw$.TL8.w.mem.CD8,
                              CD69_percent_of_TL8Mem = d190_percentMem_data_raw$CD69,
                              CXCR3_percent_of_TL8Mem = d190_percentMem_data_raw$CXCR3,
                              CXCR5_percent_of_TL8Mem = d190_percentMem_data_raw$CXCR5,
                              Ki67_percent_of_TL8Mem = d190_percentMem_data_raw$KI67,
                              PD1_percent_of_TL8Mem = d190_percentMem_data_raw$PD1,
                              PerforinGranzyme_percent_of_TL8Mem = d190_percentMem_data_raw$Perf.Granz,
                              CM_percent_of_TL8Mem = d190_percentMem_data_raw$CM,
                              EM_percent_of_TL8Mem = d190_percentMem_data_raw$EM)

ART2_percent_df <- data.frame(File = d373_percentMem_data_raw$fcs.file.name,
                              Animal = d373_percentMem_data_raw$Monkey.ID,
                              Treatment = d373_percentMem_data_raw$Treatment,
                              Sample_type = d373_percentMem_data_raw$TISSUE,
                              dpi = d373_percentMem_data_raw$DAY.PI,
                              ART = rep("ART_2", length(d190_percentMem_data_raw$FCS.FILES)),
                              TL8_percentTotal = rep(NA, length(d373_percentMem_data_raw$fcs.file.name)),
                              TL8_percentMem = d373_percentMem_data_raw$.TL8.w.mem.CD8,
                              CD69_percent_of_TL8Mem = d373_percentMem_data_raw$CD69,
                              CXCR3_percent_of_TL8Mem = d373_percentMem_data_raw$CXCR3,
                              CXCR5_percent_of_TL8Mem = d373_percentMem_data_raw$CXCR5,
                              Ki67_percent_of_TL8Mem = d373_percentMem_data_raw$KI67,
                              PD1_percent_of_TL8Mem = d373_percentMem_data_raw$PD1,
                              PerforinGranzyme_percent_of_TL8Mem = d373_percentMem_data_raw$Perf.Granz,
                              CM_percent_of_TL8Mem = d373_percentMem_data_raw$CM,
                              EM_percent_of_TL8Mem = d373_percentMem_data_raw$EM)

for (ii in 1:length(percentTot_data_raw$Monkey.ID)) {
  if (abs(percentTot_data_raw$Days.P.i.[ii] -190) < 2) {
    ART1_percent_df$TL8_percentTotal[ART1_percent_df$Animal == percentTot_data_raw$Monkey.ID[ii]] =
      percentTot_data_raw$.TL8.Total.CD8[ii]
  }else if (abs(percentTot_data_raw$Days.P.i.[ii] -373) < 2) {
    ART2_percent_df$TL8_percentTotal[ART2_percent_df$Animal == percentTot_data_raw$Monkey.ID[ii]] =
      percentTot_data_raw$.TL8.Total.CD8[ii]
  }
}

PBMC_percent_df = rbind(ART1_percent_df, ART2_percent_df)
#combining data frames

Animals = unique(PBMC_percent_df$Animal)
Num_Animals = length(Animals)

phenotypes = setdiff(colnames(PBMC_percent_df), c("File", "Animal", "Treatment",
                                             "Sample_type", "dpi", "ART"))
num_phenotypes = length(phenotypes)

Results_total = data.frame(phenotype = phenotypes,
                           median_ART1 = rep(0, num_phenotypes),
                           median_ART2 = rep(0, num_phenotypes),
                           p_val = rep(0, num_phenotypes),
                           n =rep(0, num_phenotypes),
                           median_est = rep(0, num_phenotypes))
  
#looping through and doing tests on differences for each phenotype
for (jj in 1:num_phenotypes) {
  phenotype_data_df = PBMC_percent_df %>% select(Animal, ART, phenotypes[jj]) %>%
    dcast(Animal ~ ART, value.var = phenotypes[jj])
  #getting data for just this phenotype: selecting only animal, dpi, and 
  #current phenotype column. Then selecting only day ~190 and ~373. Then
  #seperating out day 190 and day 373 into different columns
  
  phenotype_data_df = phenotype_data_df[is.finite(phenotype_data_df$ART_1 - 
                                                    phenotype_data_df$ART_2),]
  #only including animals with both time points
  
  Results_total$median_ART1[Results_total$phenotype == phenotypes[jj]] = 
    median(phenotype_data_df$ART_1)
  Results_total$median_ART2[Results_total$phenotype == phenotypes[jj]] = 
    median(phenotype_data_df$ART_2)
  
  htest_it = wilcox.test(phenotype_data_df$ART_2, 
                         phenotype_data_df$ART_1, paired = TRUE,
                         alternative = "two.sided", conf.int = TRUE)
  
  Results_total$p_val[Results_total$phenotype == phenotypes[jj]] = htest_it$p.value
  Results_total$n[Results_total$phenotype == phenotypes[jj]] = length(phenotype_data_df$Animal)
  Results_total$median_est[Results_total$phenotype == phenotypes[jj]] = htest_it$estimate

}

if (save_results == 1){
  write.csv(Results_total, "ART1d190_v_ART2d373_PBMC_TL8CD8_level_tests.csv")
}
