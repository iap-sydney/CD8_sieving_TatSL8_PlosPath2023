#Steffen Docken
#12-08-22
#Code to visualize timcourse of SL8 variant proportions

#based on code from Deborah Cromer

rm(list=ls())
graphics.off()

library(dplyr)
library(reshape2)
library(stringr)
library(ggplot2)
library(cowplot)

data_folder = 'Figures/'

paper_text_size = 10
paper_tick_size = 7
geom_text_size = 5/14*(paper_tick_size-2)

dirs = list(
  data = 'csvs/',
  plots = 'Animal_plots/'
)

animals = c('RJy13','RUh14')
files = c('_SL8_variant_ID.csv','_SL8_variant_prop.csv')
proportion_breaks= c(0,.2,.6,1)
proportion_break_labels = paste0(100*proportion_breaks,'%')

# Names for the variant structure
variant_ID_col_names = c('plot_order', 'SL8_ID', c(1:8))


# Max variant number in any animal
max_num_variants = 30

# height of grey rectangles (in variants)
grey_height = 5
white_height = 5
total_height = grey_height+white_height
grey_alpha=.2

# Read in the data
for (animal in animals){
  variant_ID_read = read.csv(paste0(data_folder, dirs$data,animal, files[1]))
  proportion_read = read.csv(paste0(data_folder, dirs$data,animal, files[2]))
  
  # Change col names of the variant structure
  colnames(variant_ID_read) = variant_ID_col_names
  variant_ID_read  = mutate(variant_ID_read, num_aa_mut = 
                              rowSums(variant_ID_read[,3:10]!=""))
  variant_ID_read$num_aa_mut[1] = 0 #indicating first variant is WT (despite all
  #aa are listed)
  
  # Melt variant structure into a column struct
  variant_info = melt(variant_ID_read, id.vars = c('plot_order', 'SL8_ID', 'num_aa_mut'), 
                      variable.name='site', value.name='aminoAcid')
  #creates new data frame that still has the variable 'SL8_ID' and makes 
  #each amino acid site data point its own row. The variable 'site' is added, 
  #which contains the amino acid position the data in each row corresponds to
  variant_info = filter(variant_info, aminoAcid!="")
  
  # Set up the proportion structure as a vertical structure (with the correct names)
  colnames(proportion_read)[1:2] = variant_ID_col_names[1:2]
  proportion_read  = mutate(proportion_read, num_aa_mut = 
                              rowSums(variant_ID_read[,3:10]!=""))
  proportion_read$num_aa_mut[1] = 0#indicating first variant is WT (despite all
  #aa are listed)
  
  proportion_info = melt(proportion_read, id.vars = c('plot_order', 'SL8_ID', 'num_aa_mut'), variable.name='day', value.name='proportion')
  # Include the day, and which ATI the day is from
  proportion_info = mutate(proportion_info,
                           # Make the variant numbers into padded strings, so only the ones present in this animal get printed
                           SL8_ID = str_pad(SL8_ID, 3, pad = '0'),
                           day = str_pad(str_split(proportion_info$day,'X',simplify = T)[,2],3,side='left',pad='0'),
                           ATI = case_when(as.numeric(day) < 21 ~ 'Primary',
                                           as.numeric(day) < 300 ~ 'ATI1',
                                           T~'ATI2') )
  #edit and add columns
  proportion_info$ATI = factor(proportion_info$ATI, levels = c('Primary', 'ATI1', 'ATI2'))
  # exclude NA values
  proportion_info = filter(proportion_info, !is.na(proportion)) #remove NA values
  
  
  # Now we do the plotting
  
  # maximum number of variants in this animal
  max_variant_current = length(unique(proportion_info$SL8_ID))
  # set up a structure that tells us how many days are in each ATI (for gray rectangles)
  ndays = proportion_info %>% group_by(ATI) %>% summarise(nx = length(unique(day)))
  #summary data frame that will contain number of days we have data for for each ATI
  
  proportion_info_plot = filter(proportion_info, num_aa_mut <2)
  variant_info_plot = filter(variant_info, num_aa_mut <2)
  
  # Make the plots
  
  variant_timecourse_plot_generator <- function(data_df, background_df){
    if("day" %in% variable.names(data_df)){#data is for SL8 proportions
      plot_it = ggplot(data_df, aes(x=day,y=plot_order-.5, 
                                    fill = 10^proportion)) + geom_tile() + 
        scale_fill_viridis_c(name = '% of VL', direction = -1, 
                             breaks=proportion_breaks, labels=proportion_break_labels, 
                             trans ='sqrt',limits=c(0,1))+
        # Put the facet labels outside and at the bottom
        theme(strip.placement = 'outside')+
        facet_wrap(~ATI, scales = 'free_x', strip.position = 'bottom')+
        labs(y='',title='', x='Day Post Innoculation')
        
    }else if("site" %in% variable.names(data_df)){#data is aa  mutations
      plot_it = ggplot(data_df, aes(x=site,y=plot_order-.5,
                                              label = aminoAcid, 
                                              color = aminoAcid)) + 
        geom_text(size = geom_text_size)+
        labs(y='',title=animal,x='')
    }
    
    plot_it = plot_it +
      geom_rect(data = background_df, aes(xmin = 0.5, xmax = nx+.5, 
                                ymin=min(white_height,max_variant_current), 
                                ymax=min(total_height,max_variant_current)), 
              fill='grey', alpha = .2, inherit.aes = FALSE)+
      geom_rect(data = background_df, aes(xmin = 0.5, xmax = nx+.5, 
                                  ymin=min(total_height+white_height,max_variant_current), 
                                  ymax=min(2*total_height,max_variant_current)), 
                fill='grey', alpha = .2, inherit.aes = FALSE)+
      geom_rect(data = background_df, aes(xmin = 0.5, xmax = nx+.5, 
                                  ymin=min(2*total_height+white_height,max_variant_current), 
                                  ymax=min(3*total_height,max_variant_current)), 
                fill='grey', alpha = .2, inherit.aes = FALSE)+
      geom_rect(data = background_df, aes(xmin = 0.5, xmax = nx+.5, 
                                  ymin=min(3*total_height+white_height,max_variant_current), 
                                  ymax=min(4*total_height,max_variant_current)), 
                fill='grey', alpha = .2, inherit.aes = FALSE)+
      geom_rect(data = background_df, aes(xmin = 0.5, xmax = nx+.5, 
                                  ymin=min(4*total_height+white_height,max_variant_current), 
                                  ymax=min(5*total_height,max_variant_current)), 
                fill='grey', alpha = .2, inherit.aes = FALSE) +
      theme_classic(base_size = paper_text_size)+
      # Remove axis ticks / lines
      theme(axis.text.y = element_blank(),axis.ticks = element_blank(), axis.line=element_blank())+
      coord_cartesian(ylim=c(0,max_num_variants))
    
    
    if("day" %in% variable.names(data_df)){#data is for SL8 proportions
      plot_it = plot_it + geom_tile()+ 
        theme(legend.key.width = unit(.25,"cm"),
              axis.text.x = element_text(size = paper_tick_size,
                                         angle = 45))
    }else if("site" %in% variable.names(data_df)){#data is aa  mutations
      plot_it = plot_it + geom_text(size = geom_text_size) + 
        theme(axis.text.x = element_blank(),
              plot.title = element_text(size = paper_text_size))
    }
    
    return(plot_it)
    
  }
  
  naa_df = data.frame(xvar= "SL8aa", nx = 8)
  
  Proportion_Heatmap = variant_timecourse_plot_generator(proportion_info_plot, 
                                                         ndays)
  
  
  Variant_aa_mut_plot = variant_timecourse_plot_generator(variant_info_plot,
                                                          naa_df)
  
  
  Full_plot_it = plot_grid(Variant_aa_mut_plot + theme(legend.position = 'none'),
                           Proportion_Heatmap, ncol = 2,
                           align = 'h', axis = 'b', rel_widths = c(1,4))
  
  ggsave(paste0(data_folder, dirs$plots, animal, 'Variant_Proportion_timecourse.pdf'), 
         Full_plot_it, width=10, height = 7, units = "cm")
  
  
}