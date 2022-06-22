#' FISH PLOT REPRESENTATION
#'
#' This script allows the user to get a fish plot representation (one for MiXCR
#' and one for Vidjil) of all rearrangement frequencies obtained in the study.
#'
#'    * Just semicolon-separated value files (.csv) containing MiXCR and 
#'    Vidjil merge output are accepted.
#'    * `RColorBrewer`, `htmlwidgets` and `timescape` previous installation
#'    in working environment is required.
#'    * The script is consisted of four functions: get_fish.df, add_WT, 
#'    get_WT and get_fishplot.
#'    * Implemented in a snakemake pipeline.
#'
#'Author: Andrea Sanchez de la Cruz - 8.6.2022


library(RColorBrewer) #rearrangement colours setting
library(htmlwidgets) #HTML objects management
library(timescape) #create fishplot object


## SNAKEMAKE I/O ##

merge_df <- snakemake@input[["merge_df"]]
mixcr_fishplot <- snakemake@output[["mixcr_fishplot"]]
vidjil_fishplot <- snakemake@output[["vidjil_fishplot"]]

# log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

## FUNCTIONS ##

#' Get rearrangement data from csv file
#' 
#' @param file CSV file (patient_merged.csv)
#' @return A data frame containing rearrangement data
get_fish.df <- function(file){
  
  #load file information in a data frame, selecting specific columns
  summary.df <- read.csv(file = file, sep = ';'
                        )[ , c('Sample','Condition','Locus',
                               'VDJ','mFreq','vFreq')]
  
  #join Locus and VDJ values and add it to a new column
  summary.df$clone_id <- paste0(summary.df$Locus, summary.df$VDJ)
  
  #data frame subset (remove Locus and VDJ cols)
  fish.df <- summary.df[ ,c('Sample','Condition','clone_id','mFreq','vFreq')]
  
  #rename columns
  colnames(fish.df) <- c('patient_id','timepoint','clone_id',
                         'Mclonal_prev','Vclonal_prev')
  
  #replace '-' values to zero values in MiXCR and Vidjil freqs
  fish.df[, c('Mclonal_prev','Vclonal_prev')] <- replace(
    fish.df[, c('Mclonal_prev','Vclonal_prev')], 
    fish.df[, c('Mclonal_prev','Vclonal_prev')] == '-',0)
  
  #convert freq cols from character to numeric
  fish.df$Mclonal_prev <- as.numeric(fish.df$Mclonal_prev)
  fish.df$Vclonal_prev <- as.numeric(fish.df$Vclonal_prev)
  
  #convert timepoint column to factor (BM < PPC < cfDNA)
  if ('PPC' %in% fish.df$timepoint){
    fish.df$timepoint <- factor(fish.df$timepoint, 
                                levels = c('BM', 'PPC', 'cfDNA'))
  } else {fish.df$timepoint <- factor(fish.df$timepoint)}
  
  #sort rows according to timepoint
  fish.df <- fish.df[order(fish.df$timepoint),]
  
  return(fish.df)
}


#' Add WT dummy rows to data frame
#' 
#' @param fish.df data frame with rearrangement information
#' @return The same data frame with additional WT rows
add_WT <- function(fish.df){
  
  #' Get WT freq vector for every timepoint
  #' 
  #' @param patient.df rearrangement data frame
  #' @param timepoint patient condition
  #' @return WT dummy row equal to fish.df
  get_WT <- function(patient.df, timepoint){

    #calculate WT dummy freq (1 - freq sum)
    M_wt <- 1 - sum(patient.df[patient.df$timepoint == timepoint,
                    'Mclonal_prev'])
    V_wt <- 1 - sum(patient.df[patient.df$timepoint == timepoint,
                    'Vclonal_prev'])
    
    #create WT vector
    patient = unique(patient.df$patient_id)
    tm_WT <- c(patient, as.character(timepoint), 'WT', M_wt, V_wt)
    
    return(tm_WT)
  }

  #get timepoints
  tm <- unique(fish.df$timepoint)
  
  #get WT frequencies for every timepoint
  WT_cols <- vapply(tm, function(x) get_WT(fish.df, x), character(5))
  #transpose to rows and convert to data frame
  WT.df <- as.data.frame(t(WT_cols))
  #rename to the same column names as fish.df
  names(WT.df) <- names(fish.df)
  
  #concatenate both data frames by row
  fish.df <- rbind(fish.df, WT.df)
  
  #convert freq columns to numeric
  fish.df$Mclonal_prev <- as.numeric(fish.df$Mclonal_prev)
  fish.df$Vclonal_prev <- as.numeric(fish.df$Vclonal_prev)
  
  return(fish.df)
}


#' Create fish plot generating specific required data frames
#' 
#' @param WT.fish.df data frame with rearrangement information and dummy WT
#' @param tool 'MiXCR' or 'Vidjil'
#' @param patool HTML filename
#' @return HTML file containing the fish plot
get_fishplot <- function(WT.fish.df, tool, patool){
  
  #rearrangement unique names and patient id
  pat = unique(WT.fish.df$patient_id)
  rearrs <- unique(WT.fish.df$clone_id)

  #associate specific RcolorBrewer color for every rearrangement
  rearr.colours <- RColorBrewer::brewer.pal(8, 'Set2')[1:length(rearrs)-1]
  #set grey color for WT and append to vector
  back.colour <- rev(RColorBrewer::brewer.pal(9, 'Set3'))[1]
  rearr.colours <- c(rearr.colours, back.colour)
  
  #select tool frequencies and rename column to clonal_prev
  if (tool == 'MiXCR'){
    #remove Vidjil column
    tool.df <- WT.fish.df[, -which(names(WT.fish.df) == 'Vclonal_prev')]
    #rename MiXCR freq column
    names(tool.df)[names(tool.df) == 'Mclonal_prev'] <- 'clonal_prev'
    
  } else if (tool == 'Vidjil'){
    #remove MiXCR column
    tool.df <- WT.fish.df[, -which(names(WT.fish.df) == 'Mclonal_prev')]
    #rename Vidjil freq column
    names(tool.df)[names(tool.df) == 'Vclonal_prev'] <- 'clonal_prev'
  }
  
  #sum frequencies of all rearrangements except WT
  rfreq <- sum(tool.df[tool.df$clone_id != 'WT',]$clonal_prev)
  
  #if there are rearrangements (sum > 0)
  if (rfreq){
    #set rearrangement data frame (Condition, VDJ and Freq)
    clonal_prev <- tool.df[,c('timepoint','clone_id','clonal_prev')]
    
    #set color data frame (specific for every rearrangement)
    clone_colours <- data.frame('clone_id' = rearrs, 'colour' = rearr.colours)
    
    #set phylogeny tree data frame
    source <- rep('WT', length(rearrs)-1) #WT
    target <- rearrs[!rearrs %in% 'WT'] #rearrangements
    tree_edges <- data.frame(source, target)
    
    #create interactive fish plot with previous data frames
    fishplot <- timescape(clonal_prev = clonal_prev,
                          tree_edges = tree_edges, 
                          clone_colours = clone_colours,
                          xaxis_title = 'Condition', 
                          yaxis_title = sprintf(
                            'Rearrangement Prevalence (%s)', tool),
                          genotype_position = "centre",
                          phylogeny_title = pat,
                          width = 800,
                          height = 350)
    
    #save fish plot in HTML file
    saveWidget(fishplot, patool)
  } else {
      writeLines(sprintf('No rearrangements found using %s for patient %s', 
                          tool, pat), con = patool)}
}


## MAIN CODE ##

#load patient data frame from file
try({fish.df <- get_fish.df(merge_df)}, TRUE)

#add WT freq to patient data frame
try({WT.df <- add_WT(fish.df)}, TRUE)

#get fish plot for MiXCR
tryCatch({
  m.fishplot <- get_fishplot(WT.df, 'MiXCR', mixcr_fishplot)
  delete_dir <- gsub('.html', '_files', mixcr_fishplot)

  if (file.exists(delete_dir)){
    system(sprintf('rm -r %s', delete_dir))
    }
  }, 
  error = function(e){
    writeLines('No rearrangements found using MiXCR for this patient.',
                con = mixcr_fishplot)
                })

#get fish plot for Vidjil
tryCatch({
  v.fishplot <- get_fishplot(WT.df, 'Vidjil', vidjil_fishplot)
  delete_dir <- gsub('.html', '_files', vidjil_fishplot)

  if (file.exists(delete_dir)){
    system(sprintf('rm -r %s', delete_dir))
    }
  }, 
  error = function(e){
    writeLines('No rearrangements found using Vidjil for this patient.',
                con = vidjil_fishplot)
                })
