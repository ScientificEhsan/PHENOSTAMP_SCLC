# PHENOSTAMP_lite Projection Version 3.0
# Code refactored by Shafqat Ehsan for rerun Phenostamp_lite_2.R July 2025
# Phenostamp_lite_3.R August 2025, Percentages added to plot
# Original author of PHENOSTAMP: Bendict Anchang @NIEHS, 
# Modified PHENOSTAMP_LITE contributers: Bendict Anchang and Jiaqi Li

#MAKE SURE YOU HAVE THE FOLLOWING FOLDERS IN YOUR WORKING DIRECTORY IN THE FOLLOWING FORMAT--------------------------
# working_directory_/
# ├── Phenostamp_lite_3.R  
# ├── stamp_compute.R                   ← helper script shipped with PHENOSTAMP
# ├── stamp_visualise.R                 ← helper script shipped with PHENOSTAMP
# ├── model/                            ← pre-trained reference objects
# │   ├── tsnedat.rdata
# │   ├── vor.rdata
# │   └── emt_nn.rdata
# ├── input_fcs_folder/                 ← Name your analyses_folder anything. This will go as input_fcs_folder
# │   ├── sample1.fcs                   ← add your fcs files here and then use one of these as one_fcs_file_path
# │   ├── sample2.fcs
# │   └── …
# └── …

#load model and packages--------------------
libs <- c("graph", "fastcluster", "party", "plot3D", "MASS", "RColorBrewer", 
          "flowCore", "cluster", "plotly", "tripack", "deldir", "sp", "Rtsne")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
missing_pkgs <- libs[!(libs %in% installed.packages()[,"Package"])]
if (length(missing_pkgs) > 0) {
  BiocManager::install(missing_pkgs)
}
for (pkg in libs) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
if (!requireNamespace("bigvis", quietly = TRUE)) {
  devtools::install_github("hadley/bigvis")
}
all_libs <- c(libs, "bigvis")
sapply(all_libs, library, character.only = TRUE)

source("stamp_compute.R")
source("stamp_visualise.R")

load("model/tsnedat.rdata")
load("model/vor.rdata")
load("model/emt_nn.rdata")


#Helper function to define cell count, state, percentage-----------
#RUN AS IS NO CHANGES NEEDED
Cell_state_count <- function(x,y,vor){
  # Predicted cell coordinates
  pred_1=x
  pred_2=y
  # coordinates of the vor centers
  vor.center=cbind(vor$tri$x,vor.center=vor$tri$y)
  colnames(vor.center)=c("x","y")
  
  # define the cell state based on the location on the map for EMT map
  cell.state = c("E1","E2","E3","pEMT1","pEMT2","pEMT3","M","MET")
  n1=NULL
  n3=NULL
  # Assign each cell to a state (triangle) 
  
  for (j in 1:length(pred_1)){
    cell.xy <- data.frame(x=pred_1[j], y=pred_2[j])
    n2 <- which.min(unname(as.matrix(dist(rbind(cell.xy, as.matrix(vor.center)), diag=FALSE))[-1, 1])) 
    n3<-c(n3,cell.state[n2])
    n1=c(n1,n2)
  }
  # Cell states of each cells
  cell_states <- data.frame(pred_1 = pred_1, pred_2 = pred_2, cell_location = n1, cell_state = n3)
  
  ## State distribution table for sample projections
  MM=matrix(0,1,length(cell.state))
  rownames(MM)="Cell.States.count"
  colnames(MM)=cell.state
  MM[1,names(table(n3))]=table(n3)
  MM2=t(apply(MM,1,function(x) round((x/sum(x))*100,2)))
  rownames(MM2)="Cell.States.per"
  MM=rbind(MM,MM2)
  #
  return(list(cell_state = cell_states, cell_ratio = MM))
}

#PHENOSTAMP LITE DEFINED----------------
###RUN AS IS 
PHENOSTAMP_lite <- 
  function(data_path, output_dir_name = "output",
           ECadh, Vimentin, CD44, CD24, Twist, MUC1,
           cex.size = 2, percent_col = "blue"){
    
    # antibody channel vector is created inside the function, like the original
    colid <- c(ECadh = ECadh, Vimentin = Vimentin, CD44 = CD44, 
               CD24 = CD24,   Twist = Twist,     MUC1 = MUC1)
    
    # create an output folder
    dir.create(output_dir_name, showWarnings = FALSE)
    
    # list the fcs files
    fcs_files <- list.files(data_path, pattern = "\\.fcs$")
    
    # process files one by one
    for (fcs in fcs_files) {
      sample_name <- sub("\\.[^\\.]*$", "", fcs) 
      print(paste0("Start: ", sample_name))
      
      # load one fcs file
      fcs.data <- read.FCS(file.path(data_path, fcs))
      dat <- exprs(fcs.data) # Extract the expression data
      
      # Safety check for empty/single-cell files
      if (is.null(dat) || nrow(dat) < 2) {
        warning(paste("Skipping file", sample_name, "because it contains fewer than 2 cells."), call. = FALSE)
        next 
      }
      
      # Get all FACS parameters
      cluster_params <- flowCore::parameters(fcs.data)
      cluster_pd     <- pData(cluster_params)
      antibody <- cluster_pd$desc
      channel <- cluster_pd$name
      
      # Get the data for 6 markers
      dat1 <- dat[, colid, drop = FALSE] 
      colnames(dat1) <- names(colid)
      dat1 <- asinh(dat1)
      
      #'@-----Predicting-emt-states
      window <- c(min(tsnedat[, 1])-5, max(tsnedat[, 1])+5, min(tsnedat[, 2])-5, max(tsnedat[, 2])+5)
      Y1 <- tsnedat[,1]
      Y2 <- tsnedat[,2]
      
      # Normalize the data using the original method
      data <- as.data.frame(dat1)
      maxs <- apply(data, 2, max)
      mins <- apply(data, 2, min)
      scaled <- as.data.frame(scale(data, center = mins, scale = maxs - mins))
      scaled[is.na(scaled)] <- 0 # Add stability for zero-variance columns
      
      # Prediction
      nnet.predict <- nnet:::predict.nnet(emt_nn, scaled)
      pr.nn1 <- nnet.predict[,1]
      pr.nn2 <- nnet.predict[,2]
      
      pr.nn11 <- pr.nn1*(max(Y1)-min(Y1))+min(Y1)
      pr.nn22 <- pr.nn2*(max(Y2)-min(Y2))+min(Y2)
      
      sample_pred <- cbind(Y1=pr.nn11, Y2=pr.nn22)
      colnames(sample_pred) <- c("Y1", "Y2")
      
      # save the predicted result to a sample-specific folder
      sample_folder <- file.path(output_dir_name, sample_name)
      dir.create(sample_folder, showWarnings = FALSE)
      
      # Assign cell states and calculate percentage
      sample_cell_state <- Cell_state_count(sample_pred[,1], sample_pred[,2], vor)
      cell_states_ratio <- sample_cell_state$cell_ratio
      
      save(sample_pred, sample_cell_state, file=file.path(sample_folder, paste0(sample_name, ".rda")))
      print(paste0("Finished: Prediction of  ", sample_name))
      
      #'@-----Ploting
      curr_dir <- getwd()
      setwd(sample_folder)
      
      #'@--Ploting-density
      sample_filename <- paste0(sample_name, "_density.tiff")
      ccast_tsne_plot2c(sample_pred[,1], sample_pred[,2], vor, file1=sample_filename, sample=sample_name, window=window)
      
      #'@--Ploting-density-with-percentage
      vor.center <- cbind(vor$tri$x, vor$tri$y)
      sample_filename2 <- paste0(sample_name, "_percentage.tiff")
      tiff(sample_filename2, height = 12, width = 12, units = 'cm', pointsize = 6, compression = "lzw", type="cairo", res=300)
      
      cex <- 1.8
      xydens <- MASS::kde2d(sample_pred[,1], sample_pred[,2], n=100)
      plot(sample_pred[,1], sample_pred[,2], main = sample_name,
           xlab = "t-SNE1", ylab = "t-SNE2", ylim = window[3:4], xlim = window[1:2],
           cex.axis = cex, cex.lab = cex, cex.main = cex, cex=0.2, pch=20, col="grey")
      plot3D::contour2D(z=xydens$z, x=xydens$x, y=xydens$y, add=TRUE, lwd = 2, colkey = TRUE)
      plot(vor, add=TRUE, lwd=2)
      
      
      # <<< This is the updated block with new label positions >>>
      text(vor.center[1,1], vor.center[1,2] - 5, labels = paste0(as.character(cell_states_ratio["Cell.States.per", 2]), "%"), cex = 1.2, col = percent_col, bg = NA)
      text(vor.center[2,1] + 15, vor.center[2,2] - 15, labels = paste0(as.character(cell_states_ratio["Cell.States.per", 1]), "%"), cex = 1.2, col = percent_col, bg = NA)
      text(vor.center[3,1] - 12, vor.center[3,2] - 8, labels = paste0(as.character(cell_states_ratio["Cell.States.per", 4]), "%"), cex = 1.2, col = percent_col, bg = NA)
      text(vor.center[4,1] - 8, vor.center[4,2] - 12, labels = paste0(as.character(cell_states_ratio["Cell.States.per", 3]), "%"), cex = 1.2, col = percent_col, bg = NA)
      text(vor.center[5,1] - 8, vor.center[5,2], labels = paste0(as.character(cell_states_ratio["Cell.States.per", 5]), "%"), cex = 1.2, col = percent_col, bg = NA)
      text(vor.center[6,1], vor.center[6,2] + 15, labels = paste0(as.character(cell_states_ratio["Cell.States.per", 6]), "%"), cex = 1.2, col = percent_col, bg = NA)
      text(vor.center[7,1] + 15, vor.center[7,2] + 15, labels = paste0(as.character(cell_states_ratio["Cell.States.per", 7]), "%"), cex = 1.2, col = percent_col, bg = NA)
      text(vor.center[8,1] + 15, vor.center[8,2] + 8, labels = paste0(as.character(cell_states_ratio["Cell.States.per", 8]), "%"), cex = 1.2, col = percent_col, bg = NA)      
      dev.off()
      
      #'@--Ploting-each-antibody
      print("  creating individual antibody plots...")
      x <- sample_pred[,1]
      y <- sample_pred[,2]
      if (!is.null(antibody)){
        for ( r in 1:length(antibody)) {
          marker_name <- antibody[r]
          if(is.na(marker_name) || nchar(marker_name) == 0) next
          
          valid_marker_name <- gsub("[^a-zA-Z0-9_.-]", "_", marker_name)
          tiff(paste0(valid_marker_name, ".tiff"), height = 15, width = 15, units = 'cm', pointsize = 6, compression = "lzw",type="cairo", res=300)
          
          z <- dat[,channel[r]]
          colvar <- asinh(z)
          
          scatter2D(x, y, colvar=colvar, pch=20, cex=cex.size, main=marker_name, colkey = FALSE, xlim=window[1:2], ylim=window[3:4])
          plot(vor, wlines="tess", add=TRUE, lwd=2)
          dev.off()
        }
      }
      
      setwd(curr_dir)
      print(paste("--> completed:", sample_name))
    }
  }

##-Analyses Begins-------------------------
#YOU NEED TO UPDATE THE FOLLOWING 4 LINES AND FOLLOW INSTRUCTIONS
analysis_name <- "Naive_treated_CTCs_mapped" #Add your analyses name
input_fcs_folder <- "./Naive_Treated_FCS/" #Add directory where all your input files are

analysis_name <- "Tarla_treatment_CTC_mapped" #Add your analyses name
input_fcs_folder <- "./Tarla_Treatment_FCS/" #Add directory where all your input files are

analysis_name <- "Subtype_CTCs_Mapped" #Add your analyses name
input_fcs_folder <- "./Subtype_FCS/" #Add directory where all your input files are

analysis_name <- "Patient_level_CTCs_mapped" #Add your analyses name
input_fcs_folder <- "./sample_fcs/" #Add directory where all your input files are



one_fcs_file_path <- file.path(input_fcs_folder, "naive.fcs") #add one fcs file sample
one_fcs_file_path <- file.path(input_fcs_folder, "post_tarla.fcs") 
one_fcs_file_path <- file.path(input_fcs_folder, "A.fcs") 
one_fcs_file_path <- file.path(input_fcs_folder, "SC162-6.fcs") 

# This part for viewing channels remains the same
temp.data <-  flowCore::read.FCS(one_fcs_file_path)
temp.params <- flowCore::parameters(temp.data)
temp.pd     <- pData(temp.params) 

View(temp.pd) 
#check the temp.pd for the numbers   ECadh= 33,Vimentin = 31, CD44= 32,CD24 = 37,Twist= 25,MUC1= 34, 
#update as needed in the following function if different

print(paste("starting analysis:", analysis_name))

PHENOSTAMP_lite(
  data_path = input_fcs_folder,
  output_dir_name = analysis_name,
  ECadh     = 33,
  Vimentin  = 31,
  CD44      = 32,
  CD24      = 37,
  Twist     = 25,
  MUC1      = 34,
  cex.size  = 2,
  percent_col = "darkgreen" 
)
print(paste("analysis complete. results in folder:", analysis_name))

