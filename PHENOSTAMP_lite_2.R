# PHENOSTAMP_lite Projection Version 2.0
# Code refactored by Shafqat Ehsan for rerun Phenostamp_lite_2.R July 2025
# Original author of PHENOSTAMP: Bendict Anchang @NIEHS, 
# Modified PHENOSTAMP_LITE contributers: Bendict Anchang and Jiaqi Li


# project_root/
# ├── Phenostamp_lite_2.R    
# ├── stamp_compute.R                   ← helper script shipped with PHENOSTAMP
# ├── stamp_visualise.R                 ← helper script shipped with PHENOSTAMP
# ├── model/                            ← pre-trained reference objects
# │   ├── tsnedat.rdata
# │   ├── vor.rdata
# │   └── emt_nn.rdata
# ├── new_fcs/                      ← example set of FCS files (any name is acceptable)
# │   ├── sample1.fcs               <- add your fcs files here
# │   ├── sample2.fcs
# │   └── …
# └── …

#------------------------load model and packages
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

#Phenostamp_lite_function defined
PHENOSTAMP_lite <- 
  function(data_path, output_dir_name, colid, cex.size = 2){
    dir.create(output_dir_name, showWarnings = FALSE)
    
    fcs_files <- list.files(data_path, pattern = "\\.fcs$", full.names = TRUE)
    
    for (fcs_path in fcs_files) {
      sample_name <- sub("\\.[^.]*$", "", basename(fcs_path)) 
      print(paste("--> processing:", sample_name))
      
      # load fcs file and extract expression data
      fcs.data <- flowCore::read.FCS(fcs_path)
      dat <- flowCore::exprs(fcs.data)
      
      # get parameter metadata for plotting
      cluster_params <- flowCore::parameters(fcs.data)
      cluster_pd     <- pData(cluster_params)
      antibody <- cluster_pd$desc
      channel <- cluster_pd$name
      
      # subset and transform the 6 marker channels
      dat1 <- dat[, colid]
      colnames(dat1) <- names(colid)
      dat1 <- asinh(dat1)
      
      # define coordinate space from the pre-trained model
      window <- c(min(tsnedat[, 1])-5, max(tsnedat[, 1])+5, min(tsnedat[, 2])-5, max(tsnedat[, 2])+5)
      Y1 <- tsnedat[,1]
      Y2 <- tsnedat[,2]
      
      # normalize new data and predict coordinates
      data_scaled <- as.data.frame(scale(as.data.frame(dat1), center = apply(dat1, 2, min), scale = apply(dat1, 2, max) - apply(dat1, 2, min)))
      nnet.predict <- nnet:::predict.nnet(emt_nn, data_scaled)
      
      # rescale predicted coordinates to the original t-sne map
      pr.nn1 <- nnet.predict[,1] * (max(Y1) - min(Y1)) + min(Y1)
      pr.nn2 <- nnet.predict[,2] * (max(Y2) - min(Y2)) + min(Y2)
      sample_pred <- cbind(Y1 = pr.nn1, Y2 = pr.nn2)
      
      # save predicted coordinates
      sample_folder <- file.path(output_dir_name, sample_name)
      dir.create(sample_folder, showWarnings = FALSE)
      save(sample_pred, file = file.path(sample_folder, paste0(sample_name, ".rda")))
      print("    prediction saved, now plotting...")
      
      # set working directory for plotting
      curr_dir <- getwd()
      setwd(sample_folder)
      
      # plot projected density map
      density_filename <- paste0(sample_name, "_density.tiff")
      ccast_tsne_plot2c(sample_pred[,1], sample_pred[,2], vor, file1 = density_filename, sample = sample_name, window = window)
      
      # plot expression of each antibody on the map
      x <- sample_pred[,1]
      y <- sample_pred[,2]
      if (!is.null(antibody)){
        for (r in 1:length(antibody)) {
          marker_name <- antibody[r]
          if (is.na(marker_name)) next # skip channels without a description
          
          # create a valid filename and save plot
          valid_marker_name <- gsub("[^a-zA-Z0-9_.-]", "_", marker_name)
          tiff(paste0(valid_marker_name, ".tiff"), height = 15, width = 15, units = 'cm', res = 300)
          colvar <- asinh(dat[, channel[r]])
          scatter2D(x, y, colvar = colvar, pch = 20, cex = cex.size, 
                    main = marker_name, colkey = FALSE, xlim = window[1:2], ylim = window[3:4])
          plot(vor, wlines = "tess", add = TRUE, lwd = 2)
          dev.off()
        }
      }
      # return to original directory
      setwd(curr_dir)
    }
  }



##-Analyses Begins-------------------------
analysis_name <- "Naive_vs_Treated_Tarla_Analysis"#Name of output folder as well. 
input_fcs_folder <- "./new_fcs_data" # path to new fcs files

one_fcs_file_path <- file.path(input_fcs_folder, "pre_tarla.fcs") 
temp.data <-  flowCore::read.FCS(one_fcs_file_path)
temp.params <- flowCore::parameters(temp.data)
temp.pd     <- pData(temp.params)
temp.dat    <- exprs(temp.data)
View(temp.pd)
# Set the channel numbers for 6 markers

marker_channels <- c(
  ECadh    = 33,
  Vimentin = 31,
  CD44     = 32,
  CD24     = 37,
  Twist    = 25,
  MUC1     = 34
)

# run the analysis with the settings above
print(paste("starting analysis:", analysis_name))

PHENOSTAMP_lite(
  data_path = input_fcs_folder,
  output_dir_name = analysis_name,
  colid = marker_channels,
  cex.size = 2
)

print(paste("analysis complete. results in folder:", analysis_name))

