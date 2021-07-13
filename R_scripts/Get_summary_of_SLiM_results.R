#!/user/bin/Rscript

options(warning.length=5000)
suppressMessages(suppressPackageStartupMessages(try(library(getopt, lib.loc="/home/rhorvath/Rpackages"))))

args <- commandArgs(TRUE)

spec = matrix(c(
  "help"                          ,  "h"  , 0, "logical",
  "Patern"                        ,  "P"  , 1, "character",
  "Method"                        ,  "M"  , 2, "character"
), byrow = TRUE, ncol = 4)

opt = getopt(spec)

## Print a help message
if (!is.null(opt$help) ){
  write(
    "\n\nScript to summarize the SLiM results. 
    
    # requires
    -P a pattern that identifies all SLiM result files to use [required argument]. (e.g. neg_sel_s_0.005_TE_burst)
    -M which analyses to run [default All]. If mulptiple analyses should be run use - (e.g. SFS-AFS). Possible analyses: All, SFS, AFS, AgeDist, DeltaFreq.

    \n", stderr())
  
  write(getopt(spec, usage = TRUE), stderr())
  q(status=1)
}

write("############################\n#####\ \ \ \ \ Starting\ \ \ \ \ #####\n############################\n", stderr())

####################################################################################################################################################################################################
###########################                              Test inputs and set defaults                                         ######################################################################
####################################################################################################################################################################################################

if (is.null(opt$Patern) ) {write("#####\ \ \ \ \ Error\nNo patern provided (-P)", stderr()); stop() }
if (is.null(opt$Method) ) {opt$Method <- "All" } else {opt$Method <- unlist(strsplit(opt$Method, split="-"))}

####################################################################################################################################################################################################
###########################                              The script                                                 ################################################################################
####################################################################################################################################################################################################


## Define functions:
function_get_SNP_SFS <- function(x){
  my_raw_data_SNP <- read.table(x, header = TRUE, sep = ",")
  my_data <- my_raw_data_SNP[my_raw_data_SNP[,3] != 1,]
  
  my_SFS <- hist(my_data[,3], breaks = seq(0,1,0.05), plot = FALSE)$counts/length(my_data[,3]) 
  
  return(my_SFS)
}
function_TE_SFS <- function(y, z) {
  my_data_active_TE <- read.table(y, header = TRUE, sep = ",")
  my_data_disabled_TE <- read.table(z, header = TRUE, sep = ",")
  my_data_TE <- rbind(my_data_active_TE[,1:3], my_data_disabled_TE)
  my_data <- my_data_TE[my_data_TE[,3] != 1,]
  
  my_SFS <- hist(my_data[,3], breaks = seq(0,1,0.05), plot = FALSE)$counts/length(my_data[,3]) 
  
  return(my_SFS)
}
function_SNP_AFS <- function(x) {
  my_data <- read.table(x, header = TRUE, sep = ",")
  my_scaled_age <- my_data[,2]/5000
  my_AFS <- hist(my_scaled_age, breaks = seq(0, 2, 0.05), plot = FALSE)$counts/length(my_data[,2]) 
  
  return(c(my_AFS, length(my_scaled_age)))
}
function_TE_AFS <- function(x, y) {
  my_data_active <- read.table(x, header = TRUE, sep = ",")
  my_data_disabled <- read.table(y, header = TRUE, sep = ",")
  my_data <- rbind(my_data_active[,1:3], my_data_disabled)
  
  my_scaled_age <- my_data[my_data[,3] != 1,2]/5000
  my_AFS <- hist(my_scaled_age, breaks = seq(0, 2, 0.05), plot = FALSE)$counts/length(my_data[my_data[,3] != 1,2]) 
  
  return(c(my_AFS, length(my_scaled_age)))
}
function_my_downsampling <- function(a, s, t) {
  my_target_age_data <- s[s[,2] == a,]
  my_add <- 1
  while (length(my_target_age_data[,1]) == 0 ) {
    my_target_age_data <- s[s[,2] %in% c((a-my_add):(a+my_add)),]
    my_add <- my_add + 1
  }
  my_n <- length(t[t[,2] == a & t[,3] != 1 ,1])
  
  if (my_n <= length(my_target_age_data[,1])) {
    my_downsampled_data <- my_target_age_data[sample(1:length(my_target_age_data[,1]), my_n, replace = FALSE ), ]
  } else {
    my_downsampled_data <- my_target_age_data[sample(1:length(my_target_age_data[,1]), my_n, replace = TRUE ), ]
  }
  
  return(my_downsampled_data)
}
function_downsample_SNP <- function(x,y){
  all_TE_ages <- unique(y[y[,3] != 1,2])
  
  my_downsampled_frequencies <- sapply(all_TE_ages, function(e) { return(function_my_downsampling(e, x, y)[,3])} )
  
  return(my_downsampled_frequencies)
  
}
function_delta_frequ_age_bin <- function(x,y,z){
  my_data_SPS <- read.table(x, header = TRUE, sep = ",")
  my_data_active_TE <- read.table(y, header = TRUE, sep = ",")
  my_data_disabled_TE <- read.table(z, header = TRUE, sep = ",")
  my_data_TE <- rbind(my_data_active_TE[,1:3], my_data_disabled_TE)
  my_data_TE_var <- my_data_TE[my_data_TE[,3] != 1,]
  my_data_TE_ord <- my_data_TE_var[order(my_data_TE_var[,2]),]
  
  
  my_delta_frequ <- unlist(sapply(1:10, function(x){
    if (x != 10) {
      return(
        mean(my_data_TE_ord[((x-1)*floor(length(my_data_TE_ord[,2])/10) + 1):(floor(length(my_data_TE_ord[,2])/10)*x), 3]) - 
          mean(unlist(function_downsample_SNP(my_data_SPS, my_data_TE_ord[((x-1)*floor(length(my_data_TE_ord[,2])/10) + 1):(floor(length(my_data_TE_ord[,2])/10)*x), ])))
      )
    } else {
      return(
        mean(my_data_TE_ord[((x-1)*floor(length(my_data_TE_ord[,2])/10) + 1):length(my_data_TE_ord[,2]), 3]) - 
          mean(unlist(function_downsample_SNP(my_data_SPS, my_data_TE_ord[((x-1)*floor(length(my_data_TE_ord[,2])/10) + 1):length(my_data_TE_ord[,2]), ])))
      )
    }
  }))
  
  return(my_delta_frequ)
  
}
function_age_dist_SNP <- function(x) {
  my_data <- read.table(x, header = TRUE, sep = ",")
  my_mean_age <- rep(NA, length(seq(0, 0.95, 0.05)))
  my_length <- rep(NA, length(seq(0, 0.95, 0.05)))
  for (i in seq(0, 0.95, 0.05)) {
    my_mean_age[i*(1/0.05) + 1] <- mean(my_data[my_data[,3] > i & my_data[,3] <= i + 0.05 & my_data[,3] != 1, 2], na.rm = TRUE)
  }
  for (j in seq(0, 0.95, 0.05)) {
    my_length[j*(1/0.05) + 1] <- length(my_data[my_data[,3] > j & my_data[,3] <= j + 0.05 & my_data[,3] != 1, 2])
  }
  
  return(c(my_mean_age, my_length))
}
function_age_dist_TE <- function(x, y) {
  my_data_active <- read.table(x, header = TRUE, sep = ",")
  my_data_disabled <- read.table(y, header = TRUE, sep = ",")
  my_data <- rbind(my_data_active[,1:3], my_data_disabled)
  
  my_mean_age <- rep(NA, length(seq(0, 0.95, 0.05)))
  my_length <- rep(NA, length(seq(0, 0.95, 0.05)))
  for (i in seq(0, 0.95, 0.05)) {
    my_mean_age[i*(1/0.05) + 1] <- mean(my_data[my_data[,3] > i & my_data[,3] <= i + 0.05 & my_data[,3] != 1, 2], na.rm = TRUE)
  }
  for (j in seq(0, 0.95, 0.05)) {
    my_length[j*(1/0.05) + 1] <- length(my_data[my_data[,3] > j & my_data[,3] <= j + 0.05 & my_data[,3] != 1, 2])
  }
  
  return(c(my_mean_age, my_length))
}


## run:
## read in data:
SNP_g0_file <- list.files("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_secound_run", pattern = paste0("SNP_generation_0_", opt$Patern, "_run_number_"), full = TRUE)
SNP_g50_file <- list.files("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_secound_run", pattern = paste0("SNP_generation_50_", opt$Patern, "_run_number_"), full = TRUE)
SNP_g250_file <- list.files("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_secound_run", pattern = paste0("SNP_generation_250_", opt$Patern, "_run_number_"), full = TRUE)
SNP_g300_file <- list.files("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_secound_run", pattern = paste0("SNP_generation_300_", opt$Patern, "_run_number_"), full = TRUE)
SNP_g500_file <- list.files("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_secound_run", pattern = paste0("SNP_generation_500_", opt$Patern, "_run_number_"), full = TRUE)
SNP_g750_file <- list.files("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_secound_run", pattern = paste0("SNP_generation_750_", opt$Patern, "_run_number_"), full = TRUE)
SNP_g1250_file <- list.files("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_secound_run", pattern = paste0("SNP_generation_1250_", opt$Patern, "_run_number_"), full = TRUE)
SNP_g2000_file <- list.files("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_secound_run", pattern = paste0("SNP_generation_2000_", opt$Patern, "_run_number_"), full = TRUE)
SNP_g5000_file <- list.files("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_secound_run", pattern = paste0("SNP_generation_5000_", opt$Patern, "_run_number_"), full = TRUE)

TE_active_g0_file <- list.files("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_secound_run", pattern= paste0("TE.active_generation_0_", opt$Patern, "_run_number_"), full = TRUE)
TE_disabled_g0_file <- list.files("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_secound_run", pattern= paste0("TE.disabled_generation_0_", opt$Patern, "_run_number_"), full = TRUE)
TE_active_g50_file <- list.files("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_secound_run", pattern= paste0("TE.active_generation_50_", opt$Patern, "_run_number_"), full = TRUE)
TE_disabled_g50_file <- list.files("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_secound_run", pattern= paste0("TE.disabled_generation_50_", opt$Patern, "_run_number_"), full = TRUE)
TE_active_g250_file <- list.files("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_secound_run", pattern= paste0("TE.active_generation_250_", opt$Patern, "_run_number_"), full = TRUE)
TE_disabled_g250_file <- list.files("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_secound_run", pattern= paste0("TE.disabled_generation_250_", opt$Patern, "_run_number_"), full = TRUE)
TE_active_g300_file <- list.files("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_secound_run", pattern= paste0("TE.active_generation_300_", opt$Patern, "_run_number_"), full = TRUE)
TE_disabled_g300_file <- list.files("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_secound_run", pattern= paste0("TE.disabled_generation_300_", opt$Patern, "_run_number_"), full = TRUE)
TE_active_g500_file <- list.files("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_secound_run", pattern= paste0("TE.active_generation_500_", opt$Patern, "_run_number_"), full = TRUE)
TE_disabled_g500_file <- list.files("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_secound_run", pattern= paste0("TE.disabled_generation_500_", opt$Patern, "_run_number_"), full = TRUE)
TE_active_g750_file <- list.files("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_secound_run", pattern= paste0("TE.active_generation_750_", opt$Patern, "_run_number_"), full = TRUE)
TE_disabled_g750_file <- list.files("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_secound_run", pattern= paste0("TE.disabled_generation_750_", opt$Patern, "_run_number_"), full = TRUE)
TE_active_g1250_file <- list.files("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_secound_run", pattern= paste0("TE.active_generation_1250_", opt$Patern, "_run_number_"), full = TRUE)
TE_disabled_g1250_file <- list.files("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_secound_run", pattern= paste0("TE.disabled_generation_1250_", opt$Patern, "_run_number_"), full = TRUE)
TE_active_g2000_file <- list.files("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_secound_run", pattern= paste0("TE.active_generation_2000_", opt$Patern, "_run_number_"), full = TRUE)
TE_disabled_g2000_file <- list.files("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_secound_run", pattern= paste0("TE.disabled_generation_2000_", opt$Patern, "_run_number_"), full = TRUE)
TE_active_g5000_file <- list.files("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_secound_run", pattern= paste0("TE.active_generation_5000_", opt$Patern, "_run_number_"), full = TRUE)
TE_disabled_g5000_file <- list.files("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_secound_run", pattern= paste0("TE.disabled_generation_5000_", opt$Patern, "_run_number_"), full = TRUE)


if ("SFS" %in% opt$Method | "All" %in% opt$Method) {

## Get SNP SFS:
SNP_g0_SFS <- sapply(SNP_g0_file, function_get_SNP_SFS)
colnames(SNP_g0_SFS) <- 1:100
SNP_g50_SFS <- sapply(SNP_g50_file, function_get_SNP_SFS)
colnames(SNP_g50_SFS) <- 1:100
SNP_g250_SFS <- sapply(SNP_g250_file, function_get_SNP_SFS)
colnames(SNP_g250_SFS) <- 1:100
SNP_g300_SFS <- sapply(SNP_g300_file, function_get_SNP_SFS)
colnames(SNP_g300_SFS) <- 1:100
SNP_g500_SFS <- sapply(SNP_g500_file, function_get_SNP_SFS)
colnames(SNP_g500_SFS) <- 1:100
SNP_g750_SFS <- sapply(SNP_g750_file, function_get_SNP_SFS)
colnames(SNP_g750_SFS) <- 1:100
SNP_g1250_SFS <- sapply(SNP_g1250_file, function_get_SNP_SFS)
colnames(SNP_g1250_SFS) <- 1:100
SNP_g2000_SFS <- sapply(SNP_g2000_file, function_get_SNP_SFS)
colnames(SNP_g2000_SFS) <- 1:100
SNP_g5000_SFS <- sapply(SNP_g5000_file, function_get_SNP_SFS)
colnames(SNP_g5000_SFS) <- 1:100


## Get TE SFS:
TE_g0_SFS <- mapply(function_TE_SFS, TE_active_g0_file, TE_disabled_g0_file)
colnames(TE_g0_SFS) <- 1:100
TE_g50_SFS <- mapply(function_TE_SFS, TE_active_g50_file, TE_disabled_g50_file)
colnames(TE_g50_SFS) <- 1:100
TE_g250_SFS <- mapply(function_TE_SFS, TE_active_g250_file, TE_disabled_g250_file)
colnames(TE_g250_SFS) <- 1:100
TE_g300_SFS <- mapply(function_TE_SFS, TE_active_g300_file, TE_disabled_g300_file)
colnames(TE_g300_SFS) <- 1:100
TE_g500_SFS <- mapply(function_TE_SFS, TE_active_g500_file, TE_disabled_g500_file)
colnames(TE_g500_SFS) <- 1:100
TE_g750_SFS <- mapply(function_TE_SFS, TE_active_g750_file, TE_disabled_g750_file)
colnames(TE_g750_SFS) <- 1:100
TE_g1250_SFS <- mapply(function_TE_SFS, TE_active_g1250_file, TE_disabled_g1250_file)
colnames(TE_g1250_SFS) <- 1:100
TE_g2000_SFS <- mapply(function_TE_SFS, TE_active_g2000_file, TE_disabled_g2000_file)
colnames(TE_g2000_SFS) <- 1:100
TE_g5000_SFS <- mapply(function_TE_SFS, TE_active_g5000_file, TE_disabled_g5000_file)
colnames(TE_g5000_SFS) <- 1:100


## Write out results
write.table(SNP_g0_SFS, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/SNP_generation_0_", opt$Patern, "_SFS"), quote = FALSE)
write.table(SNP_g50_SFS, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/SNP_generation_50_", opt$Patern, "_SFS"), quote = FALSE)
write.table(SNP_g250_SFS, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/SNP_generation_250_", opt$Patern, "_SFS"), quote = FALSE)
write.table(SNP_g300_SFS, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/SNP_generation_300_", opt$Patern, "_SFS"), quote = FALSE)
write.table(SNP_g500_SFS, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/SNP_generation_500_", opt$Patern, "_SFS"), quote = FALSE)
write.table(SNP_g750_SFS, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/SNP_generation_750_", opt$Patern, "_SFS"), quote = FALSE)
write.table(SNP_g1250_SFS, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/SNP_generation_1250_", opt$Patern, "_SFS"), quote = FALSE)
write.table(SNP_g2000_SFS, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/SNP_generation_2000_", opt$Patern, "_SFS"), quote = FALSE)
write.table(SNP_g5000_SFS, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/SNP_generation_5000_", opt$Patern, "_SFS"), quote = FALSE)

write.table(TE_g0_SFS, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/TE_generation_0_", opt$Patern, "_SFS"), quote = FALSE)
write.table(TE_g50_SFS, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/TE_generation_50_", opt$Patern, "_SFS"), quote = FALSE)
write.table(TE_g250_SFS, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/TE_generation_250_", opt$Patern, "_SFS"), quote = FALSE)
write.table(TE_g300_SFS, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/TE_generation_300_", opt$Patern, "_SFS"), quote = FALSE)
write.table(TE_g500_SFS, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/TE_generation_500_", opt$Patern, "_SFS"), quote = FALSE)
write.table(TE_g750_SFS, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/TE_generation_750_", opt$Patern, "_SFS"), quote = FALSE)
write.table(TE_g1250_SFS, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/TE_generation_1250_", opt$Patern, "_SFS"), quote = FALSE)
write.table(TE_g2000_SFS, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/TE_generation_2000_", opt$Patern, "_SFS"), quote = FALSE)
write.table(TE_g5000_SFS, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/TE_generation_5000_", opt$Patern, "_SFS"), quote = FALSE)

}

if ("AFS" %in% opt$Method | "All" %in% opt$Method) {

## Get SNP AFS:
SNP_g0_AFS <- sapply(SNP_g0_file, function_SNP_AFS)
colnames(SNP_g0_AFS) <- 1:100
SNP_g50_AFS <- sapply(SNP_g50_file, function_SNP_AFS)
colnames(SNP_g50_AFS) <- 1:100
SNP_g250_AFS <- sapply(SNP_g250_file, function_SNP_AFS)
colnames(SNP_g250_AFS) <- 1:100
SNP_g300_AFS <- sapply(SNP_g300_file, function_SNP_AFS)
colnames(SNP_g300_AFS) <- 1:100
SNP_g500_AFS <- sapply(SNP_g500_file, function_SNP_AFS)
colnames(SNP_g500_AFS) <- 1:100
SNP_g750_AFS <- sapply(SNP_g750_file, function_SNP_AFS)
colnames(SNP_g750_AFS) <- 1:100
SNP_g1250_AFS <- sapply(SNP_g1250_file, function_SNP_AFS)
colnames(SNP_g1250_AFS) <- 1:100
SNP_g2000_AFS <- sapply(SNP_g2000_file, function_SNP_AFS)
colnames(SNP_g2000_AFS) <- 1:100
SNP_g5000_AFS <- sapply(SNP_g5000_file, function_SNP_AFS)
colnames(SNP_g5000_AFS) <- 1:100


## Get TE AFS:
TE_g0_AFS <- mapply(function_TE_AFS, TE_active_g0_file, TE_disabled_g0_file)
colnames(TE_g0_AFS) <- 1:100
TE_g50_AFS <- mapply(function_TE_AFS, TE_active_g50_file, TE_disabled_g50_file)
colnames(TE_g50_AFS) <- 1:100
TE_g250_AFS <- mapply(function_TE_AFS, TE_active_g250_file, TE_disabled_g250_file)
colnames(TE_g250_AFS) <- 1:100
TE_g300_AFS <- mapply(function_TE_AFS, TE_active_g300_file, TE_disabled_g300_file)
colnames(TE_g300_AFS) <- 1:100
TE_g500_AFS <- mapply(function_TE_AFS, TE_active_g500_file, TE_disabled_g500_file)
colnames(TE_g500_AFS) <- 1:100
TE_g750_AFS <- mapply(function_TE_AFS, TE_active_g750_file, TE_disabled_g750_file)
colnames(TE_g750_AFS) <- 1:100
TE_g1250_AFS <- mapply(function_TE_AFS, TE_active_g1250_file, TE_disabled_g1250_file)
colnames(TE_g1250_AFS) <- 1:100
TE_g2000_AFS <- mapply(function_TE_AFS, TE_active_g2000_file, TE_disabled_g2000_file)
colnames(TE_g2000_AFS) <- 1:100
TE_g5000_AFS <- mapply(function_TE_AFS, TE_active_g5000_file, TE_disabled_g5000_file)
colnames(TE_g5000_AFS) <- 1:100


## Write out results
write.table(SNP_g0_AFS, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/SNP_generation_0_", opt$Patern, "_AFS"), quote = FALSE)
write.table(SNP_g50_AFS, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/SNP_generation_50_", opt$Patern, "_AFS"), quote = FALSE)
write.table(SNP_g250_AFS, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/SNP_generation_250_", opt$Patern, "_AFS"), quote = FALSE)
write.table(SNP_g300_AFS, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/SNP_generation_300_", opt$Patern, "_AFS"), quote = FALSE)
write.table(SNP_g500_AFS, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/SNP_generation_500_", opt$Patern, "_AFS"), quote = FALSE)
write.table(SNP_g750_AFS, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/SNP_generation_750_", opt$Patern, "_AFS"), quote = FALSE)
write.table(SNP_g1250_AFS, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/SNP_generation_1250_", opt$Patern, "_AFS"), quote = FALSE)
write.table(SNP_g2000_AFS, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/SNP_generation_2000_", opt$Patern, "_AFS"), quote = FALSE)
write.table(SNP_g5000_AFS, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/SNP_generation_5000_", opt$Patern, "_AFS"), quote = FALSE)

write.table(TE_g0_AFS, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/TE_generation_0_", opt$Patern, "_AFS"), quote = FALSE)
write.table(TE_g50_AFS, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/TE_generation_50_", opt$Patern, "_AFS"), quote = FALSE)
write.table(TE_g250_AFS, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/TE_generation_250_", opt$Patern, "_AFS"), quote = FALSE)
write.table(TE_g300_AFS, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/TE_generation_300_", opt$Patern, "_AFS"), quote = FALSE)
write.table(TE_g500_AFS, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/TE_generation_500_", opt$Patern, "_AFS"), quote = FALSE)
write.table(TE_g750_AFS, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/TE_generation_750_", opt$Patern, "_AFS"), quote = FALSE)
write.table(TE_g1250_AFS, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/TE_generation_1250_", opt$Patern, "_AFS"), quote = FALSE)
write.table(TE_g2000_AFS, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/TE_generation_2000_", opt$Patern, "_AFS"), quote = FALSE)
write.table(TE_g5000_AFS, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/TE_generation_5000_", opt$Patern, "_AFS"), quote = FALSE)

}

if ("AgeDist" %in% opt$Method | "All" %in% opt$Method) {
  
  ## Get age and frequency distribution:
  SNP_g0_age_dist <- sapply(SNP_g0_file, function_age_dist_SNP)
  colnames(SNP_g0_age_dist) <- 1:100
  SNP_g50_age_dist <- sapply(SNP_g50_file, function_age_dist_SNP)
  colnames(SNP_g50_age_dist) <- 1:100
  SNP_g250_age_dist <- sapply(SNP_g250_file, function_age_dist_SNP)
  colnames(SNP_g250_age_dist) <- 1:100
  SNP_g300_age_dist <- sapply(SNP_g300_file, function_age_dist_SNP)
  colnames(SNP_g300_age_dist) <- 1:100
  SNP_g500_age_dist <- sapply(SNP_g500_file, function_age_dist_SNP)
  colnames(SNP_g500_age_dist) <- 1:100
  SNP_g750_age_dist <-sapply(SNP_g750_file, function_age_dist_SNP)
  colnames(SNP_g750_age_dist) <- 1:100
  SNP_g1250_age_dist <- sapply(SNP_g1250_file, function_age_dist_SNP)
  colnames(SNP_g1250_age_dist) <- 1:100
  SNP_g2000_age_dist <- sapply(SNP_g2000_file, function_age_dist_SNP)
  colnames(SNP_g2000_age_dist) <- 1:100
  SNP_g5000_age_dist <- sapply(SNP_g5000_file, function_age_dist_SNP)
  colnames(SNP_g5000_age_dist) <- 1:100
  
  TE_g0_age_dist <- mapply(function_age_dist_TE, TE_active_g0_file, TE_disabled_g0_file)
  colnames(TE_g0_age_dist) <- 1:100
  TE_g50_age_dist <- mapply(function_age_dist_TE, TE_active_g50_file, TE_disabled_g50_file)
  colnames(TE_g50_age_dist) <- 1:100
  TE_g250_age_dist <- mapply(function_age_dist_TE, TE_active_g250_file, TE_disabled_g250_file)
  colnames(TE_g250_age_dist) <- 1:100
  TE_g300_age_dist <- mapply(function_age_dist_TE, TE_active_g300_file, TE_disabled_g300_file)
  colnames(TE_g300_age_dist) <- 1:100
  TE_g500_age_dist <- mapply(function_age_dist_TE, TE_active_g500_file, TE_disabled_g500_file)
  colnames(TE_g500_age_dist) <- 1:100
  TE_g750_age_dist <- mapply(function_age_dist_TE, TE_active_g750_file, TE_disabled_g750_file)
  colnames(TE_g750_age_dist) <- 1:100
  TE_g1250_age_dist <- mapply(function_age_dist_TE, TE_active_g1250_file, TE_disabled_g1250_file)
  colnames(TE_g1250_age_dist) <- 1:100
  TE_g2000_age_dist <- mapply(function_age_dist_TE, TE_active_g2000_file, TE_disabled_g2000_file)
  colnames(TE_g2000_age_dist) <- 1:100
  TE_g5000_age_dist <- mapply(function_age_dist_TE, TE_active_g5000_file, TE_disabled_g5000_file)
  colnames(TE_g5000_age_dist) <- 1:100
  
  ## Write out results
  write.table(SNP_g0_age_dist, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/SNP_generation_0_", opt$Patern, "_age_dist"), quote = FALSE)
  write.table(SNP_g50_age_dist, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/SNP_generation_50_", opt$Patern, "_age_dist"), quote = FALSE)
  write.table(SNP_g250_age_dist, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/SNP_generation_250_", opt$Patern, "_age_dist"), quote = FALSE)
  write.table(SNP_g300_age_dist, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/SNP_generation_300_", opt$Patern, "_age_dist"), quote = FALSE)
  write.table(SNP_g500_age_dist, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/SNP_generation_500_", opt$Patern, "_age_dist"), quote = FALSE)
  write.table(SNP_g750_age_dist, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/SNP_generation_750_", opt$Patern, "_age_dist"), quote = FALSE)
  write.table(SNP_g1250_age_dist, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/SNP_generation_1250_", opt$Patern, "_age_dist"), quote = FALSE)
  write.table(SNP_g2000_age_dist, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/SNP_generation_2000_", opt$Patern, "_age_dist"), quote = FALSE)
  write.table(SNP_g5000_age_dist, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/SNP_generation_5000_", opt$Patern, "_age_dist"), quote = FALSE)
  
  write.table(TE_g0_age_dist, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/TE_generation_0_", opt$Patern, "_age_dist"), quote = FALSE)
  write.table(TE_g50_age_dist, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/TE_generation_50_", opt$Patern, "_age_dist"), quote = FALSE)
  write.table(TE_g250_age_dist, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/TE_generation_250_", opt$Patern, "_age_dist"), quote = FALSE)
  write.table(TE_g300_age_dist, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/TE_generation_300_", opt$Patern, "_age_dist"), quote = FALSE)
  write.table(TE_g500_age_dist, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/TE_generation_500_", opt$Patern, "_age_dist"), quote = FALSE)
  write.table(TE_g750_age_dist, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/TE_generation_750_", opt$Patern, "_age_dist"), quote = FALSE)
  write.table(TE_g1250_age_dist, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/TE_generation_1250_", opt$Patern, "_age_dist"), quote = FALSE)
  write.table(TE_g2000_age_dist, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/TE_generation_2000_", opt$Patern, "_age_dist"), quote = FALSE)
  write.table(TE_g5000_age_dist, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/TE_generation_5000_", opt$Patern, "_age_dist"), quote = FALSE)
  
}

if ("DeltaFreq" %in% opt$Method | "All" %in% opt$Method) {
  
  ## Get frequency delta based on age bin:
  g0_delta_frequ <- mapply(function_delta_frequ_age_bin, SNP_g0_file, TE_active_g0_file, TE_disabled_g0_file)
  colnames(g0_delta_frequ) <- 1:100
  g50_delta_frequ <- mapply(function_delta_frequ_age_bin, SNP_g50_file, TE_active_g50_file, TE_disabled_g50_file)
  colnames(g50_delta_frequ) <- 1:100
  g250_delta_frequ <- mapply(function_delta_frequ_age_bin, SNP_g250_file, TE_active_g250_file, TE_disabled_g250_file)
  colnames(g250_delta_frequ) <- 1:100
  g300_delta_frequ <- mapply(function_delta_frequ_age_bin, SNP_g300_file, TE_active_g300_file, TE_disabled_g300_file)
  colnames(g300_delta_frequ) <- 1:100
  g500_delta_frequ <- mapply(function_delta_frequ_age_bin, SNP_g500_file, TE_active_g500_file, TE_disabled_g500_file)
  colnames(g500_delta_frequ) <- 1:100
  g750_delta_frequ <- mapply(function_delta_frequ_age_bin, SNP_g750_file, TE_active_g750_file, TE_disabled_g750_file)
  colnames(g750_delta_frequ) <- 1:100
  g1250_delta_frequ <- mapply(function_delta_frequ_age_bin, SNP_g1250_file, TE_active_g1250_file, TE_disabled_g1250_file)
  colnames(g1250_delta_frequ) <- 1:100
  g2000_delta_frequ <- mapply(function_delta_frequ_age_bin, SNP_g2000_file, TE_active_g2000_file, TE_disabled_g2000_file)
  colnames(g2000_delta_frequ) <- 1:100
  g5000_delta_frequ <- mapply(function_delta_frequ_age_bin, SNP_g5000_file, TE_active_g5000_file, TE_disabled_g5000_file)
  colnames(g5000_delta_frequ) <- 1:100
  
  ## Write out results
  write.table(g0_delta_frequ, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/Generation_0_", opt$Patern, "_delta_frequ_0.1_quant_bin"), quote = FALSE)
  write.table(g50_delta_frequ, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/Generation_50_", opt$Patern, "_delta_frequ_0.1_quant_bin"), quote = FALSE)
  write.table(g250_delta_frequ, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/Generation_250_", opt$Patern, "_delta_frequ_0.1_quant_bin"), quote = FALSE)
  write.table(g300_delta_frequ, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/Generation_300_", opt$Patern, "_delta_frequ_0.1_quant_bin"), quote = FALSE)
  write.table(g500_delta_frequ, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/Generation_500_", opt$Patern, "_delta_frequ_0.1_quant_bin"), quote = FALSE)
  write.table(g750_delta_frequ, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/Generation_750_", opt$Patern, "_delta_frequ_0.1_quant_bin"), quote = FALSE)
  write.table(g1250_delta_frequ, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/Generation_1250_", opt$Patern, "_delta_frequ_0.1_quant_bin"), quote = FALSE)
  write.table(g2000_delta_frequ, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/Generation_2000_", opt$Patern, "_delta_frequ_0.1_quant_bin"), quote = FALSE)
  write.table(g5000_delta_frequ, file = paste0("/home/rhorvath/SLiM_TE_simulations/SLiM_TE_results_summaries_secound_run/Generation_5000_", opt$Patern, "_delta_frequ_0.1_quant_bin"), quote = FALSE)
  
}


