#######################################################################################################
##################       Analyse the SLiM runs summary outputs           ##############################
##################
##################       File for Github                                 ##############################
#######################################################################################################


#######################################
#### Function needed for the analyses:

# expectad age at a specific freqeuncy (Kimura and Ohta 1973)
expected_allele_age <- function(x, N) {
  return( (-4 * N * x)/(1-x) * log(x))
}
# expectad age at a specific freqeuncy under selection (Maruyama 1974)
my_function_a <- function(x, S = 100){
  ((exp(S*x) - 1)*(exp(-S*x) - exp(-S)))/(x*(1-x))
}
my_function_b <- function(x, S = 100, y = 0.001){
  ((1 - exp(-S*(1-x)) )*(exp(-S*y) - exp(-S*x)) )/(x*(1-x))
}
my_maruyama_age <- function(S, N, y){
  return(
    ((4*N)/(S*(1 - exp(-S))))
    *integrate(my_function_a, S = S, 0, 1)$value 
    - ((4*N)/(S*(exp(-S*y) - exp(-S))))
    *integrate(my_function_b, S = S, y = y, y, 1)$value
  )
}
# get mean age at a specific freqeuncy
function_get_mean_age <- function(x) {
  my_res <- c()
  for (i in 1:20) {
    my_res <- c(my_res, mean(as.numeric(as.numeric(x[i,])), na.rm = TRUE))
  }
  return(my_res)
}
# Plot mean age
function_plot_mean_age <- function(x, ...) {
  boxplot(as.numeric(as.numeric(x[1,])), 
          as.numeric(as.numeric(x[2,])),
          as.numeric(as.numeric(x[3,])),
          as.numeric(as.numeric(x[4,])),
          as.numeric(as.numeric(x[5,])),
          as.numeric(as.numeric(x[6,])),
          as.numeric(as.numeric(x[7,])),
          as.numeric(as.numeric(x[8,])),
          as.numeric(as.numeric(x[9,])),
          as.numeric(as.numeric(x[10,])),
          as.numeric(as.numeric(x[11,])),
          as.numeric(as.numeric(x[12,])),
          as.numeric(as.numeric(x[13,])),
          as.numeric(as.numeric(x[14,])),
          as.numeric(as.numeric(x[15,])),
          as.numeric(as.numeric(x[16,])),
          as.numeric(as.numeric(x[17,])),
          as.numeric(as.numeric(x[18,])),
          as.numeric(as.numeric(x[19,])),
          as.numeric(as.numeric(x[20,])), xaxt="n", ...
  )
}
# plot expected delta age based on Maruyama (1974)
function_add_points_expected_delta_maruyama <- function(S = (4*500*-0.01), N = 500, N2 = 500, ...) {
  points(1, my_maruyama_age(S, N, 0.025) - expected_allele_age(0.025, N2), ...)
  points(2, my_maruyama_age(S, N, 0.075) - expected_allele_age(0.075, N2), ...)
  points(3, my_maruyama_age(S, N, 0.125) - expected_allele_age(0.125, N2), ...)
  points(4, my_maruyama_age(S, N, 0.175) - expected_allele_age(0.175, N2), ...)
  points(5, my_maruyama_age(S, N, 0.225) - expected_allele_age(0.225, N2), ...)
  points(6, my_maruyama_age(S, N, 0.275) - expected_allele_age(0.275, N2), ...)
  points(7, my_maruyama_age(S, N, 0.325) - expected_allele_age(0.325, N2), ...)
  points(8, my_maruyama_age(S, N, 0.375) - expected_allele_age(0.375, N2), ...)
  points(9, my_maruyama_age(S, N, 0.425) - expected_allele_age(0.425, N2), ...)
  points(10, my_maruyama_age(S, N, 0.475) - expected_allele_age(0.475, N2), ...)
  points(11, my_maruyama_age(S, N, 0.525) - expected_allele_age(0.525, N2), ...)
  points(12, my_maruyama_age(S, N, 0.575) - expected_allele_age(0.575, N2), ...)
  points(13, my_maruyama_age(S, N, 0.625) - expected_allele_age(0.625, N2), ...)
  points(14, my_maruyama_age(S, N, 0.675) - expected_allele_age(0.675, N2), ...)
  points(15, my_maruyama_age(S, N, 0.725) - expected_allele_age(0.725, N2), ...)
  points(16, my_maruyama_age(S, N, 0.775) - expected_allele_age(0.775, N2), ...)
  points(17, my_maruyama_age(S, N, 0.825) - expected_allele_age(0.825, N2), ...)
  points(18, my_maruyama_age(S, N, 0.875) - expected_allele_age(0.875, N2), ...)
  points(19, my_maruyama_age(S, N, 0.925) - expected_allele_age(0.925, N2), ...)
  points(20, my_maruyama_age(S, N, 0.975) - expected_allele_age(0.975, N2), ...)
}
# plot delta frequency
function_delta_frequ_plot <- function(x, my_xlab = "Age deciles", my_ylab = expression(paste(Delta, " frequency")), ...){
  boxplot(as.numeric(x[1,]),
          as.numeric(x[2,]),
          as.numeric(x[3,]),
          as.numeric(x[4,]),
          as.numeric(x[5,]),
          as.numeric(x[6,]),
          as.numeric(x[7,]),
          as.numeric(x[8,]),
          as.numeric(x[9,]),
          as.numeric(x[10,]), xlab = my_xlab, ylab = my_ylab, ...)
}
## spearman correlation test for negative correlation: 
my_get_number_of_sig <- function(x){
  my_p <- c()
  for (h in 1:100) {
    my_p <- c(my_p, 
              cor.test(as.numeric(x[1:10,h]), 1:10, method = "spearman", alternative = "less")$p.value
    )
  }
  return(sum(my_p < 0.05))
}


##########################
#### Read in data needed

# SFS (-M SFS)
# each column corresponds to one run 
# rows: proportion of sites at the freqeuncies 0.05, 0.10, 0.15, ..., respectively.
TE_generation_0_neutral_TE_burst_SFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_0_neutral_TE_burst_SFS")
TE_generation_50_neutral_TE_burst_SFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_50_neutral_TE_burst_SFS")
TE_generation_250_neutral_TE_burst_SFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_250_neutral_TE_burst_SFS")
TE_generation_300_neutral_TE_burst_SFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_300_neutral_TE_burst_SFS")
TE_generation_500_neutral_TE_burst_SFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_500_neutral_TE_burst_SFS")
TE_generation_750_neutral_TE_burst_SFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_750_neutral_TE_burst_SFS")
TE_generation_1250_neutral_TE_burst_SFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_1250_neutral_TE_burst_SFS")
TE_generation_2000_neutral_TE_burst_SFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_2000_neutral_TE_burst_SFS")
TE_generation_5000_neutral_TE_burst_SFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_5000_neutral_TE_burst_SFS")
SNP_generation_0_neutral_bottleneck_SFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/SNP_generation_0_neutral_bottleneck_SFS")
SNP_generation_0_neg_sel_s_0.001_bottleneck_SFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/SNP_generation_0_neg_sel_s_0.001_bottleneck_SFS")
SNP_generation_0_neg_sel_s_0.005_bottleneck_SFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/SNP_generation_0_neg_sel_s_0.005_bottleneck_SFS")
SNP_generation_0_neg_sel_s_0.01_bottleneck_SFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/SNP_generation_0_neg_sel_s_0.01_bottleneck_SFS")
SNP_generation_0_neutral_TE_burst_SFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/SNP_generation_0_neutral_TE_burst_SFS")
SNP_generation_0_neg_sel_s_0.001_TE_burst_SFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/SNP_generation_0_neg_sel_s_0.001_TE_burst_SFS")
SNP_generation_0_neg_sel_s_0.005_TE_burst_SFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/SNP_generation_0_neg_sel_s_0.005_TE_burst_SFS")
SNP_generation_0_neg_sel_s_0.01_TE_burst_SFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/SNP_generation_0_neg_sel_s_0.01_TE_burst_SFS")
SNP_generation_250_neutral_TE_burst_SFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/SNP_generation_250_neutral_TE_burst_SFS")
SNP_generation_250_neg_sel_s_0.001_TE_burst_SFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/SNP_generation_250_neg_sel_s_0.001_TE_burst_SFS")
SNP_generation_250_neg_sel_s_0.005_TE_burst_SFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/SNP_generation_250_neg_sel_s_0.005_TE_burst_SFS")
SNP_generation_250_neg_sel_s_0.01_TE_burst_SFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/SNP_generation_250_neg_sel_s_0.01_TE_burst_SFS")

# AFS (age freqeuncy spectrum; -M AFS)
# each column corresponds to one run 
# first 40 row: proportion of sites at age 250, 500, 750, ..., respectively.
# row 41: total number of observation
TE_generation_0_neutral_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_0_neutral_bottleneck_AFS")
TE_generation_50_neutral_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_50_neutral_bottleneck_AFS")
TE_generation_250_neutral_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_250_neutral_bottleneck_AFS")
TE_generation_300_neutral_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_300_neutral_bottleneck_AFS")
TE_generation_500_neutral_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_500_neutral_bottleneck_AFS")
TE_generation_750_neutral_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_750_neutral_bottleneck_AFS")
TE_generation_1250_neutral_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_1250_neutral_bottleneck_AFS")
TE_generation_2000_neutral_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_2000_neutral_bottleneck_AFS")
TE_generation_5000_neutral_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_5000_neutral_bottleneck_AFS")
TE_generation_0_neg_sel_s_0.0001_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_0_neg_sel_s_0.0001_bottleneck_AFS")
TE_generation_50_neg_sel_s_0.0001_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_50_neg_sel_s_0.0001_bottleneck_AFS")
TE_generation_250_neg_sel_s_0.0001_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_250_neg_sel_s_0.0001_bottleneck_AFS")
TE_generation_300_neg_sel_s_0.0001_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_300_neg_sel_s_0.0001_bottleneck_AFS")
TE_generation_500_neg_sel_s_0.0001_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_500_neg_sel_s_0.0001_bottleneck_AFS")
TE_generation_750_neg_sel_s_0.0001_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_750_neg_sel_s_0.0001_bottleneck_AFS")
TE_generation_1250_neg_sel_s_0.0001_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_1250_neg_sel_s_0.0001_bottleneck_AFS")
TE_generation_2000_neg_sel_s_0.0001_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_2000_neg_sel_s_0.0001_bottleneck_AFS")
TE_generation_5000_neg_sel_s_0.0001_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_5000_neg_sel_s_0.0001_bottleneck_AFS")
TE_generation_0_neg_sel_s_0.001_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_0_neg_sel_s_0.001_bottleneck_AFS")
TE_generation_50_neg_sel_s_0.001_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_50_neg_sel_s_0.001_bottleneck_AFS")
TE_generation_250_neg_sel_s_0.001_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_250_neg_sel_s_0.001_bottleneck_AFS")
TE_generation_300_neg_sel_s_0.001_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_300_neg_sel_s_0.001_bottleneck_AFS")
TE_generation_500_neg_sel_s_0.001_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_500_neg_sel_s_0.001_bottleneck_AFS")
TE_generation_750_neg_sel_s_0.001_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_750_neg_sel_s_0.001_bottleneck_AFS")
TE_generation_1250_neg_sel_s_0.001_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_1250_neg_sel_s_0.001_bottleneck_AFS")
TE_generation_2000_neg_sel_s_0.001_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_2000_neg_sel_s_0.001_bottleneck_AFS")
TE_generation_5000_neg_sel_s_0.001_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_5000_neg_sel_s_0.001_bottleneck_AFS")
TE_generation_0_neg_sel_s_0.005_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_0_neg_sel_s_0.005_bottleneck_AFS")
TE_generation_50_neg_sel_s_0.005_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_50_neg_sel_s_0.005_bottleneck_AFS")
TE_generation_250_neg_sel_s_0.005_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_250_neg_sel_s_0.005_bottleneck_AFS")
TE_generation_300_neg_sel_s_0.005_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_300_neg_sel_s_0.005_bottleneck_AFS")
TE_generation_500_neg_sel_s_0.005_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_500_neg_sel_s_0.005_bottleneck_AFS")
TE_generation_750_neg_sel_s_0.005_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_750_neg_sel_s_0.005_bottleneck_AFS")
TE_generation_1250_neg_sel_s_0.005_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_1250_neg_sel_s_0.005_bottleneck_AFS")
TE_generation_2000_neg_sel_s_0.005_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_2000_neg_sel_s_0.005_bottleneck_AFS")
TE_generation_5000_neg_sel_s_0.005_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_5000_neg_sel_s_0.005_bottleneck_AFS")
TE_generation_0_neg_sel_s_0.01_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_0_neg_sel_s_0.01_bottleneck_AFS")
TE_generation_50_neg_sel_s_0.01_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_50_neg_sel_s_0.01_bottleneck_AFS")
TE_generation_250_neg_sel_s_0.01_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_250_neg_sel_s_0.01_bottleneck_AFS")
TE_generation_300_neg_sel_s_0.01_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_300_neg_sel_s_0.01_bottleneck_AFS")
TE_generation_500_neg_sel_s_0.01_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_500_neg_sel_s_0.01_bottleneck_AFS")
TE_generation_750_neg_sel_s_0.01_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_750_neg_sel_s_0.01_bottleneck_AFS")
TE_generation_1250_neg_sel_s_0.01_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_1250_neg_sel_s_0.01_bottleneck_AFS")
TE_generation_2000_neg_sel_s_0.01_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_2000_neg_sel_s_0.01_bottleneck_AFS")
TE_generation_5000_neg_sel_s_0.01_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_5000_neg_sel_s_0.01_bottleneck_AFS")
TE_generation_0_neutral_TE_burst_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_0_neutral_TE_burst_AFS")
TE_generation_50_neutral_TE_burst_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_50_neutral_TE_burst_AFS")
TE_generation_250_neutral_TE_burst_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_250_neutral_TE_burst_AFS")
TE_generation_300_neutral_TE_burst_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_300_neutral_TE_burst_AFS")
TE_generation_500_neutral_TE_burst_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_500_neutral_TE_burst_AFS")
TE_generation_750_neutral_TE_burst_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_750_neutral_TE_burst_AFS")
TE_generation_1250_neutral_TE_burst_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_1250_neutral_TE_burst_AFS")
TE_generation_2000_neutral_TE_burst_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_2000_neutral_TE_burst_AFS")
TE_generation_5000_neutral_TE_burst_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_5000_neutral_TE_burst_AFS")
TE_generation_0_neg_sel_s_0.0001_TE_burst_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_0_neg_sel_s_0.0001_TE_burst_AFS")
TE_generation_50_neg_sel_s_0.0001_TE_burst_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_50_neg_sel_s_0.0001_TE_burst_AFS")
TE_generation_250_neg_sel_s_0.0001_TE_burst_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_250_neg_sel_s_0.0001_TE_burst_AFS")
TE_generation_300_neg_sel_s_0.0001_TE_burst_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_300_neg_sel_s_0.0001_TE_burst_AFS")
TE_generation_500_neg_sel_s_0.0001_TE_burst_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_500_neg_sel_s_0.0001_TE_burst_AFS")
TE_generation_750_neg_sel_s_0.0001_TE_burst_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_750_neg_sel_s_0.0001_TE_burst_AFS")
TE_generation_1250_neg_sel_s_0.0001_TE_burst_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_1250_neg_sel_s_0.0001_TE_burst_AFS")
TE_generation_2000_neg_sel_s_0.0001_TE_burst_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_2000_neg_sel_s_0.0001_TE_burst_AFS")
TE_generation_5000_neg_sel_s_0.0001_TE_burst_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_5000_neg_sel_s_0.0001_TE_burst_AFS")
TE_generation_0_neg_sel_s_0.001_TE_burst_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_0_neg_sel_s_0.001_TE_burst_AFS")
TE_generation_50_neg_sel_s_0.001_TE_burst_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_50_neg_sel_s_0.001_TE_burst_AFS")
TE_generation_250_neg_sel_s_0.001_TE_burst_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_250_neg_sel_s_0.001_TE_burst_AFS")
TE_generation_300_neg_sel_s_0.001_TE_burst_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_300_neg_sel_s_0.001_TE_burst_AFS")
TE_generation_500_neg_sel_s_0.001_TE_burst_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_500_neg_sel_s_0.001_TE_burst_AFS")
TE_generation_750_neg_sel_s_0.001_TE_burst_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_750_neg_sel_s_0.001_TE_burst_AFS")
TE_generation_1250_neg_sel_s_0.001_TE_burst_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_1250_neg_sel_s_0.001_TE_burst_AFS")
TE_generation_2000_neg_sel_s_0.001_TE_burst_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_2000_neg_sel_s_0.001_TE_burst_AFS")
TE_generation_5000_neg_sel_s_0.001_TE_burst_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_5000_neg_sel_s_0.001_TE_burst_AFS")
TE_generation_0_neg_sel_s_0.005_TE_burst_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_0_neg_sel_s_0.005_TE_burst_AFS")
TE_generation_50_neg_sel_s_0.005_TE_burst_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_50_neg_sel_s_0.005_TE_burst_AFS")
TE_generation_250_neg_sel_s_0.005_TE_burst_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_250_neg_sel_s_0.005_TE_burst_AFS")
TE_generation_300_neg_sel_s_0.005_TE_burst_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_300_neg_sel_s_0.005_TE_burst_AFS")
TE_generation_500_neg_sel_s_0.005_TE_burst_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_500_neg_sel_s_0.005_TE_burst_AFS")
TE_generation_750_neg_sel_s_0.005_TE_burst_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_750_neg_sel_s_0.005_TE_burst_AFS")
TE_generation_1250_neg_sel_s_0.005_TE_burst_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_1250_neg_sel_s_0.005_TE_burst_AFS")
TE_generation_2000_neg_sel_s_0.005_TE_burst_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_2000_neg_sel_s_0.005_TE_burst_AFS")
TE_generation_5000_neg_sel_s_0.005_TE_burst_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_5000_neg_sel_s_0.005_TE_burst_AFS")
TE_generation_0_neg_sel_s_0.01_TE_burst_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_0_neg_sel_s_0.01_TE_burst_AFS")
TE_generation_50_neg_sel_s_0.01_TE_burst_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_50_neg_sel_s_0.01_TE_burst_AFS")
TE_generation_250_neg_sel_s_0.01_TE_burst_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_250_neg_sel_s_0.01_TE_burst_AFS")
TE_generation_300_neg_sel_s_0.01_TE_burst_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_300_neg_sel_s_0.01_TE_burst_AFS")
TE_generation_500_neg_sel_s_0.01_TE_burst_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_500_neg_sel_s_0.01_TE_burst_AFS")
TE_generation_750_neg_sel_s_0.01_TE_burst_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_750_neg_sel_s_0.01_TE_burst_AFS")
TE_generation_1250_neg_sel_s_0.01_TE_burst_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_1250_neg_sel_s_0.01_TE_burst_AFS")
TE_generation_2000_neg_sel_s_0.01_TE_burst_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_2000_neg_sel_s_0.01_TE_burst_AFS")
TE_generation_5000_neg_sel_s_0.01_TE_burst_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_5000_neg_sel_s_0.01_TE_burst_AFS")
TE_generation_0_neutral_TE_burst_and_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_0_neutral_TE_burst_and_bottleneck_AFS")
TE_generation_50_neutral_TE_burst_and_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_50_neutral_TE_burst_and_bottleneck_AFS")
TE_generation_250_neutral_TE_burst_and_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_250_neutral_TE_burst_and_bottleneck_AFS")
TE_generation_300_neutral_TE_burst_and_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_300_neutral_TE_burst_and_bottleneck_AFS")
TE_generation_500_neutral_TE_burst_and_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_500_neutral_TE_burst_and_bottleneck_AFS")
TE_generation_750_neutral_TE_burst_and_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_750_neutral_TE_burst_and_bottleneck_AFS")
TE_generation_1250_neutral_TE_burst_and_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_1250_neutral_TE_burst_and_bottleneck_AFS")
TE_generation_2000_neutral_TE_burst_and_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_2000_neutral_TE_burst_and_bottleneck_AFS")
TE_generation_5000_neutral_TE_burst_and_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_5000_neutral_TE_burst_and_bottleneck_AFS")
TE_generation_0_neg_sel_s_0.0001_TE_burst_and_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_0_neg_sel_s_0.0001_TE_burst_and_bottleneck_AFS")
TE_generation_50_neg_sel_s_0.0001_TE_burst_and_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_50_neg_sel_s_0.0001_TE_burst_and_bottleneck_AFS")
TE_generation_250_neg_sel_s_0.0001_TE_burst_and_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_250_neg_sel_s_0.0001_TE_burst_and_bottleneck_AFS")
TE_generation_300_neg_sel_s_0.0001_TE_burst_and_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_300_neg_sel_s_0.0001_TE_burst_and_bottleneck_AFS")
TE_generation_500_neg_sel_s_0.0001_TE_burst_and_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_500_neg_sel_s_0.0001_TE_burst_and_bottleneck_AFS")
TE_generation_750_neg_sel_s_0.0001_TE_burst_and_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_750_neg_sel_s_0.0001_TE_burst_and_bottleneck_AFS")
TE_generation_1250_neg_sel_s_0.0001_TE_burst_and_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_1250_neg_sel_s_0.0001_TE_burst_and_bottleneck_AFS")
TE_generation_2000_neg_sel_s_0.0001_TE_burst_and_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_2000_neg_sel_s_0.0001_TE_burst_and_bottleneck_AFS")
TE_generation_5000_neg_sel_s_0.0001_TE_burst_and_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_5000_neg_sel_s_0.0001_TE_burst_and_bottleneck_AFS")
TE_generation_0_neg_sel_s_0.001_TE_burst_and_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_0_neg_sel_s_0.001_TE_burst_and_bottleneck_AFS")
TE_generation_50_neg_sel_s_0.001_TE_burst_and_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_50_neg_sel_s_0.001_TE_burst_and_bottleneck_AFS")
TE_generation_250_neg_sel_s_0.001_TE_burst_and_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_250_neg_sel_s_0.001_TE_burst_and_bottleneck_AFS")
TE_generation_300_neg_sel_s_0.001_TE_burst_and_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_300_neg_sel_s_0.001_TE_burst_and_bottleneck_AFS")
TE_generation_500_neg_sel_s_0.001_TE_burst_and_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_500_neg_sel_s_0.001_TE_burst_and_bottleneck_AFS")
TE_generation_750_neg_sel_s_0.001_TE_burst_and_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_750_neg_sel_s_0.001_TE_burst_and_bottleneck_AFS")
TE_generation_1250_neg_sel_s_0.001_TE_burst_and_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_1250_neg_sel_s_0.001_TE_burst_and_bottleneck_AFS")
TE_generation_2000_neg_sel_s_0.001_TE_burst_and_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_2000_neg_sel_s_0.001_TE_burst_and_bottleneck_AFS")
TE_generation_5000_neg_sel_s_0.001_TE_burst_and_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_5000_neg_sel_s_0.001_TE_burst_and_bottleneck_AFS")
TE_generation_0_neg_sel_s_0.005_TE_burst_and_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_0_neg_sel_s_0.005_TE_burst_and_bottleneck_AFS")
TE_generation_50_neg_sel_s_0.005_TE_burst_and_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_50_neg_sel_s_0.005_TE_burst_and_bottleneck_AFS")
TE_generation_250_neg_sel_s_0.005_TE_burst_and_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_250_neg_sel_s_0.005_TE_burst_and_bottleneck_AFS")
TE_generation_300_neg_sel_s_0.005_TE_burst_and_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_300_neg_sel_s_0.005_TE_burst_and_bottleneck_AFS")
TE_generation_500_neg_sel_s_0.005_TE_burst_and_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_500_neg_sel_s_0.005_TE_burst_and_bottleneck_AFS")
TE_generation_750_neg_sel_s_0.005_TE_burst_and_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_750_neg_sel_s_0.005_TE_burst_and_bottleneck_AFS")
TE_generation_1250_neg_sel_s_0.005_TE_burst_and_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_1250_neg_sel_s_0.005_TE_burst_and_bottleneck_AFS")
TE_generation_2000_neg_sel_s_0.005_TE_burst_and_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_2000_neg_sel_s_0.005_TE_burst_and_bottleneck_AFS")
TE_generation_5000_neg_sel_s_0.005_TE_burst_and_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_5000_neg_sel_s_0.005_TE_burst_and_bottleneck_AFS")
TE_generation_0_neg_sel_s_0.01_TE_burst_and_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_0_neg_sel_s_0.01_TE_burst_and_bottleneck_AFS")
TE_generation_50_neg_sel_s_0.01_TE_burst_and_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_50_neg_sel_s_0.01_TE_burst_and_bottleneck_AFS")
TE_generation_250_neg_sel_s_0.01_TE_burst_and_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_250_neg_sel_s_0.01_TE_burst_and_bottleneck_AFS")
TE_generation_300_neg_sel_s_0.01_TE_burst_and_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_300_neg_sel_s_0.01_TE_burst_and_bottleneck_AFS")
TE_generation_500_neg_sel_s_0.01_TE_burst_and_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_500_neg_sel_s_0.01_TE_burst_and_bottleneck_AFS")
TE_generation_750_neg_sel_s_0.01_TE_burst_and_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_750_neg_sel_s_0.01_TE_burst_and_bottleneck_AFS")
TE_generation_1250_neg_sel_s_0.01_TE_burst_and_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_1250_neg_sel_s_0.01_TE_burst_and_bottleneck_AFS")
TE_generation_2000_neg_sel_s_0.01_TE_burst_and_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_2000_neg_sel_s_0.01_TE_burst_and_bottleneck_AFS")
TE_generation_5000_neg_sel_s_0.01_TE_burst_and_bottleneck_AFS <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_5000_neg_sel_s_0.01_TE_burst_and_bottleneck_AFS")

# age distributions (-M AgeDist)
# each column corresponds to one run 
# first 20 row: mean age at the freqeuncies 0.05, 0.10, 0.15, ..., respectively.
# row 21 to 40: number of observation in the bins used to calculate the mean age. (e.g row 21 is the number of observation based on which row 1 was calculated) 
TE_generation_0_neutral_TE_burst_age_dist <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_0_neutral_TE_burst_age_dist")
TE_generation_50_neutral_TE_burst_age_dist <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_50_neutral_TE_burst_age_dist")
TE_generation_250_neutral_TE_burst_age_dist <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_250_neutral_TE_burst_age_dist")
TE_generation_300_neutral_TE_burst_age_dist <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_300_neutral_TE_burst_age_dist")
TE_generation_500_neutral_TE_burst_age_dist <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_500_neutral_TE_burst_age_dist")
TE_generation_750_neutral_TE_burst_age_dist <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_750_neutral_TE_burst_age_dist")
TE_generation_1250_neutral_TE_burst_age_dist <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_1250_neutral_TE_burst_age_dist")
TE_generation_2000_neutral_TE_burst_age_dist <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_2000_neutral_TE_burst_age_dist")
TE_generation_5000_neutral_TE_burst_age_dist <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_5000_neutral_TE_burst_age_dist")
TE_generation_0_neutral_bottleneck_age_dist <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_0_neutral_bottleneck_age_dist")
TE_generation_250_neutral_bottleneck_age_dist <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_250_neutral_bottleneck_age_dist")
TE_generation_500_neutral_bottleneck_age_dist <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_500_neutral_bottleneck_age_dist")
TE_generation_2000_neutral_bottleneck_age_dist <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_2000_neutral_bottleneck_age_dist")
TE_generation_0_neg_sel_s_0.001_TE_burst_age_dist <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_0_neg_sel_s_0.001_TE_burst_age_dist")
TE_generation_0_neg_sel_s_0.005_TE_burst_age_dist <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_0_neg_sel_s_0.005_TE_burst_age_dist")
TE_generation_0_neg_sel_s_0.01_TE_burst_age_dist <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_0_neg_sel_s_0.01_TE_burst_age_dist")
TE_generation_0_neg_sel_s_0.001_bottleneck_age_dist <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_0_neg_sel_s_0.001_bottleneck_age_dist")
TE_generation_0_neg_sel_s_0.005_bottleneck_age_dist <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_0_neg_sel_s_0.005_bottleneck_age_dist")
TE_generation_0_neg_sel_s_0.01_bottleneck_age_dist <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/TE_generation_0_neg_sel_s_0.01_bottleneck_age_dist")
SNP_generation_0_neutral_TE_burst_age_dist <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/SNP_generation_0_neutral_TE_burst_age_dist")
SNP_generation_0_neg_sel_s_0.005_TE_burst_age_dist <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/SNP_generation_0_neg_sel_s_0.005_TE_burst_age_dist")

# delta frequency (-M DeltaFreq)
# each column corresponds to one run 
# rows: delta freqeuncy between age-adjusted site frequency spectrum of TEs and neutral sites for each decile (1st, 2nd, 3rd, ..., respectively).
Generation_0_neutral_TE_burst_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_0_neutral_TE_burst_delta_frequ_0.1_quant_bin")
Generation_50_neutral_TE_burst_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_50_neutral_TE_burst_delta_frequ_0.1_quant_bin")
Generation_250_neutral_TE_burst_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_250_neutral_TE_burst_delta_frequ_0.1_quant_bin")
Generation_300_neutral_TE_burst_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_300_neutral_TE_burst_delta_frequ_0.1_quant_bin")
Generation_500_neutral_TE_burst_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_500_neutral_TE_burst_delta_frequ_0.1_quant_bin")
Generation_750_neutral_TE_burst_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_750_neutral_TE_burst_delta_frequ_0.1_quant_bin")
Generation_1250_neutral_TE_burst_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_1250_neutral_TE_burst_delta_frequ_0.1_quant_bin")
Generation_2000_neutral_TE_burst_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_2000_neutral_TE_burst_delta_frequ_0.1_quant_bin")
Generation_5000_neutral_TE_burst_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_5000_neutral_TE_burst_delta_frequ_0.1_quant_bin")
Generation_0_neg_sel_s_0.001_TE_burst_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_0_neg_sel_s_0.001_TE_burst_delta_frequ_0.1_quant_bin")
Generation_50_neg_sel_s_0.001_TE_burst_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_50_neg_sel_s_0.001_TE_burst_delta_frequ_0.1_quant_bin")
Generation_250_neg_sel_s_0.001_TE_burst_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_250_neg_sel_s_0.001_TE_burst_delta_frequ_0.1_quant_bin")
Generation_300_neg_sel_s_0.001_TE_burst_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_300_neg_sel_s_0.001_TE_burst_delta_frequ_0.1_quant_bin")
Generation_500_neg_sel_s_0.001_TE_burst_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_500_neg_sel_s_0.001_TE_burst_delta_frequ_0.1_quant_bin")
Generation_750_neg_sel_s_0.001_TE_burst_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_750_neg_sel_s_0.001_TE_burst_delta_frequ_0.1_quant_bin")
Generation_1250_neg_sel_s_0.001_TE_burst_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_1250_neg_sel_s_0.001_TE_burst_delta_frequ_0.1_quant_bin")
Generation_2000_neg_sel_s_0.001_TE_burst_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_2000_neg_sel_s_0.001_TE_burst_delta_frequ_0.1_quant_bin")
Generation_5000_neg_sel_s_0.001_TE_burst_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_5000_neg_sel_s_0.001_TE_burst_delta_frequ_0.1_quant_bin")
Generation_0_neg_sel_s_0.005_TE_burst_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_0_neg_sel_s_0.005_TE_burst_delta_frequ_0.1_quant_bin")
Generation_50_neg_sel_s_0.005_TE_burst_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_50_neg_sel_s_0.005_TE_burst_delta_frequ_0.1_quant_bin")
Generation_250_neg_sel_s_0.005_TE_burst_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_250_neg_sel_s_0.005_TE_burst_delta_frequ_0.1_quant_bin")
Generation_300_neg_sel_s_0.005_TE_burst_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_300_neg_sel_s_0.005_TE_burst_delta_frequ_0.1_quant_bin")
Generation_500_neg_sel_s_0.005_TE_burst_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_500_neg_sel_s_0.005_TE_burst_delta_frequ_0.1_quant_bin")
Generation_750_neg_sel_s_0.005_TE_burst_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_750_neg_sel_s_0.005_TE_burst_delta_frequ_0.1_quant_bin")
Generation_1250_neg_sel_s_0.005_TE_burst_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_1250_neg_sel_s_0.005_TE_burst_delta_frequ_0.1_quant_bin")
Generation_2000_neg_sel_s_0.005_TE_burst_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_2000_neg_sel_s_0.005_TE_burst_delta_frequ_0.1_quant_bin")
Generation_5000_neg_sel_s_0.005_TE_burst_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_5000_neg_sel_s_0.005_TE_burst_delta_frequ_0.1_quant_bin")
Generation_0_neg_sel_s_0.01_TE_burst_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_0_neg_sel_s_0.01_TE_burst_delta_frequ_0.1_quant_bin")
Generation_50_neg_sel_s_0.01_TE_burst_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_50_neg_sel_s_0.01_TE_burst_delta_frequ_0.1_quant_bin")
Generation_250_neg_sel_s_0.01_TE_burst_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_250_neg_sel_s_0.01_TE_burst_delta_frequ_0.1_quant_bin")
Generation_300_neg_sel_s_0.01_TE_burst_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_300_neg_sel_s_0.01_TE_burst_delta_frequ_0.1_quant_bin")
Generation_500_neg_sel_s_0.01_TE_burst_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_500_neg_sel_s_0.01_TE_burst_delta_frequ_0.1_quant_bin")
Generation_750_neg_sel_s_0.01_TE_burst_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_750_neg_sel_s_0.01_TE_burst_delta_frequ_0.1_quant_bin")
Generation_1250_neg_sel_s_0.01_TE_burst_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_1250_neg_sel_s_0.01_TE_burst_delta_frequ_0.1_quant_bin")
Generation_2000_neg_sel_s_0.01_TE_burst_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_2000_neg_sel_s_0.01_TE_burst_delta_frequ_0.1_quant_bin")
Generation_5000_neg_sel_s_0.01_TE_burst_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_5000_neg_sel_s_0.01_TE_burst_delta_frequ_0.1_quant_bin")
Generation_0_neg_sel_s_0.0001_TE_burst_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_0_neg_sel_s_0.0001_TE_burst_delta_frequ_0.1_quant_bin")
Generation_50_neg_sel_s_0.0001_TE_burst_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_50_neg_sel_s_0.0001_TE_burst_delta_frequ_0.1_quant_bin")
Generation_250_neg_sel_s_0.0001_TE_burst_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_250_neg_sel_s_0.0001_TE_burst_delta_frequ_0.1_quant_bin")
Generation_300_neg_sel_s_0.0001_TE_burst_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_300_neg_sel_s_0.0001_TE_burst_delta_frequ_0.1_quant_bin")
Generation_500_neg_sel_s_0.0001_TE_burst_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_500_neg_sel_s_0.0001_TE_burst_delta_frequ_0.1_quant_bin")
Generation_750_neg_sel_s_0.0001_TE_burst_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_750_neg_sel_s_0.0001_TE_burst_delta_frequ_0.1_quant_bin")
Generation_1250_neg_sel_s_0.0001_TE_burst_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_1250_neg_sel_s_0.0001_TE_burst_delta_frequ_0.1_quant_bin")
Generation_2000_neg_sel_s_0.0001_TE_burst_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_2000_neg_sel_s_0.0001_TE_burst_delta_frequ_0.1_quant_bin")
Generation_5000_neg_sel_s_0.0001_TE_burst_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_5000_neg_sel_s_0.0001_TE_burst_delta_frequ_0.1_quant_bin")
Generation_0_neutral_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_0_neutral_bottleneck_delta_frequ_0.1_quant_bin")
Generation_50_neutral_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_50_neutral_bottleneck_delta_frequ_0.1_quant_bin")
Generation_250_neutral_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_250_neutral_bottleneck_delta_frequ_0.1_quant_bin")
Generation_300_neutral_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_300_neutral_bottleneck_delta_frequ_0.1_quant_bin")
Generation_500_neutral_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_500_neutral_bottleneck_delta_frequ_0.1_quant_bin")
Generation_750_neutral_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_750_neutral_bottleneck_delta_frequ_0.1_quant_bin")
Generation_1250_neutral_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_1250_neutral_bottleneck_delta_frequ_0.1_quant_bin")
Generation_2000_neutral_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_2000_neutral_bottleneck_delta_frequ_0.1_quant_bin")
Generation_5000_neutral_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_5000_neutral_bottleneck_delta_frequ_0.1_quant_bin")
Generation_0_neg_sel_s_0.001_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_0_neg_sel_s_0.001_bottleneck_delta_frequ_0.1_quant_bin")
Generation_50_neg_sel_s_0.001_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_50_neg_sel_s_0.001_bottleneck_delta_frequ_0.1_quant_bin")
Generation_250_neg_sel_s_0.001_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_250_neg_sel_s_0.001_bottleneck_delta_frequ_0.1_quant_bin")
Generation_300_neg_sel_s_0.001_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_300_neg_sel_s_0.001_bottleneck_delta_frequ_0.1_quant_bin")
Generation_500_neg_sel_s_0.001_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_500_neg_sel_s_0.001_bottleneck_delta_frequ_0.1_quant_bin")
Generation_750_neg_sel_s_0.001_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_750_neg_sel_s_0.001_bottleneck_delta_frequ_0.1_quant_bin")
Generation_1250_neg_sel_s_0.001_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_1250_neg_sel_s_0.001_bottleneck_delta_frequ_0.1_quant_bin")
Generation_2000_neg_sel_s_0.001_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_2000_neg_sel_s_0.001_bottleneck_delta_frequ_0.1_quant_bin")
Generation_5000_neg_sel_s_0.001_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_5000_neg_sel_s_0.001_bottleneck_delta_frequ_0.1_quant_bin")
Generation_0_neg_sel_s_0.005_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_0_neg_sel_s_0.005_bottleneck_delta_frequ_0.1_quant_bin")
Generation_50_neg_sel_s_0.005_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_50_neg_sel_s_0.005_bottleneck_delta_frequ_0.1_quant_bin")
Generation_250_neg_sel_s_0.005_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_250_neg_sel_s_0.005_bottleneck_delta_frequ_0.1_quant_bin")
Generation_300_neg_sel_s_0.005_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_300_neg_sel_s_0.005_bottleneck_delta_frequ_0.1_quant_bin")
Generation_500_neg_sel_s_0.005_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_500_neg_sel_s_0.005_bottleneck_delta_frequ_0.1_quant_bin")
Generation_750_neg_sel_s_0.005_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_750_neg_sel_s_0.005_bottleneck_delta_frequ_0.1_quant_bin")
Generation_1250_neg_sel_s_0.005_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_1250_neg_sel_s_0.005_bottleneck_delta_frequ_0.1_quant_bin")
Generation_2000_neg_sel_s_0.005_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_2000_neg_sel_s_0.005_bottleneck_delta_frequ_0.1_quant_bin")
Generation_5000_neg_sel_s_0.005_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_5000_neg_sel_s_0.005_bottleneck_delta_frequ_0.1_quant_bin")
Generation_0_neg_sel_s_0.01_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_0_neg_sel_s_0.01_bottleneck_delta_frequ_0.1_quant_bin")
Generation_50_neg_sel_s_0.01_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_50_neg_sel_s_0.01_bottleneck_delta_frequ_0.1_quant_bin")
Generation_250_neg_sel_s_0.01_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_250_neg_sel_s_0.01_bottleneck_delta_frequ_0.1_quant_bin")
Generation_300_neg_sel_s_0.01_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_300_neg_sel_s_0.01_bottleneck_delta_frequ_0.1_quant_bin")
Generation_500_neg_sel_s_0.01_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_500_neg_sel_s_0.01_bottleneck_delta_frequ_0.1_quant_bin")
Generation_750_neg_sel_s_0.01_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_750_neg_sel_s_0.01_bottleneck_delta_frequ_0.1_quant_bin")
Generation_1250_neg_sel_s_0.01_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_1250_neg_sel_s_0.01_bottleneck_delta_frequ_0.1_quant_bin")
Generation_2000_neg_sel_s_0.01_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_2000_neg_sel_s_0.01_bottleneck_delta_frequ_0.1_quant_bin")
Generation_5000_neg_sel_s_0.01_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_5000_neg_sel_s_0.01_bottleneck_delta_frequ_0.1_quant_bin")
Generation_0_neg_sel_s_0.0001_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_0_neg_sel_s_0.0001_bottleneck_delta_frequ_0.1_quant_bin")
Generation_50_neg_sel_s_0.0001_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_50_neg_sel_s_0.0001_bottleneck_delta_frequ_0.1_quant_bin")
Generation_250_neg_sel_s_0.0001_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_250_neg_sel_s_0.0001_bottleneck_delta_frequ_0.1_quant_bin")
Generation_300_neg_sel_s_0.0001_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_300_neg_sel_s_0.0001_bottleneck_delta_frequ_0.1_quant_bin")
Generation_500_neg_sel_s_0.0001_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_500_neg_sel_s_0.0001_bottleneck_delta_frequ_0.1_quant_bin")
Generation_750_neg_sel_s_0.0001_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_750_neg_sel_s_0.0001_bottleneck_delta_frequ_0.1_quant_bin")
Generation_1250_neg_sel_s_0.0001_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_1250_neg_sel_s_0.0001_bottleneck_delta_frequ_0.1_quant_bin")
Generation_2000_neg_sel_s_0.0001_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_2000_neg_sel_s_0.0001_bottleneck_delta_frequ_0.1_quant_bin")
Generation_5000_neg_sel_s_0.0001_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_5000_neg_sel_s_0.0001_bottleneck_delta_frequ_0.1_quant_bin")
Generation_0_neutral_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_0_neutral_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin")
Generation_50_neutral_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_50_neutral_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin")
Generation_250_neutral_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_250_neutral_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin")
Generation_300_neutral_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_300_neutral_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin")
Generation_500_neutral_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_500_neutral_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin")
Generation_750_neutral_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_750_neutral_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin")
Generation_1250_neutral_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_1250_neutral_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin")
Generation_2000_neutral_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_2000_neutral_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin")
Generation_5000_neutral_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_5000_neutral_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin")
Generation_0_neg_sel_s_0.0001_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_0_neg_sel_s_0.0001_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin")
Generation_50_neg_sel_s_0.0001_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_50_neg_sel_s_0.0001_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin")
Generation_250_neg_sel_s_0.0001_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_250_neg_sel_s_0.0001_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin")
Generation_300_neg_sel_s_0.0001_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_300_neg_sel_s_0.0001_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin")
Generation_500_neg_sel_s_0.0001_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_500_neg_sel_s_0.0001_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin")
Generation_750_neg_sel_s_0.0001_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_750_neg_sel_s_0.0001_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin")
Generation_1250_neg_sel_s_0.0001_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_1250_neg_sel_s_0.0001_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin")
Generation_2000_neg_sel_s_0.0001_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_2000_neg_sel_s_0.0001_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin")
Generation_5000_neg_sel_s_0.0001_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_5000_neg_sel_s_0.0001_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin")
Generation_0_neg_sel_s_0.005_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_0_neg_sel_s_0.005_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin")
Generation_50_neg_sel_s_0.005_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_50_neg_sel_s_0.005_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin")
Generation_250_neg_sel_s_0.005_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_250_neg_sel_s_0.005_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin")
Generation_300_neg_sel_s_0.005_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_300_neg_sel_s_0.005_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin")
Generation_500_neg_sel_s_0.005_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_500_neg_sel_s_0.005_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin")
Generation_750_neg_sel_s_0.005_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_750_neg_sel_s_0.005_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin")
Generation_1250_neg_sel_s_0.005_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_1250_neg_sel_s_0.005_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin")
Generation_2000_neg_sel_s_0.005_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_2000_neg_sel_s_0.005_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin")
Generation_5000_neg_sel_s_0.005_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_5000_neg_sel_s_0.005_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin")
Generation_0_neg_sel_s_0.001_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_0_neg_sel_s_0.001_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin")
Generation_50_neg_sel_s_0.001_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_50_neg_sel_s_0.001_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin")
Generation_250_neg_sel_s_0.001_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_250_neg_sel_s_0.001_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin")
Generation_300_neg_sel_s_0.001_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_300_neg_sel_s_0.001_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin")
Generation_500_neg_sel_s_0.001_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_500_neg_sel_s_0.001_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin")
Generation_750_neg_sel_s_0.001_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_750_neg_sel_s_0.001_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin")
Generation_1250_neg_sel_s_0.001_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_1250_neg_sel_s_0.001_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin")
Generation_2000_neg_sel_s_0.001_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_2000_neg_sel_s_0.001_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin")
Generation_5000_neg_sel_s_0.001_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_5000_neg_sel_s_0.001_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin")
Generation_0_neg_sel_s_0.01_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_0_neg_sel_s_0.01_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin")
Generation_50_neg_sel_s_0.01_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_50_neg_sel_s_0.01_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin")
Generation_250_neg_sel_s_0.01_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_250_neg_sel_s_0.01_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin")
Generation_300_neg_sel_s_0.01_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_300_neg_sel_s_0.01_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin")
Generation_500_neg_sel_s_0.01_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_500_neg_sel_s_0.01_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin")
Generation_750_neg_sel_s_0.01_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_750_neg_sel_s_0.01_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin")
Generation_1250_neg_sel_s_0.01_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_1250_neg_sel_s_0.01_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin")
Generation_2000_neg_sel_s_0.01_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_2000_neg_sel_s_0.01_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin")
Generation_5000_neg_sel_s_0.01_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin <- read.table("/Users/roberthorvath/Desktop/Projects/2) TE SLiM simulations/8) SLiM Simulation/SLiM_output_summary_secound_run/Generation_5000_neg_sel_s_0.01_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin")


#######################
#### Analyse the data:

## The Impact of Selection, Demography and Transposition Rate on the Age Distribution
# Calculate mean age at a specific freqeuncy
G0_neutral_TE_burst_age_dist <- function_get_mean_age(TE_generation_0_neutral_TE_burst_age_dist)
G50_neutral_TE_burst_age_dist <- function_get_mean_age(TE_generation_50_neutral_TE_burst_age_dist)
G250_neutral_TE_burst_age_dist <- function_get_mean_age(TE_generation_250_neutral_TE_burst_age_dist)
G300_neutral_TE_burst_age_dist <- function_get_mean_age(TE_generation_300_neutral_TE_burst_age_dist)
G500_neutral_TE_burst_age_dist <- function_get_mean_age(TE_generation_500_neutral_TE_burst_age_dist)
G750_neutral_TE_burst_age_dist <- function_get_mean_age(TE_generation_750_neutral_TE_burst_age_dist)
G1250_neutral_TE_burst_age_dist <- function_get_mean_age(TE_generation_1250_neutral_TE_burst_age_dist)
G2000_neutral_TE_burst_age_dist <- function_get_mean_age(TE_generation_2000_neutral_TE_burst_age_dist)
G5000_neutral_TE_burst_age_dist <- function_get_mean_age(TE_generation_5000_neutral_TE_burst_age_dist)

G0_neutral_bottleneck_age_dist <- function_get_mean_age(TE_generation_0_neutral_bottleneck_age_dist)
G250_neutral_bottleneck_age_dist <- function_get_mean_age(TE_generation_250_neutral_bottleneck_age_dist)
G500_neutral_bottleneck_age_dist <- function_get_mean_age(TE_generation_500_neutral_bottleneck_age_dist)
G2000_neutral_bottleneck_age_dist <- function_get_mean_age(TE_generation_2000_neutral_bottleneck_age_dist)

# plot the allele age-frequency relationship
pdf("/Users/roberthorvath/Desktop/Figure 2.pdf" )
par(xpd=TRUE, mfrow=c(2,2))

function_plot_mean_age(TE_generation_0_neutral_TE_burst_age_dist, border = 2, ylim = c(0, 3000), ylab = "Age (generations)", xlab = "Allele frequency", main = "")
function_plot_mean_age(SNP_generation_0_neutral_TE_burst_age_dist, border = 8, add = TRUE)
axis(side = 1, at = seq(0.5, 20.5, 1), labels = seq(0, 1, 0.05), las = 2)
par(xpd=TRUE)
legend(-0.5, 4000, legend = c("Neutral SNP age"), bty = "n", pch = 1, col = 8)
legend(11, 4000, legend = c("Observed TE age", "Predicted allele age"), bty = "n", pch = c(1, 4), col = c(2,1))
par(xpd=FALSE)
for (i in 1:20) {
  points(i, expected_allele_age(i/20 - 0.025, 500), col = 1, pch = 4, cex = 0.75)
}
mtext("A", side = 3, line = 2, at = -2.9, cex = 1.2)

function_plot_mean_age(TE_generation_0_neg_sel_s_0.005_TE_burst_age_dist, border = 2, ylim = c(0, 3000), ylab = "Age (generations)", xlab = "Allele frequency", main = "")
axis(side = 1, at = seq(0.5, 20.5, 1), labels = seq(0, 1, 0.05), las = 2)
par(xpd=TRUE)
legend(0, 4000, legend = c("Observed TE age", "Predicted allele age"), bty = "n", pch = c(1,4), col = c(2,1))
par(xpd=FALSE)
for (i in 1:20) {
  points(i, my_maruyama_age((4*500*-0.005), 500, (i/20 - 0.025) ), col = 1, pch = 4, cex = 0.75)
}
mtext("B", side = 3, line = 2, at = -2.9, cex = 1.2)

plot(1:20, ylim = c(0, 2500), type="n", ylab = "Age (generations)", xlab = "Allele frequency", xaxt = "n")
lines(lowess(G0_neutral_bottleneck_age_dist, f = 0.2), col  = 1)
points(G0_neutral_bottleneck_age_dist, col  = 1, pch =16, cex = 0.5)
lines(lowess(G250_neutral_bottleneck_age_dist, f = 0.2), col  = "tan4")
points(G250_neutral_bottleneck_age_dist, col  = "tan4", pch =16, cex = 0.5)
lines(lowess(G500_neutral_bottleneck_age_dist, f = 0.2), col  = "tan3")
points(G500_neutral_bottleneck_age_dist, col  = "tan3", pch =16, cex = 0.5)
lines(lowess(G2000_neutral_bottleneck_age_dist, f = 0.2), col  = "tan1")
points(G2000_neutral_bottleneck_age_dist, col  = "tan1", pch =16, cex = 0.5)
axis(side = 1, at = seq(0.5, 20.5, 1), labels = seq(0, 1, 0.05), las = 2)
par(xpd=TRUE)
legend(0, 3500, legend = c("G 0", "G 250"), bty = "n", pch = 15, col = c(1,"tan4"))
legend(5.5, 3500, legend = c("G 500", "G 2000"), bty = "n", pch = 15, col = c("tan3","tan1"))
par(xpd=FALSE)
mtext("C", side = 3, line = 2, at = -2.9, cex = 1.2)

plot(1:20, ylim = c(0, 2500), type="n", ylab = "Age (generations)", xlab = "Allele frequency", xaxt = "n")
lines(lowess(G0_neutral_TE_burst_age_dist, f = 0.1), col  = 1)
points(G0_neutral_TE_burst_age_dist, col  = 1, pch =16, cex = 0.5)
lines(lowess(G250_neutral_TE_burst_age_dist, f = 0.1), col  = "maroon4")
points(G250_neutral_TE_burst_age_dist, col  = "maroon4", pch =16, cex = 0.5)
lines(lowess(G500_neutral_TE_burst_age_dist, f = 0.1), col  = "maroon3")
points(G500_neutral_TE_burst_age_dist, col  = "maroon3", pch =16, cex = 0.5)
lines(lowess(G1250_neutral_TE_burst_age_dist, f = 0.1), col  = "maroon1")
points(G1250_neutral_TE_burst_age_dist, col  = "maroon1", pch =16, cex = 0.5)
axis(side = 1, at = seq(0.5, 20.5, 1), labels = seq(0, 1, 0.05), las = 2)
par(xpd=TRUE)
legend(0, 3500, legend = c("G 0", "G 250"), bty = "n", pch = 15, col = c(1,"maroon4"))
legend(5.5, 3500, legend = c("G 500", "G 1250"), bty = "n", pch = 15, col = c("maroon3","maroon1"))
par(xpd=FALSE)
mtext("D", side = 3, line = 2, at = -2.9, cex = 1.2)

dev.off()

# plot the expected and observed mean age distribution of TEs at a specific frequency
pdf("/Users/roberthorvath/Desktop/Figure Sup 1.pdf", height = 10)
par(xpd=TRUE, mfrow=c(4,2), mar = c(4.1, 4.5, 4, 0.5))

function_plot_mean_age(TE_generation_0_neutral_bottleneck_age_dist, border = 2, ylim = c(0, 3000), ylab = "Age", xlab = "Allele frequency", main = "", cex.lab = 1.2)
axis(side = 1, at = seq(0.5, 20.5, 1), labels = seq(0, 1, 0.05), las = 2)
par(xpd=TRUE)
legend(0, 4000, legend = c("Observed TE age", "Predicted allele age"), bty = "n", pch = c(1,4), col = c(2,1))
par(xpd=FALSE)
for (i in 1:20) {
  points(i, expected_allele_age(i/20 - 0.025, 500), col = 1, pch = 4, cex = 0.75)
}
mtext(expression("A"[1]), side = 3, line = 2, at = -2.9, cex = 1)

function_plot_mean_age(TE_generation_0_neutral_TE_burst_age_dist, border = 2, ylim = c(0, 3000), ylab = "Age", xlab = "Allele frequency", main = "", cex.lab = 1.2)
axis(side = 1, at = seq(0.5, 20.5, 1), labels = seq(0, 1, 0.05), las = 2)
par(xpd=TRUE)
legend(0, 4000, legend = c("Observed TE age", "Predicted allele age"), bty = "n", pch = c(1,4), col = c(2,1))
par(xpd=FALSE)
for (i in 1:20) {
  points(i, expected_allele_age(i/20 - 0.025, 500), col = 1, pch = 4, cex = 0.75)
}
mtext(expression("B"[1]), side = 3, line = 2, at = -2.9, cex = 1)

function_plot_mean_age(TE_generation_0_neg_sel_s_0.001_bottleneck_age_dist, border = 2, ylim = c(0, 3000), ylab = "Age", xlab = "Allele frequency", main = "", cex.lab = 1.2)
axis(side = 1, at = seq(0.5, 20.5, 1), labels = seq(0, 1, 0.05), las = 2)
par(xpd=TRUE)
legend(0, 4000, legend = c("Observed TE age", "Predicted allele age"), bty = "n", pch = c(1,4), col = c(2,1))
par(xpd=FALSE)
for (i in 1:20) {
  points(i, my_maruyama_age((4*500*-0.001), 500, (i/20 - 0.025) ), col = 1, pch = 4, cex = 0.75)
}
mtext(expression("A"[2]), side = 3, line = 2, at = -2.9, cex = 1)

function_plot_mean_age(TE_generation_0_neg_sel_s_0.001_TE_burst_age_dist, border = 2, ylim = c(0, 3000), ylab = "Age", xlab = "Allele frequency", main = "", cex.lab = 1.2)
axis(side = 1, at = seq(0.5, 20.5, 1), labels = seq(0, 1, 0.05), las = 2)
par(xpd=TRUE)
legend(0, 4000, legend = c("Observed TE age", "Predicted allele age"), bty = "n", pch = c(1,4), col = c(2,1))
par(xpd=FALSE)
for (i in 1:20) {
  points(i, my_maruyama_age((4*500*-0.001), 500, (i/20 - 0.025) ), col = 1, pch = 4, cex = 0.75)
}
mtext(expression("B"[2]), side = 3, line = 2, at = -2.9, cex = 1)

function_plot_mean_age(TE_generation_0_neg_sel_s_0.005_bottleneck_age_dist, border = 2, ylim = c(0, 3000), ylab = "Age", xlab = "Allele frequency", main = "", cex.lab = 1.2)
axis(side = 1, at = seq(0.5, 20.5, 1), labels = seq(0, 1, 0.05), las = 2)
par(xpd=TRUE)
legend(0, 4000, legend = c("Observed TE age", "Predicted allele age"), bty = "n", pch = c(1,4), col = c(2,1))
par(xpd=FALSE)
for (i in 1:20) {
  points(i, my_maruyama_age((4*500*-0.005), 500, (i/20 - 0.025) ), col = 1, pch = 4, cex = 0.75)
}
mtext(expression("A"[3]), side = 3, line = 2, at = -2.9, cex = 1)

function_plot_mean_age(TE_generation_0_neg_sel_s_0.005_TE_burst_age_dist, border = 2, ylim = c(0, 3000), ylab = "Age", xlab = "Allele frequency", main = "", cex.lab = 1.2)
axis(side = 1, at = seq(0.5, 20.5, 1), labels = seq(0, 1, 0.05), las = 2)
par(xpd=TRUE)
legend(0, 4000, legend = c("Observed TE age", "Predicted allele age"), bty = "n", pch = c(1,4), col = c(2,1))
par(xpd=FALSE)
for (i in 1:20) {
  points(i, my_maruyama_age((4*500*-0.005), 500, (i/20 - 0.025) ), col = 1, pch = 4, cex = 0.75)
}
mtext(expression("B"[3]), side = 3, line = 2, at = -2.9, cex = 1)

function_plot_mean_age(TE_generation_0_neg_sel_s_0.01_bottleneck_age_dist, border = 2, ylim = c(0, 3000), ylab = "Age", xlab = "Allele frequency", main = "", cex.lab = 1.2)
axis(side = 1, at = seq(0.5, 20.5, 1), labels = seq(0, 1, 0.05), las = 2)
par(xpd=TRUE)
legend(0, 4000, legend = c("Observed TE age", "Predicted allele age"), bty = "n", pch = c(1,4), col = c(2,1))
par(xpd=FALSE)
for (i in 1:20) {
  points(i, my_maruyama_age((4*500*-0.01), 500, (i/20 - 0.025) ), col = 1, pch = 4, cex = 0.75)
}
mtext(expression("A"[4]), side = 3, line = 2, at = -2.9, cex = 1)

function_plot_mean_age(TE_generation_0_neg_sel_s_0.01_TE_burst_age_dist, border = 2, ylim = c(0, 3000), ylab = "Age", xlab = "Allele frequency", main = "", cex.lab = 1.2)
axis(side = 1, at = seq(0.5, 20.5, 1), labels = seq(0, 1, 0.05), las = 2)
par(xpd=TRUE)
legend(0, 4000, legend = c("Observed TE age", "Predicted allele age"), bty = "n", pch = c(1,4), col = c(2,1))
par(xpd=FALSE)
for (i in 1:20) {
  points(i, my_maruyama_age((4*500*-0.01), 500, (i/20 - 0.025) ), col = 1, pch = 4, cex = 0.75)
}
mtext(expression("B"[4]), side = 3, line = 2, at = -2.9, cex = 1)

dev.off()

# plot the difference in age between TEs and SNPs
pdf("/Users/roberthorvath/Desktop/Figure Sup 2.pdf", height = 5)
par(xpd=TRUE)

boxplot(as.numeric(TE_generation_0_neg_sel_s_0.005_TE_burst_age_dist[1,]) - as.numeric(SNP_generation_0_neg_sel_s_0.005_TE_burst_age_dist[1,]), 
        as.numeric(TE_generation_0_neg_sel_s_0.005_TE_burst_age_dist[2,]) - as.numeric(SNP_generation_0_neg_sel_s_0.005_TE_burst_age_dist[2,]),
        as.numeric(TE_generation_0_neg_sel_s_0.005_TE_burst_age_dist[3,]) - as.numeric(SNP_generation_0_neg_sel_s_0.005_TE_burst_age_dist[3,]),
        as.numeric(TE_generation_0_neg_sel_s_0.005_TE_burst_age_dist[4,]) - as.numeric(SNP_generation_0_neg_sel_s_0.005_TE_burst_age_dist[4,]),
        as.numeric(TE_generation_0_neg_sel_s_0.005_TE_burst_age_dist[5,]) - as.numeric(SNP_generation_0_neg_sel_s_0.005_TE_burst_age_dist[5,]),
        as.numeric(TE_generation_0_neg_sel_s_0.005_TE_burst_age_dist[6,]) - as.numeric(SNP_generation_0_neg_sel_s_0.005_TE_burst_age_dist[6,]),
        as.numeric(TE_generation_0_neg_sel_s_0.005_TE_burst_age_dist[7,]) - as.numeric(SNP_generation_0_neg_sel_s_0.005_TE_burst_age_dist[7,]),
        as.numeric(TE_generation_0_neg_sel_s_0.005_TE_burst_age_dist[8,]) - as.numeric(SNP_generation_0_neg_sel_s_0.005_TE_burst_age_dist[8,]),
        as.numeric(TE_generation_0_neg_sel_s_0.005_TE_burst_age_dist[9,]) - as.numeric(SNP_generation_0_neg_sel_s_0.005_TE_burst_age_dist[9,]),
        as.numeric(TE_generation_0_neg_sel_s_0.005_TE_burst_age_dist[10,]) - as.numeric(SNP_generation_0_neg_sel_s_0.005_TE_burst_age_dist[10,]),
        as.numeric(TE_generation_0_neg_sel_s_0.005_TE_burst_age_dist[11,]) - as.numeric(SNP_generation_0_neg_sel_s_0.005_TE_burst_age_dist[11,]),
        as.numeric(TE_generation_0_neg_sel_s_0.005_TE_burst_age_dist[12,]) - as.numeric(SNP_generation_0_neg_sel_s_0.005_TE_burst_age_dist[12,]),
        as.numeric(TE_generation_0_neg_sel_s_0.005_TE_burst_age_dist[13,]) - as.numeric(SNP_generation_0_neg_sel_s_0.005_TE_burst_age_dist[13,]),
        as.numeric(TE_generation_0_neg_sel_s_0.005_TE_burst_age_dist[14,]) - as.numeric(SNP_generation_0_neg_sel_s_0.005_TE_burst_age_dist[14,]),
        as.numeric(TE_generation_0_neg_sel_s_0.005_TE_burst_age_dist[15,]) - as.numeric(SNP_generation_0_neg_sel_s_0.005_TE_burst_age_dist[15,]),
        as.numeric(TE_generation_0_neg_sel_s_0.005_TE_burst_age_dist[16,]) - as.numeric(SNP_generation_0_neg_sel_s_0.005_TE_burst_age_dist[16,]),
        as.numeric(TE_generation_0_neg_sel_s_0.005_TE_burst_age_dist[17,]) - as.numeric(SNP_generation_0_neg_sel_s_0.005_TE_burst_age_dist[17,]),
        as.numeric(TE_generation_0_neg_sel_s_0.005_TE_burst_age_dist[18,]) - as.numeric(SNP_generation_0_neg_sel_s_0.005_TE_burst_age_dist[18,]),
        as.numeric(TE_generation_0_neg_sel_s_0.005_TE_burst_age_dist[19,]) - as.numeric(SNP_generation_0_neg_sel_s_0.005_TE_burst_age_dist[19,]),
        as.numeric(TE_generation_0_neg_sel_s_0.005_TE_burst_age_dist[20,]) - as.numeric(SNP_generation_0_neg_sel_s_0.005_TE_burst_age_dist[20,]), xaxt="n", border = 2, xlab = "Frequency bin", ylab = expression(paste(Delta, " age")))
par(xpd=FALSE)
abline(h = 0, lty = 2)
par(xpd=TRUE)
axis(side = 1, at = seq(0.5, 20.5, 1), labels = seq(0, 1, 0.05), las = 2)
function_add_points_expected_delta_maruyama(S = (4*500*-0.005), col = 1, pch = 4)
legend(-0.5, 1500, legend = c("Observed", "Expected"), bty = "n", pch = c(21, 4), col = c(2,1))

dev.off()

# plot the mean age distribution and proportion of TEs at a given frequency in the TE burst model
pdf("/Users/roberthorvath/Desktop/Figure Sup 3.pdf", height = 6)
par(xpd=TRUE, mfrow=c(2,1), mar = c(2, 4.6, 4.1, 2.1))

plot(1:20, ylim = c(0, 2500), type="n", ylab = "Age", xlab = "Allele frequency", xaxt = "n")
lines(lowess(G0_neutral_TE_burst_age_dist, f = 0.1), col  = 1)
points(G0_neutral_TE_burst_age_dist, col  = 1, pch =16, cex = 0.5)
lines(lowess(G50_neutral_TE_burst_age_dist, f = 0.1), col  = rainbow(8)[1])
points(G50_neutral_TE_burst_age_dist, col  = rainbow(8)[1], pch =16, cex = 0.5)
lines(lowess(G250_neutral_TE_burst_age_dist, f = 0.1), col  = rainbow(8)[2])
points(G250_neutral_TE_burst_age_dist, col  = rainbow(8)[2], pch =16, cex = 0.5)
lines(lowess(G300_neutral_TE_burst_age_dist, f = 0.1), col  = rainbow(8)[3])
points(G300_neutral_TE_burst_age_dist, col  = rainbow(8)[3], pch =16, cex = 0.5)
lines(lowess(G500_neutral_TE_burst_age_dist, f = 0.1), col  = rainbow(8)[4])
points(G500_neutral_TE_burst_age_dist, col  = rainbow(8)[4], pch =16, cex = 0.5)
lines(lowess(G750_neutral_TE_burst_age_dist, f = 0.1), col  = rainbow(8)[5])
points(G750_neutral_TE_burst_age_dist, col  = rainbow(8)[5], pch =16, cex = 0.5)
lines(lowess(G1250_neutral_TE_burst_age_dist, f = 0.1), col  = rainbow(8)[6])
points(G1250_neutral_TE_burst_age_dist, col  = rainbow(8)[6], pch =16, cex = 0.5)
lines(lowess(G2000_neutral_TE_burst_age_dist, f = 0.1), col  = rainbow(8)[7])
points(G2000_neutral_TE_burst_age_dist, col  = rainbow(8)[7], pch =16, cex = 0.5)
lines(lowess(G5000_neutral_TE_burst_age_dist, f = 0.1), col  = rainbow(8)[8])
points(G5000_neutral_TE_burst_age_dist, col  = rainbow(8)[8], pch =16, cex = 0.5)
axis(side = 1, at = seq(0.5, 20.5, 1), labels = seq(0, 1, 0.05), las = 2)
par(xpd=TRUE)
legend(-3, 3400, legend = "G 0", bty = "n", pch = 15, col = 1)
legend(-1, 3400, legend = "G 50", bty = "n", pch = 15, col = rainbow(8)[1])
legend(1.25, 3400, legend = "G 250", bty = "n", pch = 15, col = rainbow(8)[2])
legend(4, 3400, legend = "G 300", bty = "n", pch = 15, col = rainbow(8)[3])
legend(6.75, 3400, legend = "G 500", bty = "n", pch = 15, col = rainbow(8)[4])
legend(9.5, 3400, legend = "G 750", bty = "n", pch = 15, col = rainbow(8)[5])
legend(12.25, 3400, legend = "G 1250", bty = "n", pch = 15, col = rainbow(8)[6])
legend(15.25, 3400, legend = "G 2000", bty = "n", pch = 15, col = rainbow(8)[7])
legend(18.25, 3400, legend = "G 5000", bty = "n", pch = 15, col = rainbow(8)[8])
par(xpd=FALSE)

par(mar = c(5.1, 4.1, 1.5, 1.4))
barplot(as.vector(rbind(
  TE_generation_0_neutral_TE_burst_SFS[1:20,1],
  TE_generation_50_neutral_TE_burst_SFS[1:20,1],
  TE_generation_250_neutral_TE_burst_SFS[1:20,1],
  TE_generation_300_neutral_TE_burst_SFS[1:20,1],
  TE_generation_500_neutral_TE_burst_SFS[1:20,1], 
  TE_generation_750_neutral_TE_burst_SFS[1:20,1], 
  TE_generation_1250_neutral_TE_burst_SFS[1:20,1], 
  TE_generation_2000_neutral_TE_burst_SFS[1:20,1], 
  TE_generation_5000_neutral_TE_burst_SFS[1:20,1])), 
  col = c(1, rainbow(8)), space = c(1, 0, 0, 0, 0, 0, 0, 0, 0), 
  ylim = c(0, 1), xlab = "Allele frequency", ylab = "Proportion"
)
axis(1, seq(0.5, 200.5, 10), seq(0, 1, 0.05))
for (i in 1:20) {
  arrows(i + 0.5 + (i -1)*9, as.numeric(quantile(TE_generation_0_neutral_TE_burst_SFS[i,], 1)), i + 0.5 + (i -1)*9, as.numeric(quantile(TE_generation_0_neutral_TE_burst_SFS[i,], 0)), angle = 90, code = 3, length = 0.01)
  arrows(i + 1.5 + (i -1)*9, as.numeric(quantile(TE_generation_50_neutral_TE_burst_SFS[i,], 1)), i + 1.5 + (i -1)*9, as.numeric(quantile(TE_generation_50_neutral_TE_burst_SFS[i,], 0)), angle = 90, code = 3, length = 0.01)
  arrows(i + 2.5 + (i -1)*9, as.numeric(quantile(TE_generation_250_neutral_TE_burst_SFS[i,], 1)), i + 2.5 + (i -1)*9, as.numeric(quantile(TE_generation_250_neutral_TE_burst_SFS[i,], 0)), angle = 90, code = 3, length = 0.01)
  arrows(i + 3.5 + (i -1)*9, as.numeric(quantile(TE_generation_300_neutral_TE_burst_SFS[i,], 1)), i + 3.5 + (i -1)*9, as.numeric(quantile(TE_generation_300_neutral_TE_burst_SFS[i,], 0)), angle = 90, code = 3, length = 0.01)
  arrows(i + 4.5 + (i -1)*9, as.numeric(quantile(TE_generation_500_neutral_TE_burst_SFS[i,], 1)), i + 4.5 + (i -1)*9, as.numeric(quantile(TE_generation_500_neutral_TE_burst_SFS[i,], 0)), angle = 90, code = 3, length = 0.01)
  arrows(i + 5.5 + (i -1)*9, as.numeric(quantile(TE_generation_750_neutral_TE_burst_SFS[i,], 1)), i + 5.5 + (i -1)*9, as.numeric(quantile(TE_generation_750_neutral_TE_burst_SFS[i,], 0)), angle = 90, code = 3, length = 0.01)
  arrows(i + 6.5 + (i -1)*9, as.numeric(quantile(TE_generation_1250_neutral_TE_burst_SFS[i,], 1)), i + 6.5 + (i -1)*9, as.numeric(quantile(TE_generation_1250_neutral_TE_burst_SFS[i,], 0)), angle = 90, code = 3, length = 0.01)
  arrows(i + 7.5 + (i -1)*9, as.numeric(quantile(TE_generation_2000_neutral_TE_burst_SFS[i,], 1)), i + 7.5 + (i -1)*9, as.numeric(quantile(TE_generation_2000_neutral_TE_burst_SFS[i,], 0)), angle = 90, code = 3, length = 0.01)
  arrows(i + 8.5 + (i -1)*9, as.numeric(quantile(TE_generation_5000_neutral_TE_burst_SFS[i,], 1)), i + 8.5 + (i -1)*9, as.numeric(quantile(TE_generation_5000_neutral_TE_burst_SFS[i,], 0)), angle = 90, code = 3, length = 0.01)
}

dev.off()


# Conditioning on Age to Identify Allele Frequency Shifts Caused by Selective
# plot delta freqeuncy
pdf("/Users/roberthorvath/Desktop/Figure 3.pdf", height = 7)
layout(mat = matrix(c(seq(1,13,4), seq(2,14,4), seq(3,15,4), seq(4,16,4) ), nrow = 4, ncol = 4), heights = c(1, 1, 1, 1.2), widths = c(1, 0.9, 0.9, 0.9))
par(xpd=TRUE, mar = c(0.5, 3.5, 2.5, 0.5), mgp = c(2, 1, 0))

function_delta_frequ_plot(Generation_0_neutral_bottleneck_delta_frequ_0.1_quant_bin, ylim = c(-0.15, 0.15), my_xlab = "")
par(xpd=FALSE)
abline(h = 0, lty = 2)
par(xpd=TRUE)
mtext(expression("A"[1]), side = 3, line = 0.7, at = -2.5, cex = 0.9)
par(xpd=TRUE, mar = c(0.5, 2.4, 2.5, 0.5), mgp = c(2, 1, 0))

function_delta_frequ_plot(Generation_0_neg_sel_s_0.005_bottleneck_delta_frequ_0.1_quant_bin, ylim = c(-0.2, 0.2), my_xlab = "", my_ylab = "")
par(xpd=FALSE)
abline(h = 0, lty = 2)
par(xpd=TRUE)
mtext(expression("B"[1]), side = 3, line = 0.7, at = -2, cex = 0.9)

function_delta_frequ_plot(Generation_0_neutral_TE_burst_delta_frequ_0.1_quant_bin, ylim = c(-0.04, 0.04), my_xlab = "", my_ylab = "")
par(xpd=FALSE)
abline(h = 0, lty = 2)
par(xpd=TRUE)
mtext(expression("C"[1]), side = 3, line = 0.7, at = -2, cex = 0.9)

function_delta_frequ_plot(Generation_0_neg_sel_s_0.005_TE_burst_delta_frequ_0.1_quant_bin, ylim = c(-0.16, 0.05), my_xlab = "", my_ylab = "")
par(xpd=FALSE)
abline(h = 0, lty = 2)
par(xpd=TRUE)
mtext(expression("D"[1]), side = 3, line = 0.7, at = -2, cex = 0.9)
par(xpd=TRUE, mar = c(0.5, 3.5, 2.5, 0.5), mgp = c(2, 1, 0))

function_delta_frequ_plot(Generation_50_neutral_bottleneck_delta_frequ_0.1_quant_bin, ylim = c(-0.15, 0.15), my_xlab = "")
par(xpd=FALSE)
abline(h = 0, lty = 2)
par(xpd=TRUE)
mtext(expression("A"[2]), side = 3, line = 0.7, at = -2.5, cex = 0.9)
par(xpd=TRUE, mar = c(0.5, 2.4, 2.5, 0.5), mgp = c(2, 1, 0))

function_delta_frequ_plot(Generation_50_neg_sel_s_0.005_bottleneck_delta_frequ_0.1_quant_bin, ylim = c(-0.2, 0.2), my_xlab = "", my_ylab = "")
par(xpd=FALSE)
abline(h = 0, lty = 2)
par(xpd=TRUE)
mtext(expression("B"[2]), side = 3, line = 0.7, at = -2, cex = 0.9)

function_delta_frequ_plot(Generation_50_neutral_TE_burst_delta_frequ_0.1_quant_bin, ylim = c(-0.04, 0.04), my_xlab = "", my_ylab = "")
par(xpd=FALSE)
abline(h = 0, lty = 2)
par(xpd=TRUE)
mtext(expression("C"[2]), side = 3, line = 0.7, at = -2, cex = 0.9)

function_delta_frequ_plot(Generation_50_neg_sel_s_0.005_TE_burst_delta_frequ_0.1_quant_bin, ylim = c(-0.16, 0.05), my_xlab = "", my_ylab = "")
par(xpd=FALSE)
abline(h = 0, lty = 2)
par(xpd=TRUE)
mtext(expression("D"[2]), side = 3, line = 0.7, at = -2, cex = 0.9)
par(xpd=TRUE, mar = c(0.5, 3.5, 2.5, 0.5), mgp = c(2, 1, 0))

function_delta_frequ_plot(Generation_250_neutral_bottleneck_delta_frequ_0.1_quant_bin, ylim = c(-0.15, 0.15), my_xlab = "")
par(xpd=FALSE)
abline(h = 0, lty = 2)
par(xpd=TRUE)
mtext(expression("A"[3]), side = 3, line = 0.7, at = -2.5, cex = 0.9)
par(xpd=TRUE, mar = c(0.5, 2.4, 2.5, 0.5), mgp = c(2, 1, 0))

function_delta_frequ_plot(Generation_250_neg_sel_s_0.005_bottleneck_delta_frequ_0.1_quant_bin, ylim = c(-0.2, 0.2), my_xlab = "" , my_ylab = "")
par(xpd=FALSE)
abline(h = 0, lty = 2)
par(xpd=TRUE)
mtext(expression("B"[3]), side = 3, line = 0.7, at = -2, cex = 0.9)

function_delta_frequ_plot(Generation_250_neutral_TE_burst_delta_frequ_0.1_quant_bin, ylim = c(-0.04, 0.04), my_xlab = "", my_ylab = "")
par(xpd=FALSE)
abline(h = 0, lty = 2)
par(xpd=TRUE)
mtext(expression("C"[3]), side = 3, line = 0.7, at = -2, cex = 0.9)

function_delta_frequ_plot(Generation_250_neg_sel_s_0.005_TE_burst_delta_frequ_0.1_quant_bin, ylim = c(-0.16, 0.05), my_xlab = "", my_ylab = "")
par(xpd=FALSE)
abline(h = 0, lty = 2)
par(xpd=TRUE)
mtext(expression("D"[3]), side = 3, line = 0.7, at = -2, cex = 0.9)
par(mar = c(3, 3.5, 2.5, 0.5))

function_delta_frequ_plot(Generation_1250_neutral_bottleneck_delta_frequ_0.1_quant_bin, ylim = c(-0.15, 0.15))
par(xpd=FALSE)
abline(h = 0, lty = 2)
par(xpd=TRUE)
mtext(expression("1"^"st"), side = 1, line = 1, at = 1, cex = 0.6)
mtext(expression("2"^"nd"), side = 1, line = 1, at = 2, cex = 0.6)
mtext(expression("3"^"rd"), side = 1, line = 1, at = 3, cex = 0.6)
mtext(expression("4"^"th"), side = 1, line = 1, at = 4, cex = 0.6)
mtext(expression("5"^"th"), side = 1, line = 1, at = 5, cex = 0.6)
mtext(expression("6"^"th"), side = 1, line = 1, at = 6, cex = 0.6)
mtext(expression("7"^"th"), side = 1, line = 1, at = 7, cex = 0.6)
mtext(expression("8"^"th"), side = 1, line = 1, at = 8, cex = 0.6)
mtext(expression("9"^"th"), side = 1, line = 1, at = 9, cex = 0.6)
mtext(expression("10"^"th"), side = 1, line = 1, at = 10, cex = 0.6)
mtext(expression("A"[4]), side = 3, line = 0.7, at = -2.5, cex = 0.9)
par(mar = c(3, 2.4, 2.5, 0.5))

function_delta_frequ_plot(Generation_1250_neg_sel_s_0.005_bottleneck_delta_frequ_0.1_quant_bin, ylim = c(-0.2, 0.2), my_ylab = "")
par(xpd=FALSE)
abline(h = 0, lty = 2)
par(xpd=TRUE)
mtext(expression("1"^"st"), side = 1, line = 1, at = 1, cex = 0.6)
mtext(expression("2"^"nd"), side = 1, line = 1, at = 2, cex = 0.6)
mtext(expression("3"^"rd"), side = 1, line = 1, at = 3, cex = 0.6)
mtext(expression("4"^"th"), side = 1, line = 1, at = 4, cex = 0.6)
mtext(expression("5"^"th"), side = 1, line = 1, at = 5, cex = 0.6)
mtext(expression("6"^"th"), side = 1, line = 1, at = 6, cex = 0.6)
mtext(expression("7"^"th"), side = 1, line = 1, at = 7, cex = 0.6)
mtext(expression("8"^"th"), side = 1, line = 1, at = 8, cex = 0.6)
mtext(expression("9"^"th"), side = 1, line = 1, at = 9, cex = 0.6)
mtext(expression("10"^"th"), side = 1, line = 1, at = 10, cex = 0.6)
mtext(expression("B"[4]), side = 3, line = 0.7, at = -2, cex = 0.9)

function_delta_frequ_plot(Generation_1250_neutral_TE_burst_delta_frequ_0.1_quant_bin, ylim = c(-0.04, 0.04), my_ylab = "")
par(xpd=FALSE)
abline(h = 0, lty = 2)
par(xpd=TRUE)
mtext(expression("1"^"st"), side = 1, line = 1, at = 1, cex = 0.6)
mtext(expression("2"^"nd"), side = 1, line = 1, at = 2, cex = 0.6)
mtext(expression("3"^"rd"), side = 1, line = 1, at = 3, cex = 0.6)
mtext(expression("4"^"th"), side = 1, line = 1, at = 4, cex = 0.6)
mtext(expression("5"^"th"), side = 1, line = 1, at = 5, cex = 0.6)
mtext(expression("6"^"th"), side = 1, line = 1, at = 6, cex = 0.6)
mtext(expression("7"^"th"), side = 1, line = 1, at = 7, cex = 0.6)
mtext(expression("8"^"th"), side = 1, line = 1, at = 8, cex = 0.6)
mtext(expression("9"^"th"), side = 1, line = 1, at = 9, cex = 0.6)
mtext(expression("10"^"th"), side = 1, line = 1, at = 10, cex = 0.6)
mtext(expression("C"[4]), side = 3, line = 0.7, at = -2, cex = 0.9)

function_delta_frequ_plot(Generation_1250_neg_sel_s_0.005_TE_burst_delta_frequ_0.1_quant_bin, ylim = c(-0.16, 0.05), my_ylab = "")
par(xpd=FALSE)
abline(h = 0, lty = 2)
par(xpd=TRUE)
mtext(expression("1"^"st"), side = 1, line = 1, at = 1, cex = 0.6)
mtext(expression("2"^"nd"), side = 1, line = 1, at = 2, cex = 0.6)
mtext(expression("3"^"rd"), side = 1, line = 1, at = 3, cex = 0.6)
mtext(expression("4"^"th"), side = 1, line = 1, at = 4, cex = 0.6)
mtext(expression("5"^"th"), side = 1, line = 1, at = 5, cex = 0.6)
mtext(expression("6"^"th"), side = 1, line = 1, at = 6, cex = 0.6)
mtext(expression("7"^"th"), side = 1, line = 1, at = 7, cex = 0.6)
mtext(expression("8"^"th"), side = 1, line = 1, at = 8, cex = 0.6)
mtext(expression("9"^"th"), side = 1, line = 1, at = 9, cex = 0.6)
mtext(expression("10"^"th"), side = 1, line = 1, at = 10, cex = 0.6)
mtext(expression("D"[4]), side = 3, line = 0.7, at = -2, cex = 0.9)

dev.off()

# plot delta freqeuncy
pdf("/Users/roberthorvath/Desktop/Figure Sup 4.pdf", height = 7)
layout(mat = matrix(c(seq(1,13,4), seq(2,14,4), seq(3,15,4), seq(4,16,4) ), nrow = 4, ncol = 4), heights = c(1, 1, 1, 1.2), widths = c(1, 0.9, 0.9, 0.9))
par(xpd=TRUE, mar = c(0.5, 3.5, 2.5, 0.5), mgp = c(2, 1, 0))

function_delta_frequ_plot(Generation_300_neutral_bottleneck_delta_frequ_0.1_quant_bin, ylim = c(-0.15, 0.15), my_xlab = "")
par(xpd=FALSE)
abline(h = 0, lty = 2)
par(xpd=TRUE)
mtext(expression("A"[1]), side = 3, line = 0.7, at = -2.5, cex = 0.9)
par(xpd=TRUE, mar = c(0.5, 2.4, 2.5, 0.5), mgp = c(2, 1, 0))

function_delta_frequ_plot(Generation_300_neg_sel_s_0.005_bottleneck_delta_frequ_0.1_quant_bin, ylim = c(-0.2, 0.2), my_xlab = "", my_ylab = "")
par(xpd=FALSE)
abline(h = 0, lty = 2)
par(xpd=TRUE)
mtext(expression("B"[1]), side = 3, line = 0.7, at = -2, cex = 0.9)

function_delta_frequ_plot(Generation_300_neutral_TE_burst_delta_frequ_0.1_quant_bin, ylim = c(-0.04, 0.04), my_xlab = "", my_ylab = "")
par(xpd=FALSE)
abline(h = 0, lty = 2)
par(xpd=TRUE)
mtext(expression("C"[1]), side = 3, line = 0.7, at = -2, cex = 0.9)

function_delta_frequ_plot(Generation_300_neg_sel_s_0.005_TE_burst_delta_frequ_0.1_quant_bin, ylim = c(-0.16, 0.05), my_xlab = "", my_ylab = "")
par(xpd=FALSE)
abline(h = 0, lty = 2)
par(xpd=TRUE)
mtext(expression("D"[1]), side = 3, line = 0.7, at = -2, cex = 0.9)
par(xpd=TRUE, mar = c(0.5, 3.5, 2.5, 0.5), mgp = c(2, 1, 0))

function_delta_frequ_plot(Generation_500_neutral_bottleneck_delta_frequ_0.1_quant_bin, ylim = c(-0.15, 0.15), my_xlab = "")
par(xpd=FALSE)
abline(h = 0, lty = 2)
par(xpd=TRUE)
mtext(expression("A"[2]), side = 3, line = 0.7, at = -2.5, cex = 0.9)
par(xpd=TRUE, mar = c(0.5, 2.4, 2.5, 0.5), mgp = c(2, 1, 0))

function_delta_frequ_plot(Generation_500_neg_sel_s_0.005_bottleneck_delta_frequ_0.1_quant_bin, ylim = c(-0.2, 0.2), my_xlab = "", my_ylab = "")
par(xpd=FALSE)
abline(h = 0, lty = 2)
par(xpd=TRUE)
mtext(expression("B"[2]), side = 3, line = 0.7, at = -2, cex = 0.9)

function_delta_frequ_plot(Generation_500_neutral_TE_burst_delta_frequ_0.1_quant_bin, ylim = c(-0.04, 0.04), my_xlab = "", my_ylab = "")
par(xpd=FALSE)
abline(h = 0, lty = 2)
par(xpd=TRUE)
mtext(expression("C"[2]), side = 3, line = 0.7, at = -2, cex = 0.9)

function_delta_frequ_plot(Generation_500_neg_sel_s_0.005_TE_burst_delta_frequ_0.1_quant_bin, ylim = c(-0.16, 0.05), my_xlab = "", my_ylab = "")
par(xpd=FALSE)
abline(h = 0, lty = 2)
par(xpd=TRUE)
mtext(expression("D"[2]), side = 3, line = 0.7, at = -2, cex = 0.9)
par(xpd=TRUE, mar = c(0.5, 3.5, 2.5, 0.5), mgp = c(2, 1, 0))

function_delta_frequ_plot(Generation_2000_neutral_bottleneck_delta_frequ_0.1_quant_bin, ylim = c(-0.15, 0.15), my_xlab = "")
par(xpd=FALSE)
abline(h = 0, lty = 2)
par(xpd=TRUE)
mtext(expression("A"[3]), side = 3, line = 0.7, at = -2.5, cex = 0.9)
par(xpd=TRUE, mar = c(0.5, 2.4, 2.5, 0.5), mgp = c(2, 1, 0))

function_delta_frequ_plot(Generation_2000_neg_sel_s_0.005_bottleneck_delta_frequ_0.1_quant_bin, ylim = c(-0.2, 0.2), my_xlab = "" , my_ylab = "")
par(xpd=FALSE)
abline(h = 0, lty = 2)
par(xpd=TRUE)
mtext(expression("B"[3]), side = 3, line = 0.7, at = -2, cex = 0.9)

function_delta_frequ_plot(Generation_2000_neutral_TE_burst_delta_frequ_0.1_quant_bin, ylim = c(-0.04, 0.04), my_xlab = "", my_ylab = "")
par(xpd=FALSE)
abline(h = 0, lty = 2)
par(xpd=TRUE)
mtext(expression("C"[3]), side = 3, line = 0.7, at = -2, cex = 0.9)

function_delta_frequ_plot(Generation_2000_neg_sel_s_0.005_TE_burst_delta_frequ_0.1_quant_bin, ylim = c(-0.16, 0.05), my_xlab = "", my_ylab = "")
par(xpd=FALSE)
abline(h = 0, lty = 2)
par(xpd=TRUE)
mtext(expression("D"[3]), side = 3, line = 0.7, at = -2, cex = 0.9)
par(mar = c(3, 3.5, 2.5, 0.5))

function_delta_frequ_plot(Generation_5000_neutral_bottleneck_delta_frequ_0.1_quant_bin, ylim = c(-0.15, 0.15))
par(xpd=FALSE)
abline(h = 0, lty = 2)
par(xpd=TRUE)
mtext(expression("1"^"st"), side = 1, line = 1, at = 1, cex = 0.6)
mtext(expression("2"^"nd"), side = 1, line = 1, at = 2, cex = 0.6)
mtext(expression("3"^"rd"), side = 1, line = 1, at = 3, cex = 0.6)
mtext(expression("4"^"th"), side = 1, line = 1, at = 4, cex = 0.6)
mtext(expression("5"^"th"), side = 1, line = 1, at = 5, cex = 0.6)
mtext(expression("6"^"th"), side = 1, line = 1, at = 6, cex = 0.6)
mtext(expression("7"^"th"), side = 1, line = 1, at = 7, cex = 0.6)
mtext(expression("8"^"th"), side = 1, line = 1, at = 8, cex = 0.6)
mtext(expression("9"^"th"), side = 1, line = 1, at = 9, cex = 0.6)
mtext(expression("10"^"th"), side = 1, line = 1, at = 10, cex = 0.6)
mtext(expression("A"[4]), side = 3, line = 0.7, at = -2.5, cex = 0.9)
par(mar = c(3, 2.4, 2.5, 0.5))

function_delta_frequ_plot(Generation_5000_neg_sel_s_0.005_bottleneck_delta_frequ_0.1_quant_bin, ylim = c(-0.2, 0.2), my_ylab = "")
par(xpd=FALSE)
abline(h = 0, lty = 2)
par(xpd=TRUE)
mtext(expression("1"^"st"), side = 1, line = 1, at = 1, cex = 0.6)
mtext(expression("2"^"nd"), side = 1, line = 1, at = 2, cex = 0.6)
mtext(expression("3"^"rd"), side = 1, line = 1, at = 3, cex = 0.6)
mtext(expression("4"^"th"), side = 1, line = 1, at = 4, cex = 0.6)
mtext(expression("5"^"th"), side = 1, line = 1, at = 5, cex = 0.6)
mtext(expression("6"^"th"), side = 1, line = 1, at = 6, cex = 0.6)
mtext(expression("7"^"th"), side = 1, line = 1, at = 7, cex = 0.6)
mtext(expression("8"^"th"), side = 1, line = 1, at = 8, cex = 0.6)
mtext(expression("9"^"th"), side = 1, line = 1, at = 9, cex = 0.6)
mtext(expression("10"^"th"), side = 1, line = 1, at = 10, cex = 0.6)
mtext(expression("B"[4]), side = 3, line = 0.7, at = -2, cex = 0.9)

function_delta_frequ_plot(Generation_5000_neutral_TE_burst_delta_frequ_0.1_quant_bin, ylim = c(-0.04, 0.04), my_ylab = "")
par(xpd=FALSE)
abline(h = 0, lty = 2)
par(xpd=TRUE)
mtext(expression("1"^"st"), side = 1, line = 1, at = 1, cex = 0.6)
mtext(expression("2"^"nd"), side = 1, line = 1, at = 2, cex = 0.6)
mtext(expression("3"^"rd"), side = 1, line = 1, at = 3, cex = 0.6)
mtext(expression("4"^"th"), side = 1, line = 1, at = 4, cex = 0.6)
mtext(expression("5"^"th"), side = 1, line = 1, at = 5, cex = 0.6)
mtext(expression("6"^"th"), side = 1, line = 1, at = 6, cex = 0.6)
mtext(expression("7"^"th"), side = 1, line = 1, at = 7, cex = 0.6)
mtext(expression("8"^"th"), side = 1, line = 1, at = 8, cex = 0.6)
mtext(expression("9"^"th"), side = 1, line = 1, at = 9, cex = 0.6)
mtext(expression("10"^"th"), side = 1, line = 1, at = 10, cex = 0.6)
mtext(expression("C"[4]), side = 3, line = 0.7, at = -2, cex = 0.9)

function_delta_frequ_plot(Generation_5000_neg_sel_s_0.005_TE_burst_delta_frequ_0.1_quant_bin, ylim = c(-0.16, 0.05), my_ylab = "")
par(xpd=FALSE)
abline(h = 0, lty = 2)
par(xpd=TRUE)
mtext(expression("1"^"st"), side = 1, line = 1, at = 1, cex = 0.6)
mtext(expression("2"^"nd"), side = 1, line = 1, at = 2, cex = 0.6)
mtext(expression("3"^"rd"), side = 1, line = 1, at = 3, cex = 0.6)
mtext(expression("4"^"th"), side = 1, line = 1, at = 4, cex = 0.6)
mtext(expression("5"^"th"), side = 1, line = 1, at = 5, cex = 0.6)
mtext(expression("6"^"th"), side = 1, line = 1, at = 6, cex = 0.6)
mtext(expression("7"^"th"), side = 1, line = 1, at = 7, cex = 0.6)
mtext(expression("8"^"th"), side = 1, line = 1, at = 8, cex = 0.6)
mtext(expression("9"^"th"), side = 1, line = 1, at = 9, cex = 0.6)
mtext(expression("10"^"th"), side = 1, line = 1, at = 10, cex = 0.6)
mtext(expression("D"[4]), side = 3, line = 0.7, at = -2, cex = 0.9)

dev.off()

# Calculate the percentage of runs which showed a significant negative correlation between delta frequency and age
my_sig_percent_burst <- c(
  my_get_number_of_sig(Generation_0_neg_sel_s_0.01_TE_burst_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_50_neg_sel_s_0.01_TE_burst_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_250_neg_sel_s_0.01_TE_burst_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_300_neg_sel_s_0.01_TE_burst_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_500_neg_sel_s_0.01_TE_burst_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_750_neg_sel_s_0.01_TE_burst_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_1250_neg_sel_s_0.01_TE_burst_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_2000_neg_sel_s_0.01_TE_burst_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_5000_neg_sel_s_0.01_TE_burst_delta_frequ_0.1_quant_bin),
  
  my_get_number_of_sig(Generation_0_neg_sel_s_0.005_TE_burst_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_50_neg_sel_s_0.005_TE_burst_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_250_neg_sel_s_0.005_TE_burst_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_300_neg_sel_s_0.005_TE_burst_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_500_neg_sel_s_0.005_TE_burst_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_750_neg_sel_s_0.005_TE_burst_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_1250_neg_sel_s_0.005_TE_burst_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_2000_neg_sel_s_0.005_TE_burst_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_5000_neg_sel_s_0.005_TE_burst_delta_frequ_0.1_quant_bin),
  
  my_get_number_of_sig(Generation_0_neg_sel_s_0.001_TE_burst_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_50_neg_sel_s_0.001_TE_burst_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_250_neg_sel_s_0.001_TE_burst_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_300_neg_sel_s_0.001_TE_burst_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_500_neg_sel_s_0.001_TE_burst_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_750_neg_sel_s_0.001_TE_burst_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_1250_neg_sel_s_0.001_TE_burst_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_2000_neg_sel_s_0.001_TE_burst_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_5000_neg_sel_s_0.001_TE_burst_delta_frequ_0.1_quant_bin),
  
  my_get_number_of_sig(Generation_0_neutral_TE_burst_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_50_neutral_TE_burst_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_250_neutral_TE_burst_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_300_neutral_TE_burst_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_500_neutral_TE_burst_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_750_neutral_TE_burst_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_1250_neutral_TE_burst_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_2000_neutral_TE_burst_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_5000_neutral_TE_burst_delta_frequ_0.1_quant_bin)
)

my_sig_percent_bot <- c(
  my_get_number_of_sig(Generation_0_neg_sel_s_0.01_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_50_neg_sel_s_0.01_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_250_neg_sel_s_0.01_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_300_neg_sel_s_0.01_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_500_neg_sel_s_0.01_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_750_neg_sel_s_0.01_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_1250_neg_sel_s_0.01_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_2000_neg_sel_s_0.01_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_5000_neg_sel_s_0.01_bottleneck_delta_frequ_0.1_quant_bin),
  
  my_get_number_of_sig(Generation_0_neg_sel_s_0.005_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_50_neg_sel_s_0.005_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_250_neg_sel_s_0.005_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_300_neg_sel_s_0.005_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_500_neg_sel_s_0.005_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_750_neg_sel_s_0.005_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_1250_neg_sel_s_0.005_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_2000_neg_sel_s_0.005_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_5000_neg_sel_s_0.005_bottleneck_delta_frequ_0.1_quant_bin),
  
  my_get_number_of_sig(Generation_0_neg_sel_s_0.001_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_50_neg_sel_s_0.001_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_250_neg_sel_s_0.001_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_300_neg_sel_s_0.001_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_500_neg_sel_s_0.001_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_750_neg_sel_s_0.001_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_1250_neg_sel_s_0.001_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_2000_neg_sel_s_0.001_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_5000_neg_sel_s_0.001_bottleneck_delta_frequ_0.1_quant_bin),
  
  my_get_number_of_sig(Generation_0_neutral_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_50_neutral_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_250_neutral_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_300_neutral_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_500_neutral_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_750_neutral_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_1250_neutral_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_2000_neutral_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_5000_neutral_bottleneck_delta_frequ_0.1_quant_bin)
)

my_sig_percent_bot_burst <- c(
  my_get_number_of_sig(Generation_0_neg_sel_s_0.01_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_50_neg_sel_s_0.01_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_250_neg_sel_s_0.01_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_300_neg_sel_s_0.01_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_500_neg_sel_s_0.01_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_750_neg_sel_s_0.01_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_1250_neg_sel_s_0.01_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_2000_neg_sel_s_0.01_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_5000_neg_sel_s_0.01_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin),
  
  my_get_number_of_sig(Generation_0_neg_sel_s_0.005_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_50_neg_sel_s_0.005_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_250_neg_sel_s_0.005_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_300_neg_sel_s_0.005_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_500_neg_sel_s_0.005_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_750_neg_sel_s_0.005_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_1250_neg_sel_s_0.005_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_2000_neg_sel_s_0.005_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_5000_neg_sel_s_0.005_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin),
  
  my_get_number_of_sig(Generation_0_neg_sel_s_0.001_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_50_neg_sel_s_0.001_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_250_neg_sel_s_0.001_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_300_neg_sel_s_0.001_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_500_neg_sel_s_0.001_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_750_neg_sel_s_0.001_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_1250_neg_sel_s_0.001_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_2000_neg_sel_s_0.001_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_5000_neg_sel_s_0.001_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin),
  
  my_get_number_of_sig(Generation_0_neutral_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_50_neutral_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_250_neutral_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_300_neutral_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_500_neutral_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_750_neutral_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_1250_neutral_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_2000_neutral_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin),
  my_get_number_of_sig(Generation_5000_neutral_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin)
)

# plot the heatmap
library("ggplot2") 
library("gridExtra") 
library("grid")

row <- factor(c(rep("-20", 9), rep("-10", 9), rep("-2", 9), rep("0", 9)),levels = c("-20", "-10", "-2", "0"))
col <- factor(rep(c("0", "50", "250", "300", "500", "750", "1250", "2000", "5000"), 4),levels = c("0", "50", "250", "300", "500", "750", "1250", "2000", "5000"))

my_sig_percent_bot.dat <- data.frame(row, col, my_sig_percent_bot)
my_sig_percent_burst.dat <- data.frame(row, col, my_sig_percent_burst)
my_sig_percent_bot_burst.dat <- data.frame(row, col, my_sig_percent_bot_burst)

p1 <- ggplot(my_sig_percent_bot.dat, aes(col,row)) +
  geom_tile(aes(fill = my_sig_percent_bot)) +
  scale_fill_gradient(low = "blue", high = "red", limits = c(0,100)) +
  ylab(expression(paste( "4", italic("N") [ italic("e")], italic("s") ))) +
  theme(axis.title.x = element_blank(), axis.text.x=element_blank(), legend.position = "none", panel.background = element_blank())
p2 <- ggplot(my_sig_percent_burst.dat, aes(col,row)) +
  geom_tile(aes(fill = my_sig_percent_burst)) +
  scale_fill_gradient(low = "blue", high = "red", name = "Percentage", limits = c(0,100)) +
  ylab(expression(paste( "4", italic("N") [ italic("e")], italic("s") ))) +
  theme(axis.title.x = element_blank(), axis.text.x=element_blank(), panel.background = element_blank())
p3 <- ggplot(my_sig_percent_bot_burst.dat, aes(col,row)) +
  geom_tile(aes(fill = my_sig_percent_bot_burst)) +
  scale_fill_gradient(low = "blue", high = "red", limits = c(0,100)) +
  xlab("Generations") + ylab(expression(paste( "4", italic("N") [ italic("e")], italic("s") ))) +
  theme(legend.position = "none", panel.background = element_blank()) 

pdf("/Users/roberthorvath/Desktop/Figure 4.pdf")
grid.arrange(p1, p2, p3, nrow = 3, heights = c(1,1,1.15), layout_matrix = rbind(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, NA, NA),
                                                                                c(2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2),
                                                                                c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, NA, NA)))
grid.text("A", x = 0.02, y = 0.985)
grid.text("B", x = 0.02, y = 0.665)
grid.text("C", x = 0.02, y = 0.35)

dev.off()


# Table S1 
# calculate the percentage of runs which showed a significant negative correlation between delta frequency and age for the runs not included in the heatmap
my_get_number_of_sig(Generation_0_neg_sel_s_0.0001_TE_burst_delta_frequ_0.1_quant_bin)
my_get_number_of_sig(Generation_50_neg_sel_s_0.0001_TE_burst_delta_frequ_0.1_quant_bin)
my_get_number_of_sig(Generation_250_neg_sel_s_0.0001_TE_burst_delta_frequ_0.1_quant_bin)
my_get_number_of_sig(Generation_300_neg_sel_s_0.0001_TE_burst_delta_frequ_0.1_quant_bin)
my_get_number_of_sig(Generation_500_neg_sel_s_0.0001_TE_burst_delta_frequ_0.1_quant_bin)
my_get_number_of_sig(Generation_750_neg_sel_s_0.0001_TE_burst_delta_frequ_0.1_quant_bin)
my_get_number_of_sig(Generation_1250_neg_sel_s_0.0001_TE_burst_delta_frequ_0.1_quant_bin)
my_get_number_of_sig(Generation_2000_neg_sel_s_0.0001_TE_burst_delta_frequ_0.1_quant_bin)
my_get_number_of_sig(Generation_5000_neg_sel_s_0.0001_TE_burst_delta_frequ_0.1_quant_bin)

my_get_number_of_sig(Generation_0_neg_sel_s_0.0001_bottleneck_delta_frequ_0.1_quant_bin)
my_get_number_of_sig(Generation_50_neg_sel_s_0.0001_bottleneck_delta_frequ_0.1_quant_bin)
my_get_number_of_sig(Generation_250_neg_sel_s_0.0001_bottleneck_delta_frequ_0.1_quant_bin)
my_get_number_of_sig(Generation_300_neg_sel_s_0.0001_bottleneck_delta_frequ_0.1_quant_bin)
my_get_number_of_sig(Generation_500_neg_sel_s_0.0001_bottleneck_delta_frequ_0.1_quant_bin)
my_get_number_of_sig(Generation_750_neg_sel_s_0.0001_bottleneck_delta_frequ_0.1_quant_bin)
my_get_number_of_sig(Generation_1250_neg_sel_s_0.0001_bottleneck_delta_frequ_0.1_quant_bin)
my_get_number_of_sig(Generation_2000_neg_sel_s_0.0001_bottleneck_delta_frequ_0.1_quant_bin)
my_get_number_of_sig(Generation_5000_neg_sel_s_0.0001_bottleneck_delta_frequ_0.1_quant_bin)

my_get_number_of_sig(Generation_0_neg_sel_s_0.0001_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin)
my_get_number_of_sig(Generation_50_neg_sel_s_0.0001_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin)
my_get_number_of_sig(Generation_250_neg_sel_s_0.0001_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin)
my_get_number_of_sig(Generation_300_neg_sel_s_0.0001_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin)
my_get_number_of_sig(Generation_500_neg_sel_s_0.0001_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin)
my_get_number_of_sig(Generation_750_neg_sel_s_0.0001_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin)
my_get_number_of_sig(Generation_1250_neg_sel_s_0.0001_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin)
my_get_number_of_sig(Generation_2000_neg_sel_s_0.0001_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin)
my_get_number_of_sig(Generation_5000_neg_sel_s_0.0001_TE_burst_and_bottleneck_delta_frequ_0.1_quant_bin)


# Table S2
# get number of segregating TE insertions in the population
quantile(TE_generation_0_neutral_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_0_neg_sel_s_0.0001_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_0_neg_sel_s_0.001_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_0_neg_sel_s_0.005_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_0_neg_sel_s_0.01_bottleneck_AFS[41,], c(0,1))

quantile(TE_generation_50_neutral_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_50_neg_sel_s_0.0001_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_50_neg_sel_s_0.001_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_50_neg_sel_s_0.005_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_50_neg_sel_s_0.01_bottleneck_AFS[41,], c(0,1))

quantile(TE_generation_250_neutral_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_250_neg_sel_s_0.0001_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_250_neg_sel_s_0.001_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_250_neg_sel_s_0.005_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_250_neg_sel_s_0.01_bottleneck_AFS[41,], c(0,1))

quantile(TE_generation_300_neutral_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_300_neg_sel_s_0.0001_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_300_neg_sel_s_0.001_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_300_neg_sel_s_0.005_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_300_neg_sel_s_0.01_bottleneck_AFS[41,], c(0,1))

quantile(TE_generation_500_neutral_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_500_neg_sel_s_0.0001_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_500_neg_sel_s_0.001_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_500_neg_sel_s_0.005_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_500_neg_sel_s_0.01_bottleneck_AFS[41,], c(0,1))

quantile(TE_generation_750_neutral_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_750_neg_sel_s_0.0001_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_750_neg_sel_s_0.001_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_750_neg_sel_s_0.005_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_750_neg_sel_s_0.01_bottleneck_AFS[41,], c(0,1))

quantile(TE_generation_1250_neutral_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_1250_neg_sel_s_0.0001_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_1250_neg_sel_s_0.001_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_1250_neg_sel_s_0.005_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_1250_neg_sel_s_0.01_bottleneck_AFS[41,], c(0,1))

quantile(TE_generation_2000_neutral_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_2000_neg_sel_s_0.0001_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_2000_neg_sel_s_0.001_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_2000_neg_sel_s_0.005_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_2000_neg_sel_s_0.01_bottleneck_AFS[41,], c(0,1))

quantile(TE_generation_5000_neutral_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_5000_neg_sel_s_0.0001_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_5000_neg_sel_s_0.001_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_5000_neg_sel_s_0.005_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_5000_neg_sel_s_0.01_bottleneck_AFS[41,], c(0,1))

quantile(TE_generation_0_neutral_TE_burst_AFS[41,], c(0,1))
quantile(TE_generation_0_neg_sel_s_0.0001_TE_burst_AFS[41,], c(0,1))
quantile(TE_generation_0_neg_sel_s_0.001_TE_burst_AFS[41,], c(0,1))
quantile(TE_generation_0_neg_sel_s_0.005_TE_burst_AFS[41,], c(0,1))
quantile(TE_generation_0_neg_sel_s_0.01_TE_burst_AFS[41,], c(0,1))

quantile(TE_generation_50_neutral_TE_burst_AFS[41,], c(0,1))
quantile(TE_generation_50_neg_sel_s_0.0001_TE_burst_AFS[41,], c(0,1))
quantile(TE_generation_50_neg_sel_s_0.001_TE_burst_AFS[41,], c(0,1))
quantile(TE_generation_50_neg_sel_s_0.005_TE_burst_AFS[41,], c(0,1))
quantile(TE_generation_50_neg_sel_s_0.01_TE_burst_AFS[41,], c(0,1))

quantile(TE_generation_250_neutral_TE_burst_AFS[41,], c(0,1))
quantile(TE_generation_250_neg_sel_s_0.0001_TE_burst_AFS[41,], c(0,1))
quantile(TE_generation_250_neg_sel_s_0.001_TE_burst_AFS[41,], c(0,1))
quantile(TE_generation_250_neg_sel_s_0.005_TE_burst_AFS[41,], c(0,1))
quantile(TE_generation_250_neg_sel_s_0.01_TE_burst_AFS[41,], c(0,1))

quantile(TE_generation_300_neutral_TE_burst_AFS[41,], c(0,1))
quantile(TE_generation_300_neg_sel_s_0.0001_TE_burst_AFS[41,], c(0,1))
quantile(TE_generation_300_neg_sel_s_0.001_TE_burst_AFS[41,], c(0,1))
quantile(TE_generation_300_neg_sel_s_0.005_TE_burst_AFS[41,], c(0,1))
quantile(TE_generation_300_neg_sel_s_0.01_TE_burst_AFS[41,], c(0,1))

quantile(TE_generation_500_neutral_TE_burst_AFS[41,], c(0,1))
quantile(TE_generation_500_neg_sel_s_0.0001_TE_burst_AFS[41,], c(0,1))
quantile(TE_generation_500_neg_sel_s_0.001_TE_burst_AFS[41,], c(0,1))
quantile(TE_generation_500_neg_sel_s_0.005_TE_burst_AFS[41,], c(0,1))
quantile(TE_generation_500_neg_sel_s_0.01_TE_burst_AFS[41,], c(0,1))

quantile(TE_generation_750_neutral_TE_burst_AFS[41,], c(0,1))
quantile(TE_generation_750_neg_sel_s_0.0001_TE_burst_AFS[41,], c(0,1))
quantile(TE_generation_750_neg_sel_s_0.001_TE_burst_AFS[41,], c(0,1))
quantile(TE_generation_750_neg_sel_s_0.005_TE_burst_AFS[41,], c(0,1))
quantile(TE_generation_750_neg_sel_s_0.01_TE_burst_AFS[41,], c(0,1))

quantile(TE_generation_1250_neutral_TE_burst_AFS[41,], c(0,1))
quantile(TE_generation_1250_neg_sel_s_0.0001_TE_burst_AFS[41,], c(0,1))
quantile(TE_generation_1250_neg_sel_s_0.001_TE_burst_AFS[41,], c(0,1))
quantile(TE_generation_1250_neg_sel_s_0.005_TE_burst_AFS[41,], c(0,1))
quantile(TE_generation_1250_neg_sel_s_0.01_TE_burst_AFS[41,], c(0,1))

quantile(TE_generation_2000_neutral_TE_burst_AFS[41,], c(0,1))
quantile(TE_generation_2000_neg_sel_s_0.0001_TE_burst_AFS[41,], c(0,1))
quantile(TE_generation_2000_neg_sel_s_0.001_TE_burst_AFS[41,], c(0,1))
quantile(TE_generation_2000_neg_sel_s_0.005_TE_burst_AFS[41,], c(0,1))
quantile(TE_generation_2000_neg_sel_s_0.01_TE_burst_AFS[41,], c(0,1))

quantile(TE_generation_5000_neutral_TE_burst_AFS[41,], c(0,1))
quantile(TE_generation_5000_neg_sel_s_0.0001_TE_burst_AFS[41,], c(0,1))
quantile(TE_generation_5000_neg_sel_s_0.001_TE_burst_AFS[41,], c(0,1))
quantile(TE_generation_5000_neg_sel_s_0.005_TE_burst_AFS[41,], c(0,1))
quantile(TE_generation_5000_neg_sel_s_0.01_TE_burst_AFS[41,], c(0,1))

quantile(TE_generation_0_neutral_TE_burst_and_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_0_neg_sel_s_0.0001_TE_burst_and_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_0_neg_sel_s_0.001_TE_burst_and_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_0_neg_sel_s_0.005_TE_burst_and_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_0_neg_sel_s_0.01_TE_burst_and_bottleneck_AFS[41,], c(0,1))

quantile(TE_generation_50_neutral_TE_burst_and_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_50_neg_sel_s_0.0001_TE_burst_and_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_50_neg_sel_s_0.001_TE_burst_and_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_50_neg_sel_s_0.005_TE_burst_and_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_50_neg_sel_s_0.01_TE_burst_and_bottleneck_AFS[41,], c(0,1))

quantile(TE_generation_250_neutral_TE_burst_and_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_250_neg_sel_s_0.0001_TE_burst_and_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_250_neg_sel_s_0.001_TE_burst_and_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_250_neg_sel_s_0.005_TE_burst_and_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_250_neg_sel_s_0.01_TE_burst_and_bottleneck_AFS[41,], c(0,1))

quantile(TE_generation_300_neutral_TE_burst_and_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_300_neg_sel_s_0.0001_TE_burst_and_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_300_neg_sel_s_0.001_TE_burst_and_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_300_neg_sel_s_0.005_TE_burst_and_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_300_neg_sel_s_0.01_TE_burst_and_bottleneck_AFS[41,], c(0,1))

quantile(TE_generation_500_neutral_TE_burst_and_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_500_neg_sel_s_0.0001_TE_burst_and_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_500_neg_sel_s_0.001_TE_burst_and_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_500_neg_sel_s_0.005_TE_burst_and_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_500_neg_sel_s_0.01_TE_burst_and_bottleneck_AFS[41,], c(0,1))

quantile(TE_generation_750_neutral_TE_burst_and_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_750_neg_sel_s_0.0001_TE_burst_and_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_750_neg_sel_s_0.001_TE_burst_and_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_750_neg_sel_s_0.005_TE_burst_and_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_750_neg_sel_s_0.01_TE_burst_and_bottleneck_AFS[41,], c(0,1))

quantile(TE_generation_1250_neutral_TE_burst_and_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_1250_neg_sel_s_0.0001_TE_burst_and_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_1250_neg_sel_s_0.001_TE_burst_and_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_1250_neg_sel_s_0.005_TE_burst_and_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_1250_neg_sel_s_0.01_TE_burst_and_bottleneck_AFS[41,], c(0,1))

quantile(TE_generation_2000_neutral_TE_burst_and_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_2000_neg_sel_s_0.0001_TE_burst_and_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_2000_neg_sel_s_0.001_TE_burst_and_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_2000_neg_sel_s_0.005_TE_burst_and_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_2000_neg_sel_s_0.01_TE_burst_and_bottleneck_AFS[41,], c(0,1))

quantile(TE_generation_5000_neutral_TE_burst_and_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_5000_neg_sel_s_0.0001_TE_burst_and_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_5000_neg_sel_s_0.001_TE_burst_and_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_5000_neg_sel_s_0.005_TE_burst_and_bottleneck_AFS[41,], c(0,1))
quantile(TE_generation_5000_neg_sel_s_0.01_TE_burst_and_bottleneck_AFS[41,], c(0,1))


# plot the SFS of neutral sites
pdf("/Users/roberthorvath/Desktop/Figure Sup 5.pdf", height = 4 )
par(xpd=TRUE)

barplot(as.vector(rbind(
  as.numeric(apply(SNP_generation_0_neutral_bottleneck_SFS, 1, mean)),
  as.numeric(apply(SNP_generation_0_neg_sel_s_0.001_bottleneck_SFS, 1, mean)),
  as.numeric(apply(SNP_generation_0_neg_sel_s_0.005_bottleneck_SFS, 1, mean)),
  as.numeric(apply(SNP_generation_0_neg_sel_s_0.01_bottleneck_SFS, 1, mean)),
  as.numeric(apply(SNP_generation_0_neutral_TE_burst_SFS, 1, mean)),
  as.numeric(apply(SNP_generation_0_neg_sel_s_0.001_TE_burst_SFS, 1, mean)),
  as.numeric(apply(SNP_generation_0_neg_sel_s_0.005_TE_burst_SFS, 1, mean)),
  as.numeric(apply(SNP_generation_0_neg_sel_s_0.01_TE_burst_SFS, 1, mean)))), 
  col = c(8, 5:7, 1:4), space = c(1, 0, 0, 0, 0, 0, 0, 0), 
  ylim = c(0, 0.65), xlab = "Allele freqeuncy", ylab = "Proportion"
)
legend(50, 0.7, legend = c("Bottleneck S=0", "Bottleneck S=-2", "Bottleneck S=-10", "Bottleneck S=-20"), pch = 15, col = c(8, 5:7), bty = "n")
legend(130, 0.7, legend = c("TE burst S=0", "TE burst S=-2", "TE burst S=-10", "TE burst S=-20"), pch = 15, col = 1:4, bty = "n")
axis(1, seq(0.5, 180.5, 9), seq(0, 1, 0.05))

for (i in 1:20) {
  arrows(i + 0.5 + (i -1)*8, as.numeric(quantile(SNP_generation_0_neutral_bottleneck_SFS[i,], 1)), i + 0.5 + (i -1)*8, as.numeric(quantile(SNP_generation_0_neutral_bottleneck_SFS[i,], 0)), angle = 90, code = 3, length = 0.01)
  arrows(i + 1.5 + (i -1)*8, as.numeric(quantile(SNP_generation_0_neg_sel_s_0.001_bottleneck_SFS[i,], 1)), i + 1.5 + (i -1)*8, as.numeric(quantile(SNP_generation_0_neg_sel_s_0.001_bottleneck_SFS[i,], 0)), angle = 90, code = 3, length = 0.01)
  arrows(i + 2.5 + (i -1)*8, as.numeric(quantile(SNP_generation_0_neg_sel_s_0.005_bottleneck_SFS[i,], 1)), i + 2.5 + (i -1)*8, as.numeric(quantile(SNP_generation_0_neg_sel_s_0.005_bottleneck_SFS[i,], 0)), angle = 90, code = 3, length = 0.01)
  arrows(i + 3.5 + (i -1)*8, as.numeric(quantile(SNP_generation_0_neg_sel_s_0.01_bottleneck_SFS[i,], 1)), i + 3.5 + (i -1)*8, as.numeric(quantile(SNP_generation_0_neg_sel_s_0.01_bottleneck_SFS[i,], 0)), angle = 90, code = 3, length = 0.01)
  arrows(i + 4.5 + (i -1)*8, as.numeric(quantile(SNP_generation_0_neutral_TE_burst_SFS[i,], 1)), i + 4.5 + (i -1)*8, as.numeric(quantile(SNP_generation_0_neutral_TE_burst_SFS[i,], 0)), angle = 90, code = 3, length = 0.01)
  arrows(i + 5.5 + (i -1)*8, as.numeric(quantile(SNP_generation_0_neg_sel_s_0.001_TE_burst_SFS[i,], 1)), i + 5.5 + (i -1)*8, as.numeric(quantile(SNP_generation_0_neg_sel_s_0.001_TE_burst_SFS[i,], 0)), angle = 90, code = 3, length = 0.01)
  arrows(i + 6.5 + (i -1)*8, as.numeric(quantile(SNP_generation_0_neg_sel_s_0.005_TE_burst_SFS[i,], 1)), i + 6.5 + (i -1)*8, as.numeric(quantile(SNP_generation_0_neg_sel_s_0.005_TE_burst_SFS[i,], 0)), angle = 90, code = 3, length = 0.01)
  arrows(i + 7.5 + (i -1)*8, as.numeric(quantile(SNP_generation_0_neg_sel_s_0.01_TE_burst_SFS[i,], 1)), i + 7.5 + (i -1)*8, as.numeric(quantile(SNP_generation_0_neg_sel_s_0.01_TE_burst_SFS[i,], 0)), angle = 90, code = 3, length = 0.01)
}

dev.off()






