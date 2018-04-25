R code to run genomic sample size simulations
================
Elizabeth Flesch

-   [Code purpose](#code-purpose)
-   [Code instructions](#code-instructions)
-   [Within population Simulations](#within-population-simulations)
    -   [Mean kinship - within population](#mean-kinship---within-population)
    -   [IBS - within population](#ibs---within-population)
    -   [Example code to generate a basic boxplot of intrapopulation simulation results](#example-code-to-generate-a-basic-boxplot-of-intrapopulation-simulation-results)
-   [Between populations simulations](#between-populations-simulations)
    -   [Mean kinship - between populations](#mean-kinship---between-populations)
    -   [IBS - between populations](#ibs---between-populations)
    -   [Fst - between populations](#fst---between-populations)
-   [Report session information](#report-session-information)
-   [Citations](#citations)

Code purpose
============

This code is from the manuscript, "Evaluating sample size to estimate genetic management metrics in the genomics era" by Elizabeth P. Flesch, Jay J. Rotella, Jennifer M. Thomson, Tabitha A. Graves, and Robert A. Garrott. We used the following code (based on code by Nazareno et al. (2017)) to generate simulations regarding sample size for 4 different bighorn sheep populations. The provided code allows the user to produce simulations to evaluate kinship, identity by state, and *F*<sub>*S**T*</sub> estimates at different sample sizes using empirical genomic datasets.

Code instructions
=================

1.  Use a Linux system and ensure that the following software are loaded on your machine or server: PLINK v1.90 (Purcell et al., 2007) and KING v2.0 (Manichaikul et al., 2010).
2.  If running simulations on species other than one related to domestic sheep, search and replace "'--sheep'," in the code with your species of interest that is supported by PLINK.
3.  Modify function at the end of each code chunk with the required information. Your datafile must be saved in .ped format. No filtering is completed inside the simulation for the IBS and mean kinship simulations; all quality control filtering must be completed ahead of time. The Fst simulation includes filtering within it and records the number of SNPs retained at each filtering stage. These filtering thresholds can be modified if desired. We recommend testing the function with a small number of replicates (i.e. 100 replicates) first, as 10,000 replicates can take several hours depending on your dataset and system.
4.  Prior to running the code, set the working directory for where your .ped file is saved.
5.  Run code chunk.
6.  Each simulation chunk will save a file to the working directory containing a table of the simulation results. Each file is named based on the metric and dataset selected.

Within population Simulations
=============================

Mean kinship - within population
--------------------------------

``` r
## Set working directory setwd()

## Sample individuals from .fam file in R and generate new .fam file

### Function 1
resample <- function(n_ind, herd_name, fam_infile, fam_outfile) {
    # 1) Load file containing family data for all samples
    pop <- read.table(fam_infile)
    # 2) Sample from fam file a) Subset by selected herd
    popA <- pop[pop$V1 == herd_name, ]
    # b) Sample family file for required number of individuals
    dfA <- popA[sample(nrow(popA), size = n_ind, replace = FALSE), ]
    # 3) Write fam file with subsampled individuals
    write.table(dfA, file = fam_outfile, sep = " ", quote = F, row.names = F, col.names = F, 
        append = F)
}


## Function 2 gen_data should NOT include an extension, but the name of a ped file
resample_kinship <- function(reps, n_ind, herd_name, gen_data) {
    # Generate family file for input data
    system2("plink", args = c("--sheep", "--file", gen_data, "--out", "fam_infile", 
        "--make-bed", "--silent"))
    # set up an output list containing several dataframe that will contain the result
    # for the differentiation statistic; name the dataframe according to the
    # differentiation statistic; the first column of each dataframe will be the
    # sample size;
    res <- list()
    for (i in c("Kinship")) {
        res[[i]] <- data.frame(matrix(vector(), length(n_ind), reps))
        res[[i]] <- cbind(n_ind, res[[i]])
    }
    # 
    for (i in 1:length(n_ind)) {
        for (j in 1:reps) {
            # resample and save as a temporary outfile;- need to modify function to simply
            # provide data in the correct format
            resample(n_ind[i], herd_name, "fam_infile.fam", "resampled.fam")
            # Subset dataset using plink
            system2("plink", args = c("--sheep", "--file", gen_data, "--keep", "resampled.fam", 
                "--out", "filtered_data", "--make-bed", "--silent"))
            # Run KING to generate kinship files
            system2("king", args = c("-b", "filtered_data.bed", "--kinship", "--cpus", 
                "7", "--prefix", "filter_kinship", ">/dev/null"))
            # Extract the chosen differentiation statistics and save to the appropriate place
            # in the appropriate dataframe
            kinship_stats <- read.table(file = "filter_kinship.kin", header = TRUE)
            # 
            for (k in 1:length(res)) {
                stat_value_all <- kinship_stats$Kinship  #List all kinship pair values
                stat_value <- mean(stat_value_all)  #Calculate mean kinship
                res[[k]][i, j + 1] <- stat_value
            }
            # 
        }
    }
    # 
    write.table(res, file = paste0("Kinship_SD_", herd_name, "_", reps, "reps", "_", 
        gen_data, ".txt"), sep = " ", quote = F, row.names = F, col.names = T, append = F)
    return(as.data.frame(res))
    # 
}


## Run function

# Function arguments: reps (number of simulation replicates), n_ind (number of
# individuals drawn), herd_name (name of population as specified in file),
# gen_data (name of .ped file without extension)

resample_kinship(reps = 10000, n_ind = c(5, 10, 15, 20, 25), herd_name = "Fergus", 
    gen_data = "120_samples_table346_lowfilter")

resample_kinship(reps = 10000, n_ind = c(5, 10, 15, 20, 25), herd_name = "Glacier", 
    gen_data = "120_samples_table346_lowfilter")

resample_kinship(reps = 10000, n_ind = c(5, 10, 15, 20, 25), herd_name = "Taylor-Hilgard", 
    gen_data = "120_samples_table346_lowfilter")

resample_kinship(reps = 10000, n_ind = c(5, 10, 15, 20, 25), herd_name = "Wyoming", 
    gen_data = "120_samples_table346_lowfilter")


# Produce 30 sample results using 1 repetition Fergus
resample_kinship(reps = 1, n_ind = 30, herd_name = "Fergus", gen_data = "120_samples_table346_lowfilter")

# Glacier
resample_kinship(reps = 1, n_ind = 30, herd_name = "Glacier", gen_data = "120_samples_table346_lowfilter")

# Taylor-Hilgard
resample_kinship(reps = 1, n_ind = 30, herd_name = "Taylor-Hilgard", gen_data = "120_samples_table346_lowfilter")

# Wyoming
resample_kinship(reps = 1, n_ind = 30, herd_name = "Wyoming", gen_data = "120_samples_table346_lowfilter")
```

IBS - within population
-----------------------

``` r
# *****************************************************

# IBS Simulation within herds using KING

# **********************************************

## Set working directory Location on Linux laptop setwd()

## Sample individuals from .fam file in R and generate new .fam file

### Function 1-working
resample <- function(n_ind, herd_name, fam_infile, fam_outfile) {
    # 1) Load file containing family data for all samples
    pop <- read.table(fam_infile)
    # 2) Sample from fam file a) Subset by selected herd
    popA <- pop[pop$V1 == herd_name, ]
    # b) Sample family file for required number of individuals
    dfA <- popA[sample(nrow(popA), size = n_ind, replace = FALSE), ]
    # 3) Write fam file with subsampled individuals
    write.table(dfA, file = fam_outfile, sep = " ", quote = F, row.names = F, col.names = F, 
        append = F)
}


## Function 2 gen_data should NOT include an extension, but the name of a ped file
resample_kinship <- function(reps, n_ind, herd_name, gen_data) {
    # Generate family file for input data
    system2("plink", args = c("--sheep", "--file", gen_data, "--out", "fam_infile", 
        "--make-bed", "--silent"))
    # set up an output list containing several dataframe that will contain the result
    # for the differentiation statistic; name the dataframe according to the
    # differentiation statistic; the first column of each dataframe will be the
    # sample size;
    res <- list()
    for (i in c("IBS")) {
        res[[i]] <- data.frame(matrix(vector(), length(n_ind), reps))
        res[[i]] <- cbind(n_ind, res[[i]])
    }
    # 
    for (i in 1:length(n_ind)) {
        for (j in 1:reps) {
            # resample and save as a temporary outfile;
            resample(n_ind[i], herd_name, "fam_infile.fam", "resampled.fam")
            # Subset dataset using plink
            system2("plink", args = c("--sheep", "--file", gen_data, "--keep", "resampled.fam", 
                "--out", "filtered_data", "--make-bed", "--silent"))
            # Run KING to generate kinship files
            system2("king", args = c("-b", "filtered_data.bed", "--kinship", "--cpus", 
                "7", "--prefix", "filter_kinship", ">/dev/null"))
            # 2) Extract inbreeding values and calculate average
            system("awk '{ sum += $8; n++ } END { if (n > 0) print sum / n; }' <filter_kinship.kin> ibs_average.txt")
            # 3) Extract the chosen differentiation statistics and save to the appropriate
            # place in the appropriate dataframe
            ibs_average <- read.table(file = "ibs_average.txt")
            # within the res list;
            for (k in 1:length(res)) {
                res[[k]][i, j + 1] <- ibs_average[1, 1]
            }
        }
    }
    # 
    write.table(res, file = paste0("IBS0_", herd_name, "_", reps, "reps", "_", gen_data, 
        ".txt"), sep = " ", quote = F, row.names = F, col.names = T, append = F)
    return(as.data.frame(res))
    # 
}


## Run function

# Function arguments: reps (number of simulation replicates), n_ind (number of
# individuals drawn), herd_name (name of population as specified in file),
# gen_data (name of .ped file without extension)

resample_kinship(reps = 10000, n_ind = c(5, 10, 15, 20, 25), herd_name = "Fergus", 
    gen_data = "120_samples_table346_lowfilter")

resample_kinship(reps = 10000, n_ind = c(5, 10, 15, 20, 25), herd_name = "Glacier", 
    gen_data = "120_samples_table346_lowfilter")

resample_kinship(reps = 10000, n_ind = c(5, 10, 15, 20, 25), herd_name = "Taylor-Hilgard", 
    gen_data = "120_samples_table346_lowfilter")

resample_kinship(reps = 10000, n_ind = c(5, 10, 15, 20, 25), herd_name = "Wyoming", 
    gen_data = "120_samples_table346_lowfilter")


# Produce 30 sample results using 1 repetition Fergus
resample_kinship(reps = 1, n_ind = 30, herd_name = "Fergus", gen_data = "120_samples_table346_lowfilter")

# Glacier
resample_kinship(reps = 1, n_ind = 30, herd_name = "Glacier", gen_data = "120_samples_table346_lowfilter")

# Taylor-Hilgard
resample_kinship(reps = 1, n_ind = 30, herd_name = "Taylor-Hilgard", gen_data = "120_samples_table346_lowfilter")

# Wyoming
resample_kinship(reps = 1, n_ind = 30, herd_name = "Wyoming", gen_data = "120_samples_table346_lowfilter")
```

Example code to generate a basic boxplot of intrapopulation simulation results
------------------------------------------------------------------------------

``` r
## Evaluate results Fergus
kinship_Fergus_100_sim <- read.table("IBS0_Fergus_100reps_120_samples_table346_lowfilter.txt", 
    sep = " ", header = TRUE)

### Transpose the dataframe
kinship_Fergus_100_simb <- as.data.frame(t(kinship_Fergus_100_sim))
colnames(kinship_Fergus_100_simb) = kinship_Fergus_100_simb[1, ]  # the first row will be the header
kinship_Fergus_100_simb = kinship_Fergus_100_simb[-1, ]  # removing the first row.

boxplot(kinship_Fergus_100_simb, use.cols = TRUE)
```

Between populations simulations
===============================

Mean kinship - between populations
----------------------------------

``` r
# *****************************************************

# Mean Kinship Simulation between herds using KING

# *****************************************************

## Set working directory setwd()

## Sample individuals from .fam file in R and generate new .fam file

### Function 1
resample_interpop <- function(n_ind, herd_name1, herd_name2, fam_infile, fam_outfile) {
    # 1) Load file containing family data for all samples
    pop <- read.table(fam_infile)
    # 2) Sample from fam file a) Subset by selected herd #1 and herd #2
    popA <- pop[pop$V1 == herd_name1, ]
    popB <- pop[pop$V1 == herd_name2, ]
    # b) Sample family file for required number of individuals
    dfA <- popA[sample(nrow(popA), size = n_ind, replace = FALSE), ]
    dfB <- popB[sample(nrow(popB), size = n_ind, replace = FALSE), ]
    # c) Combine sampled results for each herd into one table
    df_all <- rbind(dfA, dfB)
    # 3) Write fam file with subsampled individuals
    write.table(df_all, file = fam_outfile, sep = " ", quote = F, row.names = F, 
        col.names = F, append = F)
}


## Function 2 gen_data should NOT include an extension, but the name of a ped file
resample_kinship_interpop <- function(reps, n_ind, herd_name1, herd_name2, gen_data) {
    # Generate family file for input data
    system2("plink", args = c("--sheep", "--file", gen_data, "--out", "fam_infile", 
        "--make-bed", "--silent"))
    # set up an output list containing several dataframe that will contain the result
    # for the differentiation statistic; name the dataframe according to the
    # differentiation statistic; the first column of each dataframe will be the
    # sample size;
    res <- list()
    for (i in c("Kinship")) {
        res[[i]] <- data.frame(matrix(vector(), length(n_ind), reps))
        res[[i]] <- cbind(n_ind, res[[i]])
    }
    # 
    for (i in 1:length(n_ind)) {
        for (j in 1:reps) {
            # resample and save as a temporary outfile;
            resample_interpop(n_ind[i], herd_name1, herd_name2, "fam_infile.fam", 
                "resampled.fam")
            # Subset dataset using plink
            system2("plink", args = c("--sheep", "--file", gen_data, "--keep", "resampled.fam", 
                "--out", "filtered_data", "--make-bed", "--silent"))
            # Run KING to generate kinship files - interpop data stored in .kin0 file
            system2("king", args = c("-b", "filtered_data.bed", "--kinship", "--prefix", 
                "filter_kinship", ">/dev/null"))
            # Extract the chosen differentiation statistics and save to the appropriate place
            # in the appropriate dataframe
            kinship_stats <- read.table(file = "filter_kinship.kin0", header = TRUE)
            
            # within the res list;
            for (k in 1:length(res)) {
                stat_value_all <- kinship_stats$Kinship  #List all kinship pair values
                stat_value <- mean(stat_value_all)  #Calculate mean kinship
                res[[k]][i, j + 1] <- stat_value
            }
        }
    }
    ## Herd names can be added here if desired res$Herd1 <- paste0(herd_name1)
    ## res$Herd2 <- paste0(herd_name2)
    write.table(res, file = paste0("Kinship_interpop_", herd_name1, "_", herd_name2, 
        "_", reps, "reps", "_", gen_data, ".txt"), sep = " ", quote = F, row.names = F, 
        col.names = T, append = F)
    return(as.data.frame(res))
    # 
}


## Run function

# Function arguments: reps (number of simulation replicates), n_ind (number of
# individuals drawn), herd_name1 (name of first population as specified in file),
# herd_name2 (name of second population as specified in file), gen_data (name of
# .ped file without extension)

resample_kinship_interpop(reps = 10000, n_ind = c(5, 10, 15, 20, 25), herd_name1 = "Fergus", 
    herd_name2 = "Glacier", gen_data = "120_samples_table346_lowfilter")

resample_kinship_interpop(reps = 10000, n_ind = c(5, 10, 15, 20, 25), herd_name1 = "Fergus", 
    herd_name2 = "Taylor-Hilgard", gen_data = "120_samples_table346_lowfilter")

resample_kinship_interpop(reps = 10000, n_ind = c(5, 10, 15, 20, 25), herd_name1 = "Fergus", 
    herd_name2 = "Wyoming", gen_data = "120_samples_table346_lowfilter")

resample_kinship_interpop(reps = 10000, n_ind = c(5, 10, 15, 20, 25), herd_name1 = "Glacier", 
    herd_name2 = "Taylor-Hilgard", gen_data = "120_samples_table346_lowfilter")

resample_kinship_interpop(reps = 10000, n_ind = c(5, 10, 15, 20, 25), herd_name1 = "Glacier", 
    herd_name2 = "Wyoming", gen_data = "120_samples_table346_lowfilter")

resample_kinship_interpop(reps = 10000, n_ind = c(5, 10, 15, 20, 25), herd_name1 = "Taylor-Hilgard", 
    herd_name2 = "Wyoming", gen_data = "120_samples_table346_lowfilter")


# Produce 30 sample results using 1 repetition
resample_kinship_interpop(reps = 1, n_ind = 30, herd_name1 = "Fergus", herd_name2 = "Glacier", 
    gen_data = "120_samples_table346_lowfilter")

resample_kinship_interpop(reps = 1, n_ind = 30, herd_name1 = "Fergus", herd_name2 = "Taylor-Hilgard", 
    gen_data = "120_samples_table346_lowfilter")

resample_kinship_interpop(reps = 1, n_ind = 30, herd_name1 = "Fergus", herd_name2 = "Wyoming", 
    gen_data = "120_samples_table346_lowfilter")

resample_kinship_interpop(reps = 1, n_ind = 30, herd_name1 = "Glacier", herd_name2 = "Taylor-Hilgard", 
    gen_data = "120_samples_table346_lowfilter")

resample_kinship_interpop(reps = 1, n_ind = 30, herd_name1 = "Glacier", herd_name2 = "Wyoming", 
    gen_data = "120_samples_table346_lowfilter")

resample_kinship_interpop(reps = 1, n_ind = 30, herd_name1 = "Taylor-Hilgard", herd_name2 = "Wyoming", 
    gen_data = "120_samples_table346_lowfilter")
```

IBS - between populations
-------------------------

``` r
# *****************************************************

# IBS Simulation between herds using KING

# *****************************************************

## Set working directory setwd()

## Sample individuals from .fam file in R and generate new .fam file

### Function 1
resample_interpop <- function(n_ind, herd_name1, herd_name2, fam_infile, fam_outfile) {
    # 1) Load file containing family data for all samples
    pop <- read.table(fam_infile)
    # 2) Sample from fam file a) Subset by selected herd #1 and herd #2
    popA <- pop[pop$V1 == herd_name1, ]
    popB <- pop[pop$V1 == herd_name2, ]
    # b) Sample family file for required number of individuals
    dfA <- popA[sample(nrow(popA), size = n_ind, replace = FALSE), ]
    dfB <- popB[sample(nrow(popB), size = n_ind, replace = FALSE), ]
    # c) Combine sampled results for each herd into one table
    df_all <- rbind(dfA, dfB)
    # 3) Write fam file with subsampled individuals
    write.table(df_all, file = fam_outfile, sep = " ", quote = F, row.names = F, 
        col.names = F, append = F)
}


## Function 2- interpop gen_data should NOT include an extension, but the name of
## a ped file
resample_kinship_interpop <- function(reps, n_ind, herd_name1, herd_name2, gen_data) {
    # Generate family file for input data
    system2("plink", args = c("--sheep", "--file", gen_data, "--out", "fam_infile", 
        "--make-bed", "--silent"))
    # set up an output list containing several dataframe that will contain the result
    # for the differentiation statistic; name the dataframe according to the
    # differentiation statistic; the first column of each dataframe will be the
    # sample size;
    res <- list()
    for (i in c("IBS")) {
        res[[i]] <- data.frame(matrix(vector(), length(n_ind), reps))
        res[[i]] <- cbind(n_ind, res[[i]])
    }
    # 
    for (i in 1:length(n_ind)) {
        for (j in 1:reps) {
            # resample and save as a temporary outfile;
            resample_interpop(n_ind[i], herd_name1, herd_name2, "fam_infile.fam", 
                "resampled.fam")
            # Subset dataset using plink
            system2("plink", args = c("--sheep", "--file", gen_data, "--keep", "resampled.fam", 
                "--out", "filtered_data", "--make-bed", "--silent"))
            # Run KING to generate kinship files - interpop data stored in .kin0 file
            system2("king", args = c("-b", "filtered_data.bed", "--kinship", "--prefix", 
                "filter_kinship", ">/dev/null"))
            # 2) Extract IBS values and calculate average
            system("awk '{ sum += $7; n++ } END { if (n > 0) print sum / n; }' <filter_kinship.kin0> kinship_average.txt")
            # 3) Extract the chosen differentiation statistics and save to the appropriate
            # place in the appropriate dataframe
            kinship_average <- read.table(file = "kinship_average.txt")
            # Summarize values
            for (k in 1:length(res)) {
                res[[k]][i, j + 1] <- kinship_average[1, 1]
            }
        }
    }
    write.table(res, file = paste0("IBS0_interpop_", herd_name1, "_", herd_name2, 
        "_", reps, "reps", "_", gen_data, ".txt"), sep = " ", quote = F, row.names = F, 
        col.names = T, append = F)
    return(as.data.frame(res))
    # 
}


## Run function

# Function arguments: reps (number of simulation replicates), n_ind (number of
# individuals drawn), herd_name1 (name of first population as specified in file),
# herd_name2 (name of second population as specified in file), gen_data (name of
# .ped file without extension)

# Function call
resample_kinship_interpop(reps = 10000, n_ind = c(5, 10, 15, 20, 25), herd_name1 = "Fergus", 
    herd_name2 = "Glacier", gen_data = "120_samples_table346_lowfilter")

resample_kinship_interpop(reps = 10000, n_ind = c(5, 10, 15, 20, 25), herd_name1 = "Fergus", 
    herd_name2 = "Taylor-Hilgard", gen_data = "120_samples_table346_lowfilter")

resample_kinship_interpop(reps = 10000, n_ind = c(5, 10, 15, 20, 25), herd_name1 = "Fergus", 
    herd_name2 = "Wyoming", gen_data = "120_samples_table346_lowfilter")

resample_kinship_interpop(reps = 10000, n_ind = c(5, 10, 15, 20, 25), herd_name1 = "Glacier", 
    herd_name2 = "Taylor-Hilgard", gen_data = "120_samples_table346_lowfilter")

resample_kinship_interpop(reps = 10000, n_ind = c(5, 10, 15, 20, 25), herd_name1 = "Glacier", 
    herd_name2 = "Wyoming", gen_data = "120_samples_table346_lowfilter")

resample_kinship_interpop(reps = 10000, n_ind = c(5, 10, 15, 20, 25), herd_name1 = "Taylor-Hilgard", 
    herd_name2 = "Wyoming", gen_data = "120_samples_table346_lowfilter")


# Produce 30 sample results using 1 repetition
resample_kinship_interpop(reps = 1, n_ind = 30, herd_name1 = "Fergus", herd_name2 = "Glacier", 
    gen_data = "120_samples_table346_lowfilter")

resample_kinship_interpop(reps = 1, n_ind = 30, herd_name1 = "Fergus", herd_name2 = "Taylor-Hilgard", 
    gen_data = "120_samples_table346_lowfilter")

resample_kinship_interpop(reps = 1, n_ind = 30, herd_name1 = "Fergus", herd_name2 = "Wyoming", 
    gen_data = "120_samples_table346_lowfilter")

resample_kinship_interpop(reps = 1, n_ind = 30, herd_name1 = "Glacier", herd_name2 = "Taylor-Hilgard", 
    gen_data = "120_samples_table346_lowfilter")

resample_kinship_interpop(reps = 1, n_ind = 30, herd_name1 = "Glacier", herd_name2 = "Wyoming", 
    gen_data = "120_samples_table346_lowfilter")

resample_kinship_interpop(reps = 1, n_ind = 30, herd_name1 = "Taylor-Hilgard", herd_name2 = "Wyoming", 
    gen_data = "120_samples_table346_lowfilter")
```

Fst - between populations
-------------------------

``` r
# ***********************************************************

## Fst Simulation between populations using PLINK BETA 1.90

# ***********************************************************

## Set working directory Location on Linux laptop setwd()

## Ensure that PLINK BETA 1.90 is the active PLINK program available

## Sample individuals from .fam file in R and generate new .fam file

### Function 1-working
resample_interpop <- function(n_ind, herd_name1, herd_name2, fam_infile, fam_outfile) {
    # 1) Load file containing family data for all samples
    pop <- read.table(fam_infile)
    # 2) Sample from fam file a) Subset by selected herd #1 and herd #2
    popA <- pop[pop$V1 == herd_name1, ]
    popB <- pop[pop$V1 == herd_name2, ]
    # b) Sample family file for required number of individuals
    dfA <- popA[sample(nrow(popA), size = n_ind, replace = FALSE), ]
    dfB <- popB[sample(nrow(popB), size = n_ind, replace = FALSE), ]
    # c) Combine sampled results for each herd into one table
    df_all <- rbind(dfA, dfB)
    # 3) Write fam file with subsampled individuals
    write.table(df_all, file = fam_outfile, sep = " ", quote = F, row.names = F, 
        col.names = F, append = F)
}


## Function 2- interpop working gen_data should NOT include an extension, but the
## name of a ped file
resample_Fst_interpop <- function(reps, n_ind, herd_name1, herd_name2, gen_data) {
    # Generate family file for input data
    system2("plink", args = c("--sheep", "--file", gen_data, "--out", "fam_infile", 
        "--make-bed", "--silent"))
    # set up an output list containing several dataframe that will contain the result
    # for the differentiation statistic; name the dataframe according to the
    # differentiation statistic; the first column of each dataframe will be the
    # sample size;
    res <- list()
    for (i in c("Fst")) {
        res[[i]] <- data.frame(matrix(vector(), length(n_ind), reps))
        res[[i]] <- cbind(n_ind, res[[i]])
    }
    res1 <- res
    res2 <- res
    res3 <- res
    # 
    for (i in 1:length(n_ind)) {
        for (j in 1:reps) {
            # resample and save as a temporary outfile;
            resample_interpop(n_ind[i], herd_name1, herd_name2, "fam_infile.fam", 
                "resampled.fam")
            # 1) Subset dataset using plink
            system2("plink", args = c("--sheep", "--file", gen_data, "--keep", "resampled.fam", 
                "--out", "filtered_data", "--make-bed", "--silent"))
            # 2) Generate .phe file from family file to define populations
            system("awk '{ print $1 \"\t\" $2 \"\t\" $1}' <resampled.fam> resampled.phe")
            # Filter SNPs for quality control 1) MAF < 0.0001 removed for all samples
            # included
            system2("plink", args = c("--sheep", "--bfile", "filtered_data", "--maf", 
                "0.0001", "--within", "resampled.phe", "--allow-no-sex", "--out", 
                "filtered_data_1", "--make-bed", "--silent"))
            # Produce file summarizing count of the number of alleles- TO-DO: check if this
            # should be run within family
            system2("plink", args = c("--sheep", "--bfile", "filtered_data_1", "--freq", 
                "counts", "--allow-no-sex", "--out", "filtered_data_1_summary", "--silent"))
            system("awk 'END { print NR - 1 }' <filtered_data_1_summary.frq.counts> filtered_data_1_count.txt")
            # 2) Deviate from HWe p<0.00001
            system2("plink", args = c("--sheep", "--bfile", "filtered_data_1", "--hwe", 
                "0.00001", "midp", "--allow-no-sex", "--within", "resampled.phe", 
                "--out", "filtered_data_2", "--make-bed", "--silent"))
            # Produce file summarizing count of the number of alleles
            system2("plink", args = c("--sheep", "--bfile", "filtered_data_2", "--freq", 
                "counts", "--allow-no-sex", "--out", "filtered_data_2_summary", "--silent"))
            system("awk 'END { print NR - 1 }' <filtered_data_2_summary.frq.counts> filtered_data_2_count.txt")
            # 3) LD pruning- window size 100, window increment 25, r2 threshold 0.99
            system2("plink", args = c("--sheep", "--bfile", "filtered_data_2", "--indep-pairwise", 
                "100", "25", "0.99", "-allow-no-sex", "--out", "filtered_data_3", 
                "--silent"))
            # Produce file excluding LD alleles
            system2("plink", args = c("--sheep", "--bfile", "filtered_data_2", "--extract", 
                "filtered_data_3.prune.in", "--out", "filtered_data_3b", "--make-bed", 
                "--silent"))
            # Produce file summarizing count of remaining alleles
            system2("plink", args = c("--sheep", "--bfile", "filtered_data_3b", "--freq", 
                "counts", "--allow-no-sex", "--out", "filtered_data_3b_summary", 
                "--silent"))
            system("awk 'END { print NR - 1 }' <filtered_data_3b_summary.frq.counts> filtered_data_3_count.txt")
            # 
            
            ## 3) Use Plink 1.90 version to calculate Fst per marker
            system2("plink", args = c("--sheep", "--bfile", "filtered_data_3b", "--fst", 
                "--within", "resampled.phe", "--out", "filter_Fst", "--allow-no-sex", 
                "--silent"))
            # 4) Calculate average Fst by extracting column from .fst file- also recorded in
            # log file
            system("awk '{if($5!=\"nan\"){ sum += $5; n++ }} END { if (n > 0) print sum / n; }' <filter_Fst.fst> fst_average.txt")
            # 5) Extract the chosen differentiation statistics and save to the appropriate
            # place in the appropriate dataframe
            fst_average <- read.table(file = "fst_average.txt")
            filtered_data_1_count <- read.table(file = "filtered_data_1_count.txt")
            filtered_data_2_count <- read.table(file = "filtered_data_2_count.txt")
            filtered_data_3_count <- read.table(file = "filtered_data_3_count.txt")
            # within the res list;
            for (k in 1:length(res)) {
                res[[k]][i, j + 1] <- fst_average[1, 1]
                res1[[k]][i, j + 1] <- filtered_data_1_count[1, 1]
                res2[[k]][i, j + 1] <- filtered_data_2_count[1, 1]
                res3[[k]][i, j + 1] <- filtered_data_3_count[1, 1]
            }
        }
    }
    write.table(res, file = paste0("Fst_interpop_", herd_name1, "_", herd_name2, 
        "_", reps, "reps", "_", gen_data, ".txt"), sep = " ", quote = F, row.names = F, 
        col.names = T, append = F)
    write.table(res1, file = paste0("Fst_interpop_", herd_name1, "_", herd_name2, 
        "_", reps, "reps", "_", gen_data, "_alleles_filter1", ".txt"), sep = " ", 
        quote = F, row.names = F, col.names = T, append = F)
    write.table(res2, file = paste0("Fst_interpop_", herd_name1, "_", herd_name2, 
        "_", reps, "reps", "_", gen_data, "_alleles_filter2", ".txt"), sep = " ", 
        quote = F, row.names = F, col.names = T, append = F)
    write.table(res3, file = paste0("Fst_interpop_", herd_name1, "_", herd_name2, 
        "_", reps, "reps", "_", gen_data, "_alleles_filter3", ".txt"), sep = " ", 
        quote = F, row.names = F, col.names = T, append = F)
    # 
    return(as.data.frame(res))
    # 
}


## Run function

# Function arguments: reps (number of simulation replicates), n_ind (number of
# individuals drawn), herd_name1 (name of first population as specified in file),
# herd_name2 (name of second population as specified in file), gen_data (name of
# .ped file without extension)

resample_Fst_interpop(reps = 10000, n_ind = c(5, 10, 15, 20, 25, 30), herd_name1 = "Fergus", 
    herd_name2 = "Glacier", gen_data = "120_samples_table346_lowfilter")

resample_Fst_interpop(reps = 10000, n_ind = c(5, 10, 15, 20, 25), herd_name1 = "Fergus", 
    herd_name2 = "Taylor-Hilgard", gen_data = "120_samples_table346_lowfilter")

resample_Fst_interpop(reps = 10000, n_ind = c(5, 10, 15, 20, 25), herd_name1 = "Fergus", 
    herd_name2 = "Wyoming", gen_data = "120_samples_table346_lowfilter")

resample_Fst_interpop(reps = 10000, n_ind = c(5, 10, 15, 20, 25), herd_name1 = "Glacier", 
    herd_name2 = "Taylor-Hilgard", gen_data = "120_samples_table346_lowfilter")

resample_Fst_interpop(reps = 10000, n_ind = c(5, 10, 15, 20, 25), herd_name1 = "Glacier", 
    herd_name2 = "Wyoming", gen_data = "120_samples_table346_lowfilter")

resample_Fst_interpop(reps = 10000, n_ind = c(5, 10, 15, 20, 25), herd_name1 = "Taylor-Hilgard", 
    herd_name2 = "Wyoming", gen_data = "120_samples_table346_lowfilter")


# Produce 30 sample results using 1 repetition
resample_Fst_interpop(reps = 1, n_ind = 30, herd_name1 = "Fergus", herd_name2 = "Glacier", 
    gen_data = "120_samples_table346_lowfilter")

resample_Fst_interpop(reps = 1, n_ind = 30, herd_name1 = "Fergus", herd_name2 = "Taylor-Hilgard", 
    gen_data = "120_samples_table346_lowfilter")

resample_Fst_interpop(reps = 1, n_ind = 30, herd_name1 = "Fergus", herd_name2 = "Wyoming", 
    gen_data = "120_samples_table346_lowfilter")

resample_Fst_interpop(reps = 1, n_ind = 30, herd_name1 = "Glacier", herd_name2 = "Taylor-Hilgard", 
    gen_data = "120_samples_table346_lowfilter")

resample_Fst_interpop(reps = 1, n_ind = 30, herd_name1 = "Glacier", herd_name2 = "Wyoming", 
    gen_data = "120_samples_table346_lowfilter")

resample_Fst_interpop(reps = 1, n_ind = 30, herd_name1 = "Taylor-Hilgard", herd_name2 = "Wyoming", 
    gen_data = "120_samples_table346_lowfilter")
```

Report session information
==========================

``` r
sessionInfo()
```

    ## R version 3.2.3 (2015-12-10)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 16.04.2 LTS
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] backports_1.1.2 magrittr_1.5    rprojroot_1.3-2 formatR_1.5    
    ##  [5] tools_3.2.3     htmltools_0.3.6 yaml_2.1.16     Rcpp_0.12.15   
    ##  [9] stringi_1.1.6   rmarkdown_1.8   knitr_1.19      stringr_1.2.0  
    ## [13] digest_0.6.15   evaluate_0.10.1

Citations
=========

Manichaikul, A., Mychaleckyj, J. C., Rich, S. S., Daly, K., Sale, M., & Chen, W.-M. (2010). Robust relationship inference in genome-wide association studies. Bioinformatics, 26(22), 2867–2873. <https://doi.org/10.1093/bioinformatics/btq559>

Nazareno, A. G., Bemmels, J. B., Dick, C. W., & Lohmann, L. G. (2017). Minimum sample sizes for population genomics: An empirical study from an Amazonian plant species. Molecular Ecology Resources, n/a-n/a. <https://doi.org/10.1111/1755-0998.12654>

Purcell, S., Neale, B., Todd-Brown, K., Thomas, L., Ferreira, M. A. R., Bender, D., … Sham, P. C. (2007). PLINK: A Tool Set for Whole-Genome Association and Population-Based Linkage Analyses. The American Journal of Human Genetics, 81(3), 559–575. <https://doi.org/10.1086/519795>
