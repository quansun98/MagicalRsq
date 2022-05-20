# MagicalRsq package documentation

## Introduction

MagicalRsq is a novel quality estimates of genotype imputation using XGBoost method. By default, it starts from the standard Rsq and estimated MAF given by minimac imputation software, integrates population genetics statistics and population-specific allele frequencies, and calibrates the original Rsq. You can also train MagicalRsq models yourself with whatever variant-level features you like and any amounts of variants you want. 

## Data

[Population genetics (S/HIC) features](ftp://yunlianon:anon@rc-ns-ftp.its.unc.edu/MagicalRsq/SHIC/)
[TOP-LD ancestry specific allele frequencies](ftp://yunlianon:anon@rc-ns-ftp.its.unc.edu/MagicalRsq/AF/)

## Installation

1. First, install the `devtools` package if you haven't already. 

		install.packages("devtools")

2. Load the devtools package and install through Github.

		library(devtools)
		install_github("Quanbaby/MagicalRsq")


## Quick start

Below are the steps to train MagicalRsq models in training cohort and apply to target cohort. Note that it is not necessary to have the same variants in training and target dataset.

### Integrate variant-level features

Following the steps below if you want to integrate S/HIC features and population specific allele frequency features. You can skip this if your data already contain all the variant-level features you want to leverage. The below demonstration will use one chromosome (chr22) as an example.
Note that:

* All the positions are in hg38. You need to liftover if you are using referene panels in hg19, e.g. 1000G, HRC, etc..
* If you choose to integrate S/HIC or population AF features, you need to run `data_integrate` function separately for each chromosome.
* Training data and testing data must have the same features. No need to keep feature orders the same, since you can specifiy the feature column numbers when training and calculating MagicalRsq. But we highly recommend you keep the same feature order in training and target data.

1. Start from minimac output `.info.gz` file:

		zcat data/chr22.info.gz | sed '1d' | cut -f 1,5,7 | sed 's/:/\t/g' | sed 's/^chr//g' | sed '1iCHR POS REF ALT MAF Rsq' | sed 's/ /\t/g' | bgzip > chr22.imp.sumstat.txt.gz

You should also integrate true R2 information in this file and used for training models or evaluations.

2. Download S/HIC and allele frequency features.

		wget -r ftp://yunlianon:anon@rc-ns-ftp.its.unc.edu/MagicalRsq/SHIC/
		wget -r ftp://yunlianon:anon@rc-ns-ftp.its.unc.edu/MagicalRsq/AF/

2. Run `data_integrate` function in R

		library(MagicalRsq)
		chr22_for_MRsq = data_integrate("chr22.imp.sumstat.txt.gz", chr = 22, pos_col = 2, SHIC_dir = "../SHIC/", AF_dir = "../AF/", outfile = "test_chr22.txt.gz")

You will get integrated data for chr22 written in the file `test_chr22.txt.gz`.


### Train MagicalRsq models

In order to train MagicalRsq models, you need the true R2 (true imputation quality) for all/subset of variants among subset/all of individuals. The procedure will be demonstrated with the example data `toy_chr22_50k_integrated.txt.gz` in the `data/` folder.

1. Similarly, you need to integrate imputation summary statistics with other variant-level features. Here we use both S/HIC and population specific AF as included features. You can choose to include them or not if you have your own features.

		integrated = data_integrate("data/toy_chr22_50k.txt.gz", chr = 22, pos_col = 2, 
			SHIC_dir = "../SHIC/", AF_dir = "../AF/", outfile = "toy_chr22_50k_integrated.txt.gz") 

2. Train MagicalRsq models with integrated data. You can choose how many variants you want to include in the training models. We recommend using 10k - 1,000k variants would be enough. Too many variants will increase computational burden, making the function take longer to finish.

		toy_model = train_MagicalRsq(file = "toy_chr22_50k_integrated.txt.gz", outfile = "toy_model", nvar_use = 100000) 

Note that you can use whole genome data (separated by chromosomes) to train a whole-genome model. Input files could be a vector of file names.

	toy2_model = train_MagicalRsq(file = paste0("toy2_chr",1:22,".txt.gz"), outfile = "toy2", nvar_use = 500000)

### Calculate MagicalRsq values

After getting the integrated data, you can calculate MagicalRsq values for each variant in the integrated data file.

        chr22_MagicalRsq  = calc_MagicalRsq(file = "test_chr22.txt.gz", model = "toy_model",
                FeatureCols = 6:84, keptCols = 1:7, outfile = "test_MagicalRsq.txt.gz")

You will get variant-level MagicalRsq values, together with the original Rsq and estimated MAF written in the file `test_MagicalRsq.txt.gz`


### Evaluate Rsq/MagicalRsq performance (true R2 required)

After calculation of MagicalRsq, you can evaluate the performance of Rsq and/or MagicalRsq compared to true R2. You may want to specify the column numbers of trueR2, Rsq and MagicalRsq. You can only choose to disregard extremely rare variants for evaluation purpose by setting MAF or (MAC and N samples), and specifying the column number of MAF.

	chr22_eval = eval_MagicalRsq(file = "test_MagicalRsq.txt.gz", MAF = 0.001, trueR2Col = 8, MAFCol = 5, RsqCol = 6, MagicalRsqCol = 7)
 
All 22 chromosomes can also be evaluated together:

	all_eval = eval_MagicalRsq(file = paste0("test2_chr",1:22,".txt.gz"))



