# DeepRsq package documentation

## Introduction

DeepRsq is a novel quality estimates of genotype imputation using XGBoost method. By default, it starts from the standard Rsq and estimated MAF given by minimac imputation software, integrates population genetics statistics and population-specific allele frequencies, and calibrates the original Rsq. You can also train deepRsq models yourself with whatever variant-level features you like and any amounts of variants you want. 

## Data

* Population genetics (S/HIC) features can be downloaded at <ftp://yunlianon:anon@rc-ns-ftp.its.unc.edu/deepRsq/SHIC/>
* Population specific allele frequencies can be downloaded at <ftp://yunlianon:anon@rc-ns-ftp.its.unc.edu/deepRsq/AF/>

## Installation

1. First, install the `devtools` package if you haven't already. 

		install.packages("devtools")

2. Load the devtools package and install through Github.

		library(devtools)
		install_github("Quanbaby/deepRsq")


## Quick start

If you don't have true R2 information, and just want to calculated deepRsq values based on pre-calculated models. Here are the steps.

### Integrate variant-level features

Following the steps below if you want to integrate S/HIC features and population specific allele frequency features. The below demonstration will use one chromosome (chr22) as an example.
Note that:

* All the positions are in hg38. You need to liftover if you are using referene panels in hg19, e.g. 1000G, HRC, etc..
* If you choose to integrate S/HIC or population AF features, you need to run `data_integrate` function separately for each chromosome.

1. Start from minimac output `.info.gz` file:

		zcat data/chr22.info.gz | sed '1d' | cut -f 1,5,7 | sed 's/:/\t/g' | sed 's/^chr//g' | sed '1iCHR POS REF ALT MAF Rsq' | sed 's/ /\t/g' | bgzip > chr22.imp.sumstat.txt.gz

You can also integrate true R2 information in this file and used for training models or evaluations.

2. Download S/HIC and allele frequency features.

		wget -r ftp://yunlianon:anon@rc-ns-ftp.its.unc.edu/deepRsq/SHIC/
		wget -r ftp://yunlianon:anon@rc-ns-ftp.its.unc.edu/deepRsq/AF/

2. Run `data_integrate` function in R

		library(deepRsq)
		chr22_for_deepRsq = data_integrate("chr22.imp.sumstat.txt.gz", chr = 22, pos_col = 2, SHIC_dir = "../SHIC/", AF_dir = "../AF/", outfile = "chr22_integrated.txt.gz")

You will get integrated data for chr22 written in the file `chr22_integrated.txt.gz`.


### Train deepRsq models yourself

You can also train deepRsq models yourself, if you have the true R2 (true imputation quality) for all/subset of variants. The procedure will be demonstrated with the example data `BioMe_EUR_Rsq_trueR2_chr22_50k.txt.gz` in the `data/` folder.

1. Similarly, you need to integrate imputation summary statistics with other variant-level features. Here we use both S/HIC and population specific AF as included features. You can choose to include them or not.

		integrated = data_integrate("data/BioMe_EUR_Rsq_trueR2_chr22_50k.txt.gz", chr = 22, pos_col = 2, 
			SHIC_dir = "../SHIC/", AF_dir = "../AF/", outfile = "BioMe_EUR_Rsq_trueR2_chr22_integrated_50k.txt.gz") 

2. Train deepRsq models with integrated data. You can choose how many variants you want to include in the training models. We recommend using 10k - 1,000k variants would be enough. Too many variants will increase computational burden, making the function take longer to finish.

		toy_model = train_deepRsq(file = "BioMe_EUR_Rsq_trueR2_chr22_integrated_50k.txt.gz", outfile = "toy_model", nvar_use = 100000) 

Note that you can use whole genome data (separated by chromosomes) to train a whole-genome model. Input files could be a vector of file names.

		JHS_model = train_deepRsq(file = paste0("JHS_1000G_chr",1:22,".txt.gz"), outfile = "JHS_trained", nvar_use = 500000)

### Calculate deepRsq values

After getting the integrated data, you can calculate deepRsq values for each variant in the integrated data file.

                chr22_deepRsq  = calc_deepRsq(file = "chr22_integrated.txt.gz", model = "model/BioMe_EUR_uncommoon_200k",
                        FeatureCols = 5:85, keptCols = 1:6, outfile = "chr22_deepRsq.txt.gz")

You will get variant-level deepRsq values, together with the original Rsq and estimated MAF written in the file `chr22_deepRsq.txt.gz`


### Evaluate Rsq/deepRsq performance (true R2 required)

After calculation of deepRsq, you can evaluate the performance of Rsq and/or deepRsq compared to true R2. You may want to specify the column numbers of trueR2, Rsq and deepRsq. You can only choose to disregard extremely rare variants for evaluation purpose by setting MAF or (MAC and N samples), and specifying the column number of MAF.

		chr22_eval = eval_deepRsq(file = "MESA_deepRsq_JHS_model.txt.gz", MAF = 0.001, trueR2Col = 8, MAFCol = 5, RsqCol = 6, deepRsqCol = 7)

All 22 chromosomes can also be evaluated together:

		all_eval = eval_deepRsq(file = paste0("BioMe_deepRsq_MESA_model_chr",1:22,".txt.gz"))



