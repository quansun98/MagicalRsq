#' Calculate MagicalRsq values using pre-trained MagicalRsq models
#'
#' This function is used to calculate MagicalRsq values given pre-trained MagicalRsq models from `train_MagicalRsq` function.
#' Variant level features in the input data should be exactly the same as the input data for `train_MagicalRsq` function.
#'
#' @param file Required. A vector of strings with each element being one file name. Could be a single file or multiple files (if separated by chromosomes).
#' @param header Logical. The parameter is passing to the `fread` function. TRUE by default. Strongly recommended containing a header in your data file.
#' @param model Required. The pre-trained MagicalRsq model used as references for calculating MagicalRsq values. Output from the `train_MagicalRsq` function.
#' @param FeatureCols Vector of column numbers of all the features used for model training. By default, it will be 5:83 including 79 features covering
#' population genetics, allele frequencies across populations, etc., provided online. More user-defined features are also feasible. 
#' Note that the default values are different from `train_MagicalRsq` function since we assume true R2 information is not available in the input data here.
#' @param keptCols Vector of column numbers that needs to be carried over into the output files, e.g. variants ID, chromosome, positions, alleles, estimated MAF. etc.
#' Default values 1:6, including Chr, Pos, Ref, Alt, MAF and Rsq.
#' @param MAF_cate Minor allee frequency category. Must be chosen from c("common","lowfreq","rare). If not specified, it will use all the variants to train models.
#' @param MAFCol Column number of minor allele frequency. Optional. Must be specified if `MAF_cate` is used.
#' @param outfile Name of the output file containing the calculated MagicalRsq values for each variant. Strongly recommended to specify the output file name to save the data.
#' @param seed Random seed. Default value 123.
#'
#' @return Dataframe of the calculated MagicalRsq values, along with the information carried over from the input data (keptCols of the input data).
#'
#' @examples
#'
#' toy_MagicalRsq = calc_MagicalRsq(file = "test_chr22.txt.gz", 
#' model = toy_model, FeatureCols = 6:84, keptCols = 1:7, 
#' outfile = "test_MagicalRsq.txt.gz")
#'
#' ## If data are separated by chromosomes 
#' toy2_MagicalRsq = calc_MagicalRsq(file = paste0("toy2_chr",1:22,".txt.gz"), 
#' model = toy2,  outfile = "toy2_MagicalRsq.txt.gz")
#'
#' @importFrom data.table fread fwrite
#'
#' @export

calc_MagicalRsq = function(file, header = T, model, FeatureCols = 5:83, keptCols = 1:6, 
  MAF_cate = NULL, MAFCol = NULL, outfile = NULL, seed = 123){

set.seed(seed)

if(is.null(model)){
stop("No MagicalRsq model input")
}

bst = xgb.load(model)

n_file = length(file)

if(n_file == 1){
message("1 file detected")
dat = fread(file, header = header)
} else if(n_file == 0){
stop("No file detected. Please check your input files")
} else if(n_file >= 2){
message(paste(n_file,"files detected"))
dat = fread(file[1], header = header)
# multiple files
for(i in 2:n_file){
tmp = fread(file[i], header = header)
dat = rbind(dat,tmp)
}
}

dat = as.data.frame(dat)

# MAF filter

if(!is.null(MAF_cate)){
if(is.null(MAFCol)){
stop("MAF category specified, but MAF column doesn't specified")
}else if(MAF_cate == "common"){
dat = dat[which(dat[,MAFCol] >= 0.05), ]
}else if(MAF_cate == "lowfreq"){
dat = dat[which((dat[,MAFCol] > 0.005) & (dat[,MAFCol] < 0.05)), ]
}else if(MAF_cate == "rare"){
dat = dat[which(dat[,MAFCol] <= 0.005),]
}else{
stop("MAF category must have values chosen from 'common','lowfreq' and 'rare'")
}
}



output = dat[,keptCols]
output$MagicalRsq = predict(bst, as.matrix(dat[,FeatureCols]))

if(!is.null(outfile)){
message(paste("Saving data to", outfile))
fwrite(output, file = outfile, quote = F, sep = "\t", col.names = T, row.names = F)
}

return(output)

}


