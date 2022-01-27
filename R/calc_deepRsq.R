#' Calculate deepRsq values using pre-trained deepRsq models
#'
#' This function is used to calculate deepRsq values given pre-trained deepRsq models from `train_deepRsq` function.
#' Variant level features in the input data should be exactly the same as the input data for `train_deepRsq` function.
#'
#' @param file Required. A vector of strings with each element being one file name. Could be a single file or multiple files (if separated by chromosomes).
#' @param header Logical. The parameter is passing to the `fread` function. TRUE by default. Strongly recommended containing a header in your data file.
#' @param model Required. The pre-trained deepRsq model used as references for calculating deepRsq values. Output from the `train_deepRsq` function.
#' @param FeatureCols Vector of column numbers of all the features used for model training. By default, it will be 5:85 including 81 features covering
#' population genetics, allele frequencies across populations, etc., provided online. More user-defined features are also feasible. 
#' Note that the default values are different from `train_deepRsq` function since we assume true R2 information is not available in the input data here.
#' @param keptCols Vector of column numbers that needs to be carried over into the output files, e.g. variants ID, chromosome, positions, alleles, estimated MAF. etc.
#' Default values 1:6, including Chr, Pos, Ref, Alt, MAF and Rsq.
#' @param outfile Name of the output file containing the calculated deepRsq values for each variant. Strongly recommended to specify the output file name to save the data.
#' @param seed Random seed. Default value 123.
#'
#' @return Dataframe of the calculated deepRsq values, along with the information carried over from the input data (keptCols of the input data).
#'
#' @examples
#'
#' MESA_deepRsq = calc_deepRsq(file = "MESA_1000G_allchr.txt.gz", 
#' model = JHS_model, FeatureCols = 6:86, keptCols = 1:7, 
#' outfile = "MESA_deepRsq_JHS_model.txt.gz")
#'
#' ## If data are separated by chromosomes 
#' WHI_deepRsq = calc_deepRsq(file = paste0("WHI_1000G_chr",1:22,".txt.gz"), 
#' model = JHS_model,  outfile = "WHI_deepRsq_JHS_model.txt.gz")
#'
#' @importFrom data.table fread fwrite
#'
#' @export

calc_deepRsq = function(file, header = T, model, FeatureCols = 5:85, keptCols = 1:6, 
  outfile = NULL, seed = 123){

set.seed(seed)

if(is.null(model)){
stop("No deepRsq model input")
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
tmp = fread(file, header = header)
dat = rbind(dat,tmp)
}
}

dat = as.data.frame(dat)

output = dat[,keptCols]
output$deepRsq = predict(bst, as.matrix(dat[,FeatureCols]))

if(!is.null(outfile)){
message(paste("Saving data to", outfile))
fwrite(output, file = outfile, quote = F, sep = "\t", col.names = T, row.names = F)
}

return(output)

}


