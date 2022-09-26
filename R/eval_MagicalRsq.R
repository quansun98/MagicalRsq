#' Evaluate MagicalRsq performance (compared to true R2)
#' 
#' This function is used to evaluate MagicalRsq values calculated by `calc_MagicalRsq` in the presence of true R2.
#' 
#' @param file Required. A vector of strings with each element being one file name. Could be a single file or multiple files (if separated by chromosomes).
#' @param header Logical. The parameter is passing to the `fread` function. TRUE by default. Strongly recommended containing a header in your data file.
#' @param minMAC Disgard variants with extremely low minor allele counts for the evaluation purpose. Must be used together with `N` and should not be used together with `MAF`.
#' @param N Total sample size in terms of individuals in the dataset. Used to calculate MAF as the variants inclusion criteria. Must be used together with `MAC` should not be used together with `MAF`. 
#' @param minMAF Disgard variants with extremely low minor allele frequencies for the evaluation purpose. If specified some cutoff c, it will only evaluate variants with `MAF > c`. Should not be used together with `minMAC` or `N`. 
#' Note that if `minMAC`,`N`, and `minMAF` are all provided, `minMAF` will be used as the primary threshold, and the other two parameters won't take any effects.
#' @param trueR2Col Column number of true R2. Default value 5.
#' @param MAFCol Column number of estimated minor allele frequency. Default value 6.
#' @param RsqCol Column number of Rsq. Default value 7.
#' @param MagicalRsqCol Column number of MagicalRsq. Default value 8.
#'
#' @return A named one-row dataframe with three metrics: R2 (squared Pearson correlation with true R2), RMSE (root mean square error) and MAE (mean absolute error) for both Rsq and MagicalRsq. Six columns in total.
#'
#' @examples
#'
#' test_eval = eval_MagicalRsq(file = "test_MagicalRsq.txt.gz", minMAC = 5, N = 2000)
#'
#' test2_eval = eval_MagicalRsq(file = paste0("test2_chr",1:22,".txt.gz"), 
#' minMAF = 0.001, trueR2Col = 8, MAFCol = 5, RsqCol = 6, MagicalRsqCol = 7)
#'
#' @export

eval_MagicalRsq = function(file, header = T, minMAC = NULL, N = NULL, minMAF = NULL, trueR2Col = 5, MAFCol = 6, RsqCol = 7, MagicalRsqCol = 8){

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

# Filter out low MAF variants if specified
if((!is.null(minMAC) && !is.null(N)) || (!is.null(minMAF))){
message("Evaluation not including extremely rare variants")
if(is.null(minMAF)){
minMAF = as.numeric(minMAC)/2/as.numeric(N)
}

message(paste("Including variants with MAF > ", minMAF))
dat = dat[dat[,MAFCol] > minMAF,]
}

# evaluate MagicalRsq and Rsq performance
final = as.data.frame(matrix(0, nrow=1, ncol = 6))
colnames(final) = c("R2_Rsq","R2_MagicalRsq","RMSE_Rsq","RMSE_MagicalRsq","MAE_Rsq","MAE_MagicalRsq")

final$R2_Rsq = cor(dat[,trueR2Col], dat[,RsqCol])^2
final$R2_MagicalRsq = cor(dat[,trueR2Col], dat[,MagicalRsqCol])^2
final$RMSE_Rsq = sqrt(mean((dat[,trueR2Col] - dat[,RsqCol])^2))
final$RMSE_MagicalRsq = sqrt(mean((dat[,trueR2Col] - dat[,MagicalRsqCol])^2))
final$MAE_Rsq = mean(abs(dat[,trueR2Col] - dat[,RsqCol]))
final$MAE_MagicalRsq = mean(abs(dat[,trueR2Col] - dat[,MagicalRsqCol]))

return(final)

}


