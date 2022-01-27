#' Train deepRsq model
#'
#' This function is used to train deepRsq models with data files containing all the features used.
#' 
#' @param file Required. A vector of strings with each element being one file name. Could be a single file or multiple files (if separated by chromosomes).
#' @param header Logical. The parameter is passing to the `fread` function. TRUE by default. Strongly recommended containing a header in your data file.
#' @param outfile Name of the deepRsq model output. Strings. Strongly recommended to specify the output file name to save the model.
#' @param nvar_use Number of variants to use for the model training and validation purpose. If not specified, all variants in the data file will be used. 
#' Note that the model training will take much longer time if using full set of variants. Recommended values ranging from 10,000 to 1,000,000. 
#' @param p_train Proportion of variants used in training part. The rest of the variants would be used for validation purpose. Default value is 70%. 
#' Note that the total number of variants for validation is `nvar_use`*(1 - `p_train`).
#' @param trueR2Col Column number of true R2. Default value 5.
#' @param FeatureCols Vector of column numbers of all the features used for model training. By default, it will be 6:86 including 81 features covering 
#' population genetics, allele frequencies across populations, etc., provided online. More user-defined features are also feasible. 
#' @param seed Random seed. Default value 123.
#' @param verbose Parameter passing to the `xgb.train` function: If 0, xgboost will stay silent. If 1, it will print information about performance. If 2, some additional information will be printed out. 
#' @param nrounds Parameter passing to the `xgb.train` function: max number of boosting iterations. Default value 3000.
#' @param early_stopping_rounds Parameter passing to the `xgb.train` function: training with a validation set will stop if the performance doesn't improve for `k` rounds. Default value 50.
#' @param print_every_n Parameter passing to the `xgb.train` function: Print each n-th iteration evaluation messages when `verbose > 0`. Default is 20.
#' 
#' @return The trained deepRsq XGBoost model.
#'
#' @examples
#' 
#' MESA_model = train_deepRsq(file = "MESA_1000G_allchr.txt.gz", 
#' outfile = "MESA_trained", nvar_use = 500000)
#'
#' ## If data are separated by chromosomes 
#' JHS_model = train_deepRsq(file = paste0("JHS_1000G_chr",1:22,".txt.gz"), 
#' outfile = "JHS_trained", nvar_use = 100000)
#'
#' @importFrom data.table fread 
#'
#' @export

train_deepRsq = function(file, header = T, outfile = NULL, nvar_use = NULL, p_train = 0.7,
 trueR2Col = 5, FeatureCols = 6:86, seed = 123, nrounds = 3000, verbose = 1, early_stopping_rounds = 50,
 print_every_n = 20) {

set.seed(seed)

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

# use a subset of variants to train model

if(is.null(nvar_use)){
message("Not specified number of variants used for training. 
Will use full variants to train models. 
Note that it may take very long time due to computational burden.")
nvar_use = dim(dat)[1]
}

## random subset nvar_train variants
nvar_use = floor(as.numeric(nvar_use))
index = sample(1:nvar_use)
dat = dat[index[1:nvar_use],]

# Split into training (p_train) and testing (1-p_train) data. 
index = sample(1:nvar_use)
train_size = floor(as.numeric(p_train)*nvar_use)
train_index = index[1:train_size]
test_index = index[(train_size+1):nvar_use]

train = dat[train_index,]
test = dat[test_index,]

dtrain = xgb.DMatrix(data = as.matrix(train[,FeatureCols]),label = train[,trueR2Col])
dtest = xgb.DMatrix(data = as.matrix(test[,FeatureCols]),label = test[,trueR2Col])

bst = xgb.train(data = dtrain, booster="gbtree",
                  watchlist =list(eval=dtest),
                  max.depth = 6,
                  nrounds = nrounds,
                  early_stopping_rounds = early_stopping_rounds,
                  maximize = FALSE,
                  verbose  = verbose,
                  print_every_n = print_every_n,
                  objective = "reg:logistic")

if(!is.null(outfile)){
message(paste("Saving model to",outfile))
xgb.save(bst, outfile)
}

return(bst)

}


