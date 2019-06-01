# GTEX data analysis
rm(list=ls())

library(glmnet)
library(clustermq)

load("all_tissues.preprocessed.data.RData.gz")

# Intersect the genes across all data sets

intersect_names <- colnames(datasets[["Esophagus_Muscularis"]]$expr.normal)
for (name in names(datasets)) intersect_names <- intersect(intersect_names,colnames(datasets[[name]]$expr.normal))
length(intersect_names)

for (name in names(datasets)) {
  datasets[[name]]$expr.normal <- datasets[[name]]$expr.normal[,intersect_names]
  datasets[[name]]$expr.unconfounded <- datasets[[name]]$expr.normal[,intersect_names]
}

# Save the data of the tissues individually (so that tissues can be loaded separately with less memory consumption)

for (name in names(datasets)) {
  normal <- datasets[[name]]$expr.normal
  covariates <- datasets[[name]]$covariates
  save(normal,covariates,file = name)
}

########## 

#for (gamma_max in c(.25,.5,1,2,4,8,16)){
#for (gamma_type in c("anchor_stab","gamma_large")){
gamma_max <- .5
gamma_type <- "anchor_stab"

#dataset_name <- names(datasets)[1]

which_reproducible <- function(dataset_name,j,gamma_max,gamma_type){
  load(dataset_name)
  covariates$ID <- NULL
  library(glmnet)
  
  fit_const <- lm(normal ~ 1)
  fit <- lm(normal ~ t(covariates))
 
  # find lambda via cross-validation 
  gamma <- 0
  newdata <- fit_const$fitted.values + fit$residuals + sqrt(gamma)*(fit$fitted.values-fit_const$fitted.values)
  indices <- 1:nrow(newdata)
  fit_glmnet <- cv.glmnet(x = newdata[indices,-c(j)],newdata[indices,j])
  lambda_cv <- fit_glmnet$lambda.1se
  
  gamma_vec <- (0:10/10*gamma_max)
  mat <- NULL
  stability <- function(i){  
    set.seed(i)
    if (gamma_type == "anchor_stab") { gamma <- gamma_vec[i] } else { gamma <- gamma_max } 
    
    # create artificial data
    newdata <- fit_const$fitted.values + fit$residuals + sqrt(gamma)*(fit$fitted.values-fit_const$fitted.values)
    
    # fit anchor regression
    fit_glmnet_new <- glmnet(x = newdata[indices,-c(j)],newdata[indices,j],lambda = lambda_cv)
    
    # fit lasso on normalized data
    newdata <- fit_const$fitted.values + fit$residuals
    fit_glmnet_0 <- glmnet(x = newdata[indices,-c(j)],newdata[indices,j],lambda = lambda_cv)
    
    return(c(as.vector(coef(fit_glmnet_new)),as.vector(coef(fit_glmnet_0))))
  }
  mat <- lapply(X = 1:length(gamma_vec),FUN=stability)
  mat <- simplify2array(mat)
  
  # compute minimal coefficients over anchor path (in the lasso case, the path is constant)
  min_vec <- apply(abs(mat),1,min)
  
  # split in anchor half and lasso half
  min_vec_anchor <-  min_vec[1:(length(min_vec)/2)]
  min_vec_lasso <- min_vec[(length(min_vec)/2+1):length(min_vec)]
  # remove intercept
  min_vec_anchor <- min_vec_anchor[-c(1)]
  min_vec_lasso <- min_vec_lasso[-c(1)]
  
  # sort nonzero effects and prepare return object
  ret_anchor <- order(as.vector(min_vec_anchor),decreasing = TRUE)[1:sum(min_vec_anchor>0)]
  ret_lasso  <- order(as.vector(min_vec_lasso),decreasing = TRUE)[1:sum(min_vec_lasso>0)]
  sort_anchor <- sort(as.vector(min_vec_anchor),decreasing = TRUE)[1:sum(min_vec_anchor>0)]
  sort_lasso  <- sort(as.vector(min_vec_lasso),decreasing = TRUE)[1:sum(min_vec_lasso>0)]
  
  return(list(lasso=ret_lasso,anchor=ret_anchor,sort_anchor=sort_anchor,sort_lasso=sort_lasso)) 
}

#gamma_max <- 2
#j_seq <- 1:10
set.seed(1)
j_seq <- sample(1:ncol(datasets$Adipose_Subcutaneous$expr.normal),200,replace = FALSE)

grid <- expand.grid(names(datasets),j_seq,stringsAsFactors = FALSE)
start <- Sys.time()
result_list <- Q(fun = which_reproducible, dataset_name = grid[,1],j = grid[,2],const=list(gamma_type = gamma_type, gamma_max=gamma_max),n_jobs = 250)
print(Sys.time() - start)

#result_name <- paste("result_","gamma_type","_",gamma_max,sep="")
#save(result_list,file = result_name)

#################### Create plots ####################################################

#name1 <- names(datasets)[1]
#name2 <- names(datasets)[2]
compare <- function(vec1,vec2){
  intersect_first_k <- function(k){ length(intersect(vec1[1:k],vec2[1:k]))}
  return(sapply(1:20,intersect_first_k))
}

#  comparison lasso-anchor anchor-anchor lasso-lasso
lasso_anchor <- NULL
anchor_anchor <- NULL
lasso_lasso <- NULL
# go through target genes
for (j in j_seq){
  # compare results on two different tissues
  for (name1 in names(datasets))
    for (name2 in setdiff(names(datasets),name1)){
      first_result <- result_list[[which((name1 == grid[,1]) & (j == grid[,2]))]]
      second_result <- result_list[[which( (name2 == grid[,1]) & (j == grid[,2]))]]
      
      lasso_anchor <- cbind(lasso_anchor,compare(first_result$lasso,second_result$anchor))
      anchor_anchor <- cbind(anchor_anchor,compare(first_result$anchor,second_result$anchor))
      lasso_lasso <- cbind(lasso_lasso,compare(first_result$lasso,second_result$lasso))
    }
}
# common effects, summed over the other tissues
lasso_lasso <- rowMeans(lasso_lasso)*12
lasso_anchor <- rowMeans(lasso_anchor)*12
anchor_anchor <- rowMeans(anchor_anchor)*12



pdf_name <- paste("replicability_",gamma_type,"_",gamma_max,".pdf",sep="")
pdf(pdf_name,width=8,height=6) 
matplot(1:20,t(rbind(anchor_anchor,lasso_anchor,lasso_lasso)),type="l",lwd=4,xlab="K", ylab="number of replicable features on a different tissue",cex=2.0)
legend("topleft", c( "anchor regression - anchor regression", "lasso - anchor regression","lasso - lasso"),col=c(1,2,3),cex=1.0,fill=c(1,2,3))
dev.off()
#}

# aggregated across j
aggregate_lasso <- NULL
aggregate_anchor <- NULL
aggregate_sort_lasso <- NULL
aggregate_sort_anchor <- NULL


# combine results for all target genes
for (name1 in names(datasets)) {
  for (j in j_seq){
    first_result <- result_list[[which((name1 == grid[,1]) & (j == grid[,2]))]]
    
    aggregate_lasso[[name1]] <- c(aggregate_lasso[[name1]],paste(first_result$lasso,j))
    aggregate_anchor[[name1]] <- c(aggregate_anchor[[name1]],paste(first_result$anchor,j))
    aggregate_sort_lasso[[name1]] <- c(aggregate_sort_lasso[[name1]],first_result$sort_lasso)
    aggregate_sort_anchor[[name1]] <- c(aggregate_sort_anchor[[name1]],first_result$sort_anchor)
    
  }
}

compare <- function(vec1,vec2){
  intersect_first_k <- function(k){ length(intersect(vec1[1:k],vec2[1:k]))}
  return(sapply(1:100,intersect_first_k))
}

#  comparison lasso-anchor anchor-anchor lasso-lasso
lasso_anchor <- NULL
anchor_anchor <- NULL
lasso_lasso <- NULL
# compare results on two different tissues
for (name1 in names(datasets)) { 
  for (name2 in setdiff(names(datasets),name1)){
    
    first_result_lasso <- aggregate_lasso[[name1]][order(aggregate_sort_lasso[[name1]],decreasing=TRUE)]
    first_result_anchor <- aggregate_anchor[[name1]][order(aggregate_sort_anchor[[name1]],decreasing=TRUE)]
    second_result_lasso <- aggregate_lasso[[name2]][order(aggregate_sort_lasso[[name2]],decreasing=TRUE)]
    second_result_anchor <- aggregate_anchor[[name2]][order(aggregate_sort_anchor[[name2]],decreasing=TRUE)]
    
    lasso_anchor <- cbind(lasso_anchor,compare(first_result_lasso,second_result_anchor))
    anchor_anchor <- cbind(anchor_anchor,compare(first_result_anchor,second_result_anchor))
    lasso_lasso <- cbind(lasso_lasso,compare(first_result_lasso,second_result_lasso))
  } }
lasso_lasso <- rowMeans(lasso_lasso)
lasso_anchor <- rowMeans(lasso_anchor)
anchor_anchor <- rowMeans(anchor_anchor)

pdf_name <- paste("replicability_acrossj_",gamma_type,"_",gamma_max,".pdf",sep="")
pdf(pdf_name,width=8,height=6) 
matplot(1:100,t(rbind(anchor_anchor,lasso_anchor,lasso_lasso)),type="l",lwd=4,xlab="K", ylab="number of replicable features on a different tissue",cex=2.0)
legend("topleft", c( "anchor regression - anchor regression", "lasso - anchor regression","lasso - lasso"),col=c(1,2,3),cex=1.0,fill=c(1,2,3))
dev.off()

#}}
