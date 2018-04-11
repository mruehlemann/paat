#### convert tree to matrix

branchNodeAdjacency <- function(x) {
    m <- matrix(0, nrow=length(x$tip.label), ncol=nrow(x$edge))
    from <- x$edge[,1]
    to <- x$edge[,2]
    g <- seq_along(x$tip.label)
    while (any(!is.na(g))) {
        i <- match(g, to)
        m[cbind(seq_along(i),i)] <- 1
        g <- from[i]
    }
    colnames(m) <- paste0("B", seq.int(ncol(m)))
    rownames(m) <- x$tip.label
    m
}

### arcsin-sqrt-transform

trans.arcsine <- function(x){
  asin(sign(x) * sqrt(abs(x)))
}

### abundance testing using negbin GLM with < 5% zeroes and zero-truncated negbin GLM with > 5%

test_abundance<-function(abundance_data, sample_specs,testingvar="Group",covars,zerothresh=0.05){
require(MASS)
returnmat<-data.frame(row.names=colnames(abundance_data),method=rep(NA,ncol(abundance_data)),beta=NA,se=NA,p=NA,zero=NA,count=NA)
for(taxon in colnames(abundance_data)){
x<-as.numeric(abundance_data[,taxon])
zero=mean(x==0)
count=mean(x[x!=0])
if(zero>zerothresh){
f1=tryCatch({summary(glm.nb(as.formula(paste("count ~ ", paste(c(covars, testingvar),collapse="+"))),data=data.frame(sample_specs,count=x),subset=x>0))$coefficients[testingvar,c(1,2,4)]},error = function(x) NA,finally={})
returnmat[taxon,"method"]<-"GLMtrunc"
}else{
f1=tryCatch({summary(glm.nb(as.formula(paste("count ~ ", paste(c(covars, testingvar),collapse="+"))),data=data.frame(sample_specs,count=x)))$coefficients[testingvar,c(1,2,4)]},error = function(x) NA,finally={})
returnmat[taxon,"method"]<-"GLM"
}
returnmat[taxon,c("zero","count")]<-c(zero,count)
returnmat[taxon,c("beta","se","p")]<-as.numeric(f1)
}
return(returnmat)
}

### abundance testing using arcsin sqrt transformed data

test_abundance_asinsqrt<-function(abundance_data, sample_specs,testingvar="Group",covars){
returnmat<-data.frame(row.names=colnames(abundance_data),method=rep(NA,ncol(abundance_data)),beta=NA,se=NA,p=NA,zero=NA,count=NA)
for(taxon in colnames(abundance_data)){
x<-as.numeric(abundance_data[,taxon])
zero=mean(x==0)
count=mean(x[x!=0])
f1=tryCatch({summary(lm(as.formula(paste("count ~ ", paste(c(covars, testingvar),collapse="+"))),data=data.frame(sample_specs,count=trans.arcsine(x))))$coefficients[testingvar,c(1,2,4)]},error = function(x) NA,finally={})
returnmat[taxon,"method"]<-"LMasinsqrt"
returnmat[taxon,c("zero","count")]<-c(zero,count)
returnmat[taxon,c("beta","se","p")]<-as.numeric(f1)
}
return(returnmat)
}

