### filter abundance tables

filter_abundance<-function(abutab, groups, min.mean=0, group.min.mean=0, max.mean=1, group.max.mean=1, min.presence=0, group.min.presence=0, tot=subset_depth){
abutab.means<-colMeans(abutab)/tot
abutab.presence<-colMeans(abutab>0)
abutab.group.means<-apply(aggregate(. ~ groups,data=data.frame(abutab),"mean")[,-1],2,min)/tot
abutab.group.presence<-apply(aggregate(. ~ groups,data=data.frame(ifelse(abutab>0,1,0)),"mean")[,-1],2,min)
abutab.sub<-abutab[,abutab.presence>min.presence & abutab.group.presence>group.min.presence & abutab.means>min.mean & abutab.means<=max.mean & abutab.group.means>group.min.mean & abutab.group.means<=group.max.mean]
return(abutab.sub)
}


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
f1=tryCatch({summary(glm.nb(as.formula(paste("count ~ ", paste(c(covars, testingvar),collapse="+"))),data=data.frame(sample_specs,count=x),subset=x>0))$coefficients[length(covars)+2,c(1,2,4)]},error = function(x) NA,finally={})
returnmat[taxon,"method"]<-"GLMtrunc"
}else{
f1=tryCatch({summary(glm.nb(as.formula(paste("count ~ ", paste(c(covars, testingvar),collapse="+"))),data=data.frame(sample_specs,count=x)))$coefficients[length(covars)+2,c(1,2,4)]},error = function(x) NA,finally={})
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
f1=tryCatch({summary(lm(as.formula(paste("count ~ ", paste(c(covars, testingvar),collapse="+"))),data=data.frame(sample_specs,count=trans.arcsine(x))))$coefficients[length(covars)+2,c(1,2,4)]},error = function(x) NA,finally={})
returnmat[taxon,"method"]<-"LMasinsqrt"
returnmat[taxon,c("zero","count")]<-c(zero,count)
returnmat[taxon,c("beta","se","p")]<-as.numeric(f1)
}
return(returnmat)
}

### branch annotation

annotate_branches<-function(phyloseq, results, phylomat, annot_thresh=.5, show_min=.2, annot_min=.05){
results.out<-results
results.out[,c("tag","tax","ntips")]<-NA
for(signal in rownames(results)){
tips_in_branch<-names(which(phylomat[,signal]==1))
tipcounts.signal<-as.matrix(otu_table(phyloseq))[tips_in_branch,,drop=F]
tipcounts.weights<-rowSums(tipcounts.signal)/sum(tipcounts.signal)
df<-data.frame(tax_table(phyloseq),stringsAsFactors=F)[tips_in_branch,]
df[is.na(df)]<-FALSE
if(ncol(df)==7){df[,7]<-ifelse(df[,7]==F,F,paste(df[,6],df[,7],sep="_"))}
if(length(tips_in_branch)>1){
ta<-apply(df,2,function(x){df2<-data.frame(tax=as.character(x),w=tipcounts.weights);if(all(df2$tax==F)){return(FALSE)};zz=aggregate(w ~ tax,data=df2,sum);notshown=sum(zz$w>annot_min & zz$w<=show_min & zz$tax!=F);zz<-zz[zz$w>show_min & zz$tax!="FALSE",];zz<-zz[order(zz$w,decreasing=T),];
   ifelse(sum(zz$w)>annot_thresh,ifelse(notshown>0,paste0(paste(zz$tax,collapse=","),"(+",notshown,")"),paste(zz$tax,collapse=",")),"")})
}else{
ta<-df
}
ta<-ta[is.na(ta)==F & ta!=F & ta!=""]
results.out[signal,c("tag","tax","ntips")]<-c(names(which.max(tipcounts.weights)),ta[length(ta)],length(tips_in_branch))
}
return(results.out)
}


### plotting 

plot_annotated_tree<-function(phyloseq, results, phymat, tips, expansion=10){
require(ggtree)
phyloseq.sig<-prune_taxa(taxa_names(phyloseq) %in% c(tips, results$tag), phyloseq)
phyloseq.sig.tree<-phy_tree(phyloseq.sig)
treeplot<-ggtree(phyloseq.sig.tree,layout="circular",branch.length="none")+scale_x_continuous(expand=c(0,expansion))+geom_tiplab2()
for(signal in rownames(results)){
entry<-results[signal,]
tips_in_branch<-names(which(phymat[,signal]==1))
tips_in_plot=tips_in_branch[tips_in_branch %in% taxa_names(phyloseq.sig)]
color=ifelse(entry$beta<0,"steelblue","firebrick2")
thisalpha=0.35+(abs(entry$beta)/max(abs(results$beta)))*0.5
if(length(tips_in_plot)>1){
nodenum<-MRCA(phyloseq.sig.tree,tips_in_plot) 
treeplot<-treeplot+geom_hilight(node=nodenum,fill=color,alpha=thisalpha)
treeplot<-treeplot+geom_cladelabel(node=nodenum, label=paste0(signal,"-",entry$tax,";",entry$ntips), align=F,offset=expansion/2,hjust=ifelse(treeplot$data$angle[nodenum]<90 | treeplot$data$angle[nodenum]>270,0,1))
}else{
nodenum<-which(treeplot$data$label==tips_in_plot)
treeplot<-treeplot+geom_point2(aes_string(subset=paste("node ==",nodenum)), fill=color,size=3, shape=ifelse(length(tips_in_branch)>1,24,21),alpha=thisalpha)
treeplot<-treeplot+geom_cladelabel(node=nodenum, label=paste0(signal,"-",entry$tax,";",entry$ntips), align=T,offset=expansion/2,hjust=ifelse(treeplot$data$angle[nodenum]<90 | treeplot$data$angle[nodenum]>270,0,1))
}
}
return(treeplot)
}

