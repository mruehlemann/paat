library(phyloseq)
library(vegan)
library(ape)
source("https://github.com/mruehlemann/paat/raw/master/Rscripts/paat_functions.R")

load("gevers.Robj")
ps

subset_depth=3000 
ps.even<-rarefy_even_depth(ps,sample.size=subset_depth,replace=F,rngseed=666)

testvar="diagnosis"
set1=c("no")
set2=c("CD")
covar=c("sex","age")

testdataset<-prune_samples(unlist(sample_data(ps.even)[,testvar]) %in% c(set1,set2), ps.even)
model.specs<-data.frame(sample_data(testdataset))[,c(covar,testvar)]
model.specs$Group<-factor(ifelse(model.specs[,testvar]==set1,"G1","G2"),levels=c("G1","G2"))

phylomat<-branchNodeAdjacency(phy_tree(testdataset))

phylocount<-t(as(otu_table(testdataset),"matrix")) %*% phylomat

presence_thresh=0.25
abu_thresh_upper=0.5
abu_thresh_lower=0.001
phylocount.filter<-filter_abundance(abutab=phylocount, groups=model.specs$Group, min.mean=abu_thresh_lower, max.mean=abu_thresh_upper, group.min.presence=presence_thresh)

phylomat.filter<-phylomat[,colnames(phylocount.filter)]
phylomat.filter.cor<-cor(phylomat.filter)
diag(phylomat.filter.cor)<-0
phylomat.filter.cor[upper.tri(phylomat.filter.cor)]<-0
phylomat.filter.cor<-ifelse(phylomat.filter.cor>0,1,0)

phylocount.bc<-1-as.matrix(vegdist(t(phylocount.filter)))
diag(phylocount.bc)<-0
phylocount.filter<-phylocount.filter[,apply((phylocount.bc * phylomat.filter.cor),2,function(x) any(x>0.95))==F]

phylomat.jc<-1-as.matrix(vegdist(t(phylomat.filter),"jaccard"))
excl<-names(which(apply(phylomat.jc*phylomat.filter.cor,2,max)<0.8 & colSums(phylomat.filter)>(nrow(phylomat.filter)*0.25)))
phylocount.filter<-phylocount.filter[,!(colnames(phylocount.filter) %in% unique(c(excl,names(which(apply(phylomat.filter.cor[c(excl),,drop=F]==1,2,any))))))]

outmat<-test_abundance(phylocount.filter,model.specs,testingvar="Group",covars=covar,zerothresh=0.05)

phylocount.filter.rel<-phylocount.filter/subset_depth
outmat<-test_abundance_asinsqrt(phylocount.filter.rel,model.specs,testingvar="Group",covars=covar)
outmat$p.adj<-p.adjust(outmat$p,"fdr")
outmat.phylo<-outmat[outmat$p.adj<0.05 & is.na(outmat$p.adj)==F,]
nrow(outmat.phylo)

phylomat.filter.cor2<-phylomat.filter.cor
diag(phylomat.filter.cor2)<-1
phylomat.filter.cor2.sub<-phylomat.filter.cor2[rownames(outmat.phylo),rownames(outmat.phylo)]
outmat.phylo$cluster<-factor(apply(phylomat.filter.cor2.sub,1,function(x) min(which(x==1)))*sign(outmat.phylo$beta))
outmat.phylo$cluster<-as.numeric(factor(outmat.phylo$cluster,levels=unique(outmat.phylo$cluster)))
unique(outmat.phylo$cluster)

outmat.phylo.refined<-outmat.phylo[0,]
for(cl in unique(outmat.phylo$cluster)){
branch_in_cluster<-outmat.phylo[outmat.phylo$cluster==cl,]
    branch_in_cluster<-branch_in_cluster[order(abs(branch_in_cluster$beta),decreasing=T),]
    while(nrow(branch_in_cluster)>0){
      outmat.phylo.refined[nrow(outmat.phylo.refined)+1,]<-branch_in_cluster[1,]
      branch_in_cluster<-branch_in_cluster[-1,]
      branch_in_cluster<-branch_in_cluster[apply(cor(phylomat.filter[,rownames(branch_in_cluster),drop=F],phylomat.filter[,rownames(outmat.phylo.refined),drop=F]),1,function(x) any(x>0))==F,]
    }
  }
nrow(outmat.phylo.refined)

phylocount.final<-phylocount.filter[,rownames(outmat.phylo.refined)]
phylomat.final<-phylomat.filter[,rownames(outmat.phylo.refined)]

outmat.annotated<-annotate_branches(phyloseq=testdataset,results=outmat.phylo.refined, phylomat=phylomat.final)

otu.mat<-as.data.frame(t(otu_table(testdataset)))
otu.mat.sub<-filter_abundance(abutab=otu.mat, groups=model.specs$Group, min.mean=abu_thresh_lower, group.min.presence=presence_thresh)
tips_to_keep<-colnames(otu.mat.sub)
treep<-plot_annotated_tree(phyloseq=testdataset, results=outmat.annotated, phymat=phylomat.final, tips=tips_to_keep)
ggsave(treep, file=paste0("examples/gevers_",paste0(c(set1),collapse=""),".vs.",paste0(c(set2),collapse=""),".paat.pdf"), height=16,width=20)
