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


res.otu<-test_abundance_asinsqrt(abundance_data=otu.mat.sub/subset_depth, sample_specs=model.specs, testingvar="Group", covars=covar)
res.otu$p.adj<-p.adjust(res.otu$p,"fdr")
outmat.otu<-res.otu[res.otu$p.adj<0.05 & is.na(res.otu$p.adj)==F,]
nrow(outmat.otu)

phylomat.otu<-ifelse(cor(otu.mat.sub)==1,1,0)
outmat.otu.annotated<-annotate_branches(phyloseq=testdataset,results=outmat.otu, phylomat=phylomat.otu)

treep.otu<-plot_annotated_tree(phyloseq=testdataset, results=outmat.otu.annotated, phymat=phylomat.otu, tips=c(tips_to_keep,outmat.annotated$tag))
ggsave(treep.otu, file=paste0("examples/gevers_",paste0(c(set1),collapse=""),".vs.",paste0(c(set2),collapse=""),".otu.pdf"), height=16,width=20)



tax.gen<-unique(tax_table(testdataset)[,1:6])
tt.gen<-tax_table(testdataset)[,1:6]
alltax<-apply(tax.gen,1,paste,collapse=";")
cl<-data.frame(id=rownames(tt.gen), cluster=sapply(apply(tt.gen,1,paste,collapse=";"),function(x) match(x, alltax)))
library(reshape2)
ii<-data.frame(dcast(id ~ cluster,data=cl,fill=0),row.names=1)
ii[ii>0]<-1
ii<-ii[colnames(otu.mat),]
colnames(ii)<-as.character(alltax)

gen.mat<-as.matrix(otu.mat) %*% as.matrix(ii)
gen.mat.sub<-filter_abundance(abutab=gen.mat, groups=model.specs$Group, min.mean=abu_thresh_lower, group.min.presence=presence_thresh)
res.gen<-test_abundance_asinsqrt(abundance_data=gen.mat.sub/subset_depth, sample_specs=model.specs, testingvar="Group", covars=covar)
res.gen$p.adj<-p.adjust(res.gen$p,"fdr")
outmat.gen<-res.gen[res.gen$p.adj<0.05 & is.na(res.gen$p.adj)==F,]

otu_to_siggen<-names(which(rowSums(ii[,rownames(outmat.gen)])>0))
res.gen<-res.otu[rownames(res.otu) %in% otu_to_siggen,]

outmat.gen.annotated<-res.gen
outmat.gen.annotated$tag<-rownames(outmat.gen.annotated)
outmat.gen.annotated$tax<-apply(tt.gen[rownames(outmat.gen.annotated),],1,function(x) x[max(which(is.na(x)==F))])
outmat.gen.annotated$ntips<-colSums(ii)[match(apply(tt.gen[rownames(outmat.gen.annotated),],1,paste, collapse=";"),alltax)]

treep.gen<-plot_annotated_tree(phyloseq=testdataset, results=outmat.gen.annotated, phymat=phylomat.otu, tips=c(tips_to_keep,outmat.annotated$tag))
ggsave(treep.gen, file=paste0("examples/gevers_",paste0(c(set1),collapse=""),".vs.",paste0(c(set2),collapse=""),".gen.pdf"), height=16,width=20)

