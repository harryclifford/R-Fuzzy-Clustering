
####################################################
# Script to run fuzzy clustering on Normalized EData
# With Automated Cluster Number Selection
####################################################

#minimum cluster number = 2, maximum = 25


#####################################################


## user-defined info

# working directory
working_dir <- "/net/isi-scratch/"

# input expression file (from working directory)
# (with row and col names. Rows are genes, columns are samples)
infile <- "normalized_counts.txt"

## relationship of samples
# give table of sample names, where each row contains samples to be grouped
# and each row name is name of group
sample_relationships <- rbind(c("P0_R1","P0_R2","P0_R3"),c("P4_R1","P4_R2","P4_R3"),c("P8_R1","P8_R2","P8_R3"),c("P14_R1","P14_R2","P14_R3"),c("P21_R1","P21_R2","P21_R3"))

sample_relationships <- rbind(c("P4_MMB1.IgL","P4_MMB1.L4","P4_MMB1.SgL","P4_MMB4.IgL","P4_MMB4.L4","P4_MMB4.SgL"),
                                c("P6_MMB5.IgL","P6_MMB5.L4","P6_MMB5.SgL","P6_MMB6.IgL","P6_MMB6.L4","P6_MMB6.SgL"),
                                c("P8_MMB10.IgL","P8_MMB10.L4","P8_MMB10.SgL","P8_MMB9.IgL","P8_MMB9.L4","P8_MMB9.SgL"),
                                c("P10_MMB13.IgL","P10_MMB13.L4","P10_MMB13.SgL","P10_MMB14.IgL","P10_MMB14.L4","P10_MMB14.SgL"),
                                c("P14_MMB17.IgL","P14_MMB17.L4","P14_MMB17.SgL","P14_MMB18.IgL","P14_MMB18.L4","P14_MMB18.SgL")
                                )
rownames(sample_relationships) <- c("P4","P6","P8","P10","P14")


## alpha core thresholds to report
# when reporting which genes form cluster cores, minimum thresholds for membership values
# give in form of vector of numbers between 0 and 1
acore_vector <- seq(0.5,0.9,0.1)

# maximum number of threads
nthreads <- 15


#####################################################


# ensure installation of packages

library("Mfuzz")
library("marray")
library("fpc")
library("animation")
library("vegan")
library("GMD")
library("foreach")
library("doMC")


#####################################################
#####################################################
#####################################################


## Main script


# read in data
setwd(working_dir)
raw_data <- read.table(infile,stringsAsFactors=F,header=T,row.names=1)

# combines data of same group, using mean
combined_data <- as.data.frame(matrix( rep(0,nrow(raw_data)*nrow(sample_relationships)) , nrow=nrow(raw_data), ncol=nrow(sample_relationships) ))
colnames(combined_data) <- rownames(sample_relationships)
rownames(combined_data) <- rownames(raw_data)

for(i in 1:nrow(sample_relationships)){
    tmp_data <- raw_data[,sample_relationships[i,]]
    combined_data[,i] <- apply(tmp_data,1,mean)
}


# filter and refill (using weighted k-nearest neighbour method) NA data
expr_data <- ExpressionSet(assayData=as.matrix(combined_data))
expr_data.r <- filter.NA(expr_data,thres=0.25)
expr_data.f <- fill.NA(expr_data.r,mode="wknn")

# data standardized to have mean value of zero and SD of one
expr_data.s <- standardise(expr_data.f)

# filters any rows of zero variance and rows of NAs 
exprs(expr_data.s) <- exprs(expr_data.s)[ complete.cases(exprs(expr_data.s)) ,]
exprs(expr_data.s) <- exprs(expr_data.s)[! apply(exprs(expr_data.s), 1, var) == 0 ,]

# records genes removed
total_genes <- nrow(exprs(expr_data))
genes_remaining <- nrow(exprs(expr_data.s))

write.table(row.names(exprs(expr_data.s)),file="background.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)

gene_removal_1 <- total_genes - nrow(exprs(expr_data.r))
gene_removal_2 <- nrow(exprs(expr_data.r)) - genes_remaining

qc_report <- cbind(c("all_genes","genes_removed_due_to_NAs","genes_removed_due_to_zero_variance","total_genes_remaining"),c(total_genes,gene_removal_1,gene_removal_2,genes_remaining))
write.table(qc_report,file="QC_REPORT.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)

# produces colorbar for use with cluster results
png("COLORBAR.png",width=100)
suppressWarnings(mfuzzColorBar())
dev.off()




############################################################

### AUTO-DETERMINES CLUSTER NUMBER(S)


# data prep
dir.create("cluster_number",showWarnings=FALSE)
d <- exprs(expr_data.s)


## silhouette width (method two of stackoverflow)

# methods
pamk.best <- pamk(d)
pamdata <- pam(d, pamk.best$nc)

# plot 1
pdf("cluster_number/silhouettes_details.png")
plot(silhouette(pamdata),main="Silhouette Width for each Cluster")
dev.off()
ani.options(autobrowse=F)
# (has to generate pdf then convert, due to bug in plotting silhouettes)
im.convert("cluster_number/silhouettes_details.png",
            output="cluster_number/silhouettes_details.png",extra.opts="-density 100")

silwidths <- pamk.best$crit[2:10]
cols <- rep("black",9)
cols[pamk.best$nc-1] <- "red"

# plot 2
png("cluster_number/silhouette_selection.png")
plot(2:10,silwidths,type="b",col=cols,pch=19,
        ylim=c(min(silwidths)*0.9,max(silwidths)*1.1),
        main="Average silhouette width with\nincreasing number of clusters",
        xaxt="n",xlab="number of clusters",ylab="Average silhouette width"
    )
axis(1,at=2:10)
legend("topright","optimum",col="red",pch=19)
dev.off()

# optimum cluster number
optimum_sils <- pamk.best$nc


## calinski criterion - takes a very long time to run, and cannot be done in parallel

# methods
fit <- cascadeKM(scale(d, center = TRUE,  scale = TRUE), 1, 10, iter = 1000)

# plot 1
png("cluster_number/optimum_calinski.png",width=700)
plot(fit,sortg =T,grpmts.plot=T)
dev.off()

# optimum cluster number
optimum_group <- names(which.max(fit$results[2,][which(fit$results[2,]!=Inf)]))
optimum_calinski <- as.numeric(gsub("[^0-9]","",optimum_group))


## cluster number selection

cluster_vector <- c(optimum_sils,optimum_calinski)
# minimum of 2, maximum of 25
cluster_vector[which(cluster_vector>25)] <- 25
cluster_vector[which(cluster_vector<2)] <- 2
cluster_vector <- sort(unique(cluster_vector))


# justify cluster number selection through elbow

fuzzfactor.m <- round(mestimate(expr_data.s),2)
dist.obj <- dist(d)

registerDoMC(nthreads)
exp_vars <- foreach(i=2:25) %dopar% {
    hardcluster <- mfuzz(expr_data.s,c=i,m=fuzzfactor.m)$cluster
    css.obj <- css(dist.obj,hardcluster)
    exp_var <- css.obj$totbss/css.obj$tss
    return(exp_var)
}


# plots elbow and highlights cluster_vector choices

col_vector <- rep("black",25)
col_vector[cluster_vector] <- "red"

png("cluster_number/cluster_number_chosen.png")
plot(1:25,c(0,unlist(exp_vars)*100),pch=19,col=col_vector,xlab="Number of Clusters",ylab="Variance Explained (%)",
        main="Elbow Plot")
legend("bottomright","Cluster Number(s) Selected",col="red",pch=19)
dev.off()


## saves data so far
save.image(file="cluster_number/cluster_number.RData")


############################################################


# loops through cluster numbers
mfrow_list <- c( "","1,2","1,3","2,2","2,3","2,3","2,4","2,4","3,3","3,4","3,4","3,4","4,4","4,4","4,4","4,4","5,4","5,4","5,4","5,4","5,5","5,5","5,5","5,5","5,5" )
for(cluster_number in cluster_vector){
    
    # creates output directory
    cluster_dir <- paste(cluster_number,"clusters",sep="")
    dir.create(cluster_dir,showWarnings=FALSE)
    
    # calculates and plots clusters
    cl <- mfuzz(expr_data.s,c=cluster_number,m=fuzzfactor.m)
    tmp_mfrow <- as.numeric(unlist(strsplit(mfrow_list[cluster_number],",")))
    
    png(filename=paste(cluster_dir,"cluster_plots.png",sep="/"),width=400*tmp_mfrow[2],height=400*tmp_mfrow[1])
    mfuzz.plot(expr_data.s,cl=cl,mfrow=tmp_mfrow,time.labels=colnames(combined_data),new.window=FALSE)
    dev.off()
    
    # investigates overlap between clusters
    O <- overlap(cl)
    
    tmp_prcomp <- prcomp(cl[[1]],scale=TRUE)
    xlabel <- paste("(",round(summary(tmp_prcomp)$importance[2,1]*100),"% of Variance)",sep="")
    ylabel <- paste("(",round(summary(tmp_prcomp)$importance[2,2]*100),"% of Variance)",sep="")
    
    png(filename=paste(cluster_dir,"clusters_PCA_plot.png",sep="/"))
    bin <- overlap.plot(cl,over=O,thres=0.05)
    title(xlab=xlabel,line=4)
    title(ylab=ylabel,line=2)
    title(main="PCA of Cluster Centers")
    title(main="(overlap visualised through lines, where width represents strength of overlap)",line=0.8,cex.main=0.8)
    dev.off()
    
    # writes genes forming the alpha cores of the soft clusters
    # loops through different threshold writing results for each
    acore_dir <- paste(cluster_dir,"alpha_core_genes",sep="/")
    dir.create(acore_dir,showWarnings=FALSE)
    for(tmp_acore in acore_vector){
        
        # obtains gene memberships
        tmp_acore_dir <- paste(acore_dir,"/",tmp_acore,"_membership_threshold_genes",sep="")
        dir.create(tmp_acore_dir,showWarnings=FALSE)
        cluster_genes <- acore(expr_data.s,cl,min.acore=tmp_acore)
        
        # loops through each cluster, writing to table
        for(i in 1:length(cluster_genes)){ write.table(cluster_genes[[i]],file=paste(tmp_acore_dir,"/cluster",i,"genes.txt",sep=""),quote=FALSE,sep="\t",row.names=FALSE) }
        
    }
    
}




