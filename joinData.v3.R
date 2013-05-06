
library(gplots)
library(RColorBrewer)
library(ConsensusClusterPlus)
library(ggplot2)

homeDir <- "/data/matheus/labbioinfo/data/LGG"
setwd(homeDir)

source("library.R")

batch.names <- list.files('decompressed', full.name=T)

batches <- list()

for(b in batch.names){
  filenames <- list.files(b, pattern='.*LGG.*', full.name=T)

  files <- list()
  for(f in filenames){
    name <- gsub("(\\S+\\.)(\\S+)(\\.txt)", "\\2", f)
    files[[name]] <- read.table(f, skip=1, header=T, sep='\t')
  }
  
  merged <- files[[1]]
  merged <- cbind(merged[, -which(colnames(merged) %in% 'Beta_value')], 'Beta_value'= merged[,'Beta_value'])
  colnames(merged) <- sub("Beta_value", paste("Beta_value", names(files)[1], sep="."), colnames(merged))
  
  f <- names(files)[2]
  for(f in names(files[-1])){
    merged <- merge(merged, files[[f]][,c('Composite.Element.REF','Beta_value')], by='Composite.Element.REF')
    colnames(merged) <- sub("^Beta_value$", paste("Beta_value", f, sep="."), colnames(merged))
  }
  
  rm(files)
  
  batches[[b]] <- merged
  
}
rm(merged)
## Memory usage 
# batches 1.08   GB

data <- batches[[1]]

for(b in names(batches[-1])){
  cols <- which(colnames(batches[[b]]) %in% 'Composite.Element.REF')
  cols <- c(cols, grep('Beta_value', colnames(batches[[b]])))
  data <- merge(data, batches[[b]][,cols], by='Composite.Element.REF')
}

dim(data)

colnames(data) <- gsub("Beta_value.(\\S+)", "\\1", colnames(data))

rownames(data) <- with(data, paste(Gene_Symbol, Composite.Element.REF, sep="/"))


#save(file="data.RData", list=c("data"))
#load("data.RData")
rm(batches)





manifest <- data.frame()
for( b in batch.names){
  filenames <- list.files(b,pattern='.*LGG.*', full.name=T)
  
  manifest <- rbind(manifest,cbind(
    batch = gsub(pattern="(\\S+Level_\\d+.)(\\d+)(.\\d+.\\d+)(\\S+)", replacement="\\2", x=filenames),   # batches
    level = gsub("(\\S+Level_)(\\d+)(\\S+)", "\\2", filenames),   # level
    sample = gsub("(\\S+.)(TCGA-\\S+)(.txt)", "\\2", filenames)   # sample name
  ))
}
rownames(manifest) <- manifest$sample

tmp <- as.data.frame(matrix(unlist(apply(manifest, 1, function(x) strsplit(x['sample'], "-") )), ncol=7, byrow=T))
colnames(tmp) <- c("project", "tss", "participant", "sample.type.vial", "portion.analyte", "plate", "center")
head(tmp)
tmp$sample.type <- substr(as.character(tmp$sample.type.vial), 1, 2)
tmp$vial <-  substr(as.character(tmp$sample.type.vial), 3, 3)
tmp$portion <- substr(as.character(tmp$portion.analyte), 1, 2)
tmp$analyte <- substr(as.character(tmp$portion.analyte), 3, 3)
manifest <- cbind(manifest, tmp)
head(manifest)
if(is.character(manifest$sample.type)){
  manifest$sample.type <- as.numeric(manifest$sample.type)
}else if(is.factor(manifest$sample.type)){
  manifest$sample.type <- as.numeric(unfactor(manifest$sample.type))
}


showMemoryUse()


## Clinical Information

patient <- read.table('clinical/clinical_patient_lgg.txt', header=T, sep='\t')
rownames(patient) <- patient$bcr_patient_barcode
table(patient$person_neoplasm_cancer_status)
head(patient)
# remove x and y chromosomes
table(data$Chromosome)
data$Chromosome <- as.character(data$Chromosome)
data <- subset(data, Chromosome != "X" & Chromosome != "Y")


smp.bvalue <- colnames(data)[grep("TCGA", colnames(data))]
smp.normal <- as.character(subset(manifest, sample.type >= 10 & sample.type < 20)$sample)
smp.tumor <- as.character(subset(manifest, sample.type < 10 )$sample)

dim(data)
sd.tumor.cut <- 0.3
data$sd.tumor <- apply(data[, smp.tumor], 1, sd, na.rm=T)
p <- ggplot(data, aes(sd.tumor)) + geom_histogram() 
p + geom_vline(xintercept=sd.tumor.cut)
ggsave(filename='histogram.sd.png')
hi.sd.tumor.100 <- rownames(data)[tail(order(data$sd.tumor, na.last=F), n=100)]
row.sd.tumor.cut <- rownames(data)[which(data$sd.tumor>sd.tumor.cut)]
length(row.sd.tumor.cut)
## Heatmap
tmp.batches <- manifest[smp.tumor, "batch"]
colann <- matrix(brewer.pal(7, 'Greens')[tmp.batches], ncol=1)
tmp.patients <- patient[gsub("(\\S+-\\S\\S-\\d+)(\\S*)", "\\1", smp.tumor),c('visual_changes',
  'mental_status_changes', 'headache_history', "gender", "tumor_location", 'family_history_of_cancer',
   'family_history_of_primary_brain_tumor')]
for( i in 1:ncol(tmp.patients)){
  colann <- cbind(colann, brewer.pal(length(levels(tmp.patients[,i])), 'Greens')[tmp.patients[,i]])
}
apply(tmp.patients, 2, function(x){print(levels(x))})
rownames(colann) <- smp.tumor
colnames(colann) <- c("Batches", colnames(tmp.patients))
png('heatmap.png', res=300, height=2000, width=2000)
heatmap.3(as.matrix(data[row.sd.tumor.cut, smp.tumor]),
          scale='none',
          hclustfun = function(x,...){hclust(x, method='average')},
          distfun = function(x,...){dist(x, method='euclidean')},
          dendrogram = "both",
          col=greenred(1000),
          tracecol=F,
          na.color = "grey",
          margins = c(10,5),
          #labRow=data[hi.sd.100, "Gene_Symbol"],
          labCol=smp.tumor,
          ColSideColors=colann,
          key=T,
          keysize=1.8,
          cexRow=0.8,
          main="hi sd",
          )

dev.off()

ccluster = ConsensusClusterPlus(as.matrix(data[row.sd.tumor.cut, smp.bvalue]), maxK=20, reps=1000,
  title="TCGA",clusterAlg="hc",distance="euclidean", plot="png")

## PCA
# pca1 <- prcomp(na.exclude(t(data[,smp.bvalue])))
# plot(pca1$x[,1], pca1$x[,2], col=colann[,1])
# dim(pca1$x)







