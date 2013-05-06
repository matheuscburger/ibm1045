
batches <- list.files('decompressed', full.name=T)

b1 <- batches[1]

filenames <- list.files(b1, pattern='.*LGG.*', full.name=T)

files <- list()
for(f in filenames){
	files[[f]] <- read.table(f, skip=1, header=T, sep='\t')
}


merged <- files[[1]]
colnames(merged) <- sub("Beta_value", paste("Beta_value", gsub("(\\S+\\.)(\\S+)(\\.txt)", "\\2", names(files)[1]), sep="."), colnames(merged))
f <- names(files)[2]
for(f in names(files[-1])){
	name <- gsub("(\\S+\\.)(\\S+)(\\.txt)", "\\2", f)
	merged <- merge(merged, files[[f]], by='Composite.Element.REF')
	merged <- merged[,c('Composite.Element.REF',colnames(merged)[grep('.x', colnames(merged))], 'Beta_value.y')]
	cnames <- colnames(merged)
	cnames[which(cnames == "Beta_value")] <- paste("Beta_value",name, sep=".")
	cnames <- sub(".x", "", cnames)
	colnames(merged) <- cnames
}

