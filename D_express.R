### COMBINING FEATURE COUNTS INTO ONE SINGLE TABLE ###

fetch_combine <- function(path, pattern){
  setwd(path)
  files <- list.files(path, pattern)
  nfiles <- length(files)
  file <- read.table(files[1], 
                     sep = "\t", 
                     header = TRUE)
  
  data_counts <- as.data.frame(file[1][1])
  for(i in 1:nfiles){
    file_data <- read.table(files[i], 
                            sep = "\t", 
                            header = TRUE)
    
    data_counts <-(cbind(data_counts,
                         file_data[,ncol(file_data)]))
  }
  colnames(data_counts) <- c("GeneID", files)
  #write.csv(as.data.frame(data_counts),
  #          file="combined_counts.csv")                               # uncomment to save the combined counts into a csv file
  data_matrix<-as.matrix(data_counts[,-1])
  rownames(data_matrix)<-data_counts[,1]
  data_matrix
}

### USING THE COMBINED COUNTS TO PERFORM DESEQ AND IDENTIFY DIFFERENTIALLY EXPRESSED GENES ###

library(DESeq2)
D_express <- function(in_matrix, condition, p_val="."){
  countdata <- in_matrix
  
  (coldata <- data.frame(row.names=colnames(countdata), 
                         condition))
  
  dds <- DESeqDataSetFromMatrix(countData=countdata, 
                                colData=coldata, 
                                design=~condition)
  dds <- DESeq(dds)
  res<-results(dds)
  if (isTRUE(p_val != ".")){
    res_filt <- res[ which(res$padj < p_val), ]                          # filtering based on p-val cut off of 1e-10 and saving it into a variable
    resOrdered <- res_filt[order(-res_filt$log2FoldChange),]             # ordering based on fold change values
    print("Filtered DEseq:")
    print(resOrdered)
  }
  else{
    print("DEseq RESULTS:")
    print(res)
  }
  
  #write.csv(as.data.frame(resOrdered),
  #          file="A_vs_B.csv")                                        # uncomment to save the DEseq results as a csv file
  
  dds
}

### USING THE DESEQ RESULTS TO PLOT A HEATMAP ###

library("pheatmap")
heat_mapping <- function(dds, cutoff=0){
  
  if (cutoff>0){
    select <- order(rowMeans(counts(dds,normalized=TRUE)), 
                    decreasing=TRUE)[1:cutoff]
  }
  else{
    select <- order(rowMeans(counts(dds,normalized=TRUE)), 
                    decreasing=TRUE)
  }
  
  nt <- normTransform(dds)                                             # defaults to log2(x+1)
  log2.norm.counts <- assay(nt)[select,]
  
  df <- as.data.frame(colData(dds)[,c("condition", "sizeFactor")]) 
  df$sizeFactor <- NULL                                                # removing sizeFactor column out of colData(dds)
  
  pheatmap(log2.norm.counts, 
           cluster_rows=FALSE, 
           show_rownames=TRUE,
           cluster_cols=FALSE, 
           annotation_col=df, 
           filename = "heatmap_deseq.pdf", 
           fontsize = 6)
}

### MAIN ###

file_format = "*.txt"                                                  # change for different formats
run_path = "/Users/sridhs01/Desktop/Waxman/RNA_counts/count.geneID"    # change for different working directories

datain <- fetch_combine(run_path, file_format)
colnames(datain)                                                       # prints only the first 10 lines, change if need to view more

condition <- factor(c("PCMV", "PCMV", "PCMV", "PF1", "PF1", "PF1"))    # consider the colnames printed in the above statement before setting an order to the conditions
heat_mapping(D_express(datain, condition, 0.5),3)                      # change p-adj val threshold and number of genes cutoff as per requirement
                                                                       # number of genes cutoff is the number of genes included in the heatmap.
                                                                       # It is not mandatory to mention this option.
