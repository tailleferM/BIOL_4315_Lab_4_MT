BiocManager::install("clusterProfiler")
BiocManager::install("EnhancedVolcano")
BiocManager::install('biomaRt')
BiocManager::install('Rsubread')

lapply(c(
  "docopt", "DT", "pheatmap","GenomicFeatures", "DESeq2",
  "edgeR", "systemPipeR", "systemPipeRdata", "BiocStyle", "GO.db", "dplyr",
  "tidyr", "stringr", "Rqc", "QuasR", "DT", "ape", "clusterProfiler","biomaRt","EnhancedVolcano"
), require, character.only = TRUE)

##### 5.1 Accessing the RNAseq data and sample data #####

#e.g.1 get the path 
fq_path <- systemPipeRdata::pathList()$fastqdir

#e.g.2 list the files 

fq_files <- list.files(fq_path)

print(fq_files)

# import meta data
meta_data <- read.table(
  "/Users/mtaillefer00/Library/R/arm64/4.4/library/systemPipeRdata/extdata/param/targetsPE.txt", 
  header = TRUE)
#replace string with correct path
meta_data$FileName1 <- gsub('./data/',fq_path,meta_data$FileName1)
meta_data$FileName2 <- gsub('./data/',fq_path,meta_data$FileName2)
#display data table
library('DT')
DT::datatable(meta_data)

##### 5.2 Read processing #####

qc_res <- rqc(path = fq_path, pattern = '.fastq.gz', openBrowser = FALSE)

rqcCycleQualityBoxPlot(qc_res[0:12])+theme_bw()
rqcCycleQualityBoxPlot(qc_res[13:24])+theme_bw()
rqcCycleQualityBoxPlot(qc_res[25:36])+theme_bw()

rqcCycleBaseCallsLinePlot(qc_res[0:12])
rqcCycleBaseCallsLinePlot(qc_res[13:24])
rqcCycleBaseCallsLinePlot(qc_res[25:36])

#question 3

meta_data$processed1 <- sub("^/Users/mtaillefer00/Library/R/arm64/4.4/library/systemPipeRdata/extdata/fastq/", "/Users/mtaillefer00/Documents/BIOL_4315_Lab_4_MT/processed_reads/", meta_data$FileName1) 
meta_data$processed1 <- sub("(\\.fastq\\.gz)$", "_processed\\1", meta_data$processed1)  
meta_data$processed2 <- sub("^/Users/mtaillefer00/Library/R/arm64/4.4/library/systemPipeRdata/extdata/fastq/", "/Users/mtaillefer00/Documents/BIOL_4315_Lab_4_MT/processed_reads/", meta_data$FileName2) 
meta_data$processed2 <- sub("(\\.fastq\\.gz)$", "_processed\\1", meta_data$processed2) 

for(i in 1:nrow(meta_data)){
  QuasR::preprocessReads(filename = meta_data$FileName1[i],outputFilename = meta_data$processed1[i], truncateEndBases = 3, Lpattern = 'GCCCGGGTAA', nBases = 1)
} 

##### 5.3 Alignments #####

#5.3.1 Indexing the reference genome

#create the output directory
dir.create("./hisat2_index", recursive = TRUE)

at_genome <- "data/GCF_000001735.4_TAIR10.1_genomic.fna"

#Use system2 to run hisat2 from within R

tryCatch({system2(command = "hisat2-build", 
                  args = c("-p","8", at_genome,
                           "./hisat2_index/tair10_1_index"),
                  stdout = TRUE, stderr = TRUE)}, error = function(e) {
                    paste("hisat2-build", "indexing failed with error:", e$message)
                  })

#5.3.2 Mapping the reads to the indexed genome
#question 4

dir.create("./sam_files",recursive = TRUE)
dir.create("./bam_files",recursive = TRUE)

#setup the file names
meta_data$processedUnzip1 <- sub("\\.gz$", "", meta_data$processed1)
meta_data$sam1 <- paste('/Users/mtaillefer00/Documents/BIOL_4315_Lab_4_MT/sam_files/',meta_data$SampleName, sep = '')
meta_data$sam1 <- paste0(meta_data$sam1, '.sam')
meta_data$hisatLog <- paste('/Users/mtaillefer00/Documents/BIOL_4315_Lab_4_MT/bam_files/',meta_data$SampleName, sep = '')
meta_data$hisatLog <- paste0(meta_data$hisatLog, '_hisat2.log')
meta_data$bam1 <- paste('/Users/mtaillefer00/Documents/BIOL_4315_Lab_4_MT/bam_files/',meta_data$SampleName, sep = '')
meta_data$bam1 <- paste0(meta_data$bam1, '.bam')
meta_data$bam1sorted <- paste('/Users/mtaillefer00/Documents/BIOL_4315_Lab_4_MT/bam_files/',meta_data$SampleName, sep = '')
meta_data$bam1sorted <- paste0(meta_data$bam1sorted, '.sorted.bam')

for( i in 1:nrow(meta_data)){
  hisat2 <- tryCatch({system2("hisat2", 
                    args = c('-x', 'hisat2_index/tair10_1_index', 
                             '-U',meta_data$processedUnzip1[i],
                             '-p', '8',
                             '-S', meta_data$sam1[i]),
                    stdout = TRUE, stderr = TRUE)}, 
                    error = function(e) {
                      paste("hisat2", "indexing failed with error:", e$message)
                      })
  
 tryCatch({system2('samtools', 
                           args = c('view', '-bS', meta_data$sam1[i], '-o', meta_data$bam1[i]),
                   stdout = TRUE, stderr = TRUE)}, 
                   error = function(e) {
                             paste("samtools", "sam>bam failed with error:", e$message)
                           })
  
  tryCatch({system2('samtools', 
                           args = c('sort', meta_data$bam1[i], '-o', meta_data$bam1sorted[i]),
                    stdout = TRUE, stderr = TRUE)}, 
                    error = function(e) {
                      paste("samtools", "sorting failed with error:", e$message)
                    })
  
  tryCatch({system2('samtools', 
                           args = c('index', meta_data$bam1sorted[i]),
                    stdout = TRUE, stderr = TRUE)}, 
                    error = function(e) {
                      paste("samtools", "indexing failed with error:", e$message)
                  })
  writeLines(hisat2,meta_data$hisatLog[i])
}

#5.3.4 Reading alignment stats

# get the folder
hisat2_logs_dir <- "bam_files/"

# 1. Get a list of all HISAT2 log files
log_files <- list.files(hisat2_logs_dir, pattern = "*.log", full.names = TRUE)

#preper an "empty" vector of the apropreate size, faster than appending. 
percent_aligned <- 1:length(log_files)

#loop through log files
for (i in seq_along(percent_aligned)) {
  
  percent_aligned[i] <- readLines(log_files[i])[length(readLines(log_files[i]))]
  
}

#bind vectors as dataframe
align_df <- data.frame(sort(meta_data$SampleName),percent_aligned)

# extract the numeric percent value
align_df <- align_df %>% 
  mutate(percent_aligned = as.numeric(
    stringr::str_split_i(align_df$percent_aligned, "%",1))) %>% 
  rename(samplename = sort.meta_data.SampleName.)

head(align_df)

#question 5
#display align_df
DT::datatable(align_df)

#plot the data
ggplot(align_df, aes(x = "", y = percent_aligned)) +  
geom_boxplot() +
  xlab("") +  
  ylab("percent aligned") +
  ggtitle("Boxplot of alignment")

#list bam files
bfiles <- list.files("./bam_files", pattern = ".sorted.bam$", full.names = TRUE)

#Counting how many reads corespond to each gene
gene_count_list <- Rsubread::featureCounts(files = bfiles, annot.ext = "data/GCF_000001735.4_TAIR10.1_genomic.gtf", 
                                           isGTFAnnotationFile = TRUE, # <--- input annotation is GTF
                                           allowMultiOverlap = FALSE,  # <--- don't allow reads that overlap with multiple loci
                                           isPairedEnd = FALSE, nthreads = 8,
                                           minMQS = 10, # <--- minimum mapping quality score of 10 (like a phred score for the hisat2 alignment)
                                           GTF.featureType = "exon",  # <--- Count reads overlapping 'exon' features
                                           GTF.attrType = "gene_id" # <--- Groups exon by gene id
)

glimpse(gene_count_list$annotation)[1:5,]
glimpse(gene_count_list$stat)[,1:5]


count_table <- gene_count_list$counts

#drop .sorted.bam from column titles
colnames(count_table) <- gsub(".sorted\\.bam$", "", colnames(count_table))

#filter out genes not mapped to any reads
filtered_table <- count_table[rowSums(count_table) != 0, ]

#interactable table
DT::datatable(filtered_table)

#5.5.1 Data prep for DESeq2

#Get the metadata table that would accompany the count table
coldata <- meta_data %>% dplyr::select(SampleName,SampleLong,Factor) %>% 
  dplyr::mutate(SampleLong=str_split_i(SampleLong, "\\.",1)) %>% #getting the groups name (Avirulent, Mock and Virulent)
  dplyr::rename(condition = SampleLong) %>%
  dplyr::mutate(condition = factor(condition)) %>% #for group comp
  dplyr::mutate(Factor = factor(Factor)) #for sample comp

base::rownames(coldata) <- coldata$SampleName
coldata <- coldata %>% mutate(SampleName = factor(SampleName))

coldata$type <- factor(rep("single-read", nrow(coldata)))

#Make sure samples are in the same order between the two. 
coldata <- coldata[base::match(base::colnames(count_table), rownames(coldata)),]

#Check that they are indeed same order
all(rownames(coldata) == base::colnames(count_table))

#creating a dds object where the conditions are avr vs mock vs vir

dds1 <- DESeqDataSetFromMatrix(countData = count_table,
                               colData = coldata,
                               design = ~ condition)
#creating a dds object where the conditions are the samples themslevs 
dds2 <- DESeqDataSetFromMatrix(countData = count_table,
                               colData = coldata,
                               design = ~ Factor)

#5.5.2 Sample correlation based on transcript abundance table

#correlating the samples
d <- cor(assay(rlog(dds1)), method = "spearman")
#turning correlation to a distance (1 - correlation) and clustring
hc <- hclust(dist(1 - d))

#hierarchal clustering
plot.phylo(as.phylo(hc), type = "p", edge.col = "blue", edge.width = 2,
           show.node.label = TRUE, no.margin = TRUE)

#5.5.3 Analyzing differential gene expression (DGE) with DESeq2

dds1_results <- DESeq(dds1)
dds2_results <- DESeq(dds2)
res1 <- DESeq2::results(dds1_results)
res1
res2 <- DESeq2::results(dds2_results)
res2
dds2

#5.5.3.1 Comparing Vir, Mock, and Avr broad overview
res_vir_mock <- DESeq2::results(dds1_results, contrast = c("condition", "Vir", "Mock"), alpha = 0.2)

# Avr vs Mock  
res_avr_mock <- DESeq2::results(dds1_results, contrast = c("condition", "Avr", "Mock"), alpha = 0.2)

# Vir vs Avr
res_vir_avr <- DESeq2::results(dds1_results, contrast = c("condition", "Vir", "Avr"), alpha = 0.2)

# Function to filter and count DE genes
filter_and_count <- function(res_obj, comparison_name, fc_threshold = 2) {
  # Remove NAs
  res_filtered <- res_obj[!is.na(res_obj$padj) & !is.na(res_obj$log2FoldChange), ]
  
  # Apply filters: |log2FC| >= log2(2) = 1 and padj <= alpha (already set in results())
  sig_genes <- res_filtered[abs(res_filtered$log2FoldChange) >= log2(fc_threshold), ]
  
  # Count up and down regulated
  up_regulated <- sum(sig_genes$log2FoldChange > 0)
  down_regulated <- sum(sig_genes$log2FoldChange < 0)
  
  return(data.frame(
    Comparison = comparison_name,
    Up_regulated = up_regulated,
    Down_regulated = down_regulated
  ))
}

# Apply filtering and counting to all comparisons
results_summary <- rbind(
  filter_and_count(res_vir_mock, "Vir vs Mock"),
  filter_and_count(res_avr_mock, "Avr vs Mock"),
  filter_and_count(res_vir_avr, "Vir vs Avr")
)

# Print summary
print("Summary of DE genes (FC >= 2, alpha = 0.2):")
print(results_summary)

plot_data <- results_summary %>%
  pivot_longer(cols = c(Up_regulated, Down_regulated), 
               names_to = "Regulation", 
               values_to = "Count") %>%
  mutate(Regulation = factor(Regulation, levels = c("Up_regulated", "Down_regulated")))

# Create horizontal stacked bar plot
p <- ggplot(plot_data, aes(x = Comparison, y = Count, fill = Regulation)) +
  geom_bar(stat = "identity", position = "stack") +
  coord_flip() +  # Makes it horizontal
  labs(
    title = "Differentially Expressed Genes by Comparison",
    subtitle = "Fold Change >= 2, alpha = 0.2",
    x = "Comparison",
    y = "Number of Genes",
    fill = "Regulation"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10))
p

# Quesiont 7
comp <- systemPipeR::readComp("/Users/mtaillefer00/Library/R/arm64/4.4/library/systemPipeRdata/extdata/param/targetsPE.txt")
comp[[1]]
compPair <- comp[[1]]

results <- list()
for (pair in compPair){
  groups <- unlist(strsplit(pair,'-'))
  gr1 <- groups[1]
  gr2 <- groups[2]
  compName <- paste(gr1, 'vs', gr2)
  filtered_data <- filter_and_count(res, compName)
  results[[compName]] <- filtered_data
}

results_summary2 <- do.call(rbind, results)

plot_data2 <- results_summary2 %>%
  pivot_longer(cols = c(Up_regulated, Down_regulated),
               names_to = "Regulation",
               values_to = "Count") %>%
  mutate(Regulation = factor(Regulation, levels = c("Up_regulated", "Down_regulated")))

p_2 <- ggplot(plot_data2, aes(x = Comparison, y = Count, fill = Regulation)) +
  geom_bar(stat = "identity", position = "stack") +
  coord_flip() +
  labs(
    title = "Differentially Expressed Genes by Sample Comparison dd2",
    subtitle = "Fold Change >= 2,, alpha = 0.2",
    x = "Comparison",
    y = "Number of Genes",
    fill = "Regulation"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10))

p_2

#5.5.3.3 Adding gene descriptions and getting specific with volcano plots

m <- biomaRt::useMart("plants_mart", dataset = "athaliana_eg_gene",
                      host = "https://plants.ensembl.org")
desc <- biomaRt::getBM(attributes = c("tair_locus", "description"), mart = m)
desc <- desc[!duplicated(desc[, 1]), ]
desc <- desc %>% rename( gene_id = tair_locus)

annotate_results <- function(res_obj, desc_df) {
  res_df <- as.data.frame(res_obj)
  res_df$gene_id <- rownames(res_df)
  res_df <- left_join(res_df, desc_df, by = "gene_id") %>%
    mutate(description = str_split_i(description,"\\[",1)) # <-- remove the source of the annotation from the description
  
  return(res_df)
}

# Add annotations to all results
res_vir_mock_annot <- annotate_results(res_vir_mock, desc)
res_avr_mock_annot <- annotate_results(res_avr_mock, desc)
res_vir_avr_annot <- annotate_results(res_vir_avr, desc)


volcano1 <- EnhancedVolcano(res_vir_mock_annot,
                            lab = res_vir_mock_annot$description,
                            x = 'log2FoldChange',
                            y = 'pvalue',
                            title = 'Vir vs Mock',
                            pCutoff = 0.05,           # pvalue threshold
                            FCcutoff = 1.0,
                            pointSize = 4.0,
                            labSize = 4.0,
                            labCol = 'black',
                            labFace = 'bold',
                            boxedLabels = TRUE,
                            colAlpha = 4/5,
                            legendPosition = 'right',
                            legendLabSize = 14,
                            legendIconSize = 4.0,
                            drawConnectors = TRUE,
                            widthConnectors = 1.0,
                            colConnectors = 'black') + ggplot2::scale_y_continuous(
                              breaks=seq(0,4, 1))

volcano1

volcano2 <- EnhancedVolcano(res_vir_avr_annot,
                                lab = res_vir_avr_annot$description,
                                x = 'log2FoldChange',
                                y = 'pvalue',
                                title = 'Vir vs avir',
                                pCutoff = 0.05,           # pvalue threshold
                                FCcutoff = 1.0,
                                pointSize = 4.0,
                                labSize = 4.0,
                                labCol = 'black',
                                labFace = 'bold',
                                boxedLabels = TRUE,
                                colAlpha = 4/5,
                                legendPosition = 'right',
                                legendLabSize = 14,
                                legendIconSize = 4.0,
                                drawConnectors = TRUE,
                                widthConnectors = 1.0,
                                colConnectors = 'black')  + ggplot2::scale_y_continuous(
                                  breaks=seq(0,7, 1))

volcano2

volcano3 <- EnhancedVolcano(res_avr_mock_annot,
                            lab = res_avr_mock_annot$description,
                            x = 'log2FoldChange',
                            y = 'pvalue',
                            title = 'Avr vs Mock',
                            pCutoff = 0.05,           # pvalue threshold
                            FCcutoff = 1.0,
                            pointSize = 4.0,
                            labSize = 4.0,
                            labCol = 'black',
                            labFace = 'bold',
                            boxedLabels = TRUE,
                            colAlpha = 4/5,
                            legendPosition = 'right',
                            legendLabSize = 14,
                            legendIconSize = 4.0,
                            drawConnectors = TRUE,
                            widthConnectors = 1.0,
                            colConnectors = 'black') + ggplot2::scale_y_continuous(
                              breaks=seq(0,4, 1))

volcano3

#5.5.4 Clustered heatmap based on significant DEGs
