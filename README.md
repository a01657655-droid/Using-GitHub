# RNA-Seq Differential Expression Analysis

## Overview
This repository contains an RNA-Seq analysis pipeline using edgeR and limma packages in R to identify differentially expressed genes in Saccharomyces cerevisiae under different conditions.

## Author
Joana Paulina Pineda Corella (A01657655)

## Analysis Workflow
1. Data loading and preprocessing
2. Quality control and normalization
3. Differential expression analysis using edgeR
4. Visualization of results
5. Generalized Linear Model (GLM) analysis for complex comparisons

# 1. Installation and Setup
```r
# Install Bioconductor manager and required packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("edgeR")  # For differential expression analysis of RNA-seq data
BiocManager::install("limma")   # For linear models and empirical Bayes methods

library(limma)
library(edgeR)
```

- BiocManager is used to install and manage Bioconductor packages
- edgeR is specifically designed for differential expression analysis of RNA-seq count data
- limma provides powerful linear modeling capabilities often used in conjunction with edgeR

# 2. Project Directory Setup
```r
# Create project directory structure
dir.create("RNAseq_project", showWarnings = FALSE)
setwd("RNAseq_project")
dir.create("results_diffexp", showWarnings = FALSE)
```
- Organized directory structure ensures reproducible analysis
- Separate results folder keeps output organized and easily accessible

# 3. Data Loading and Inspection 
```r
# Read count data from remote source
counts <- read.table("https://raw.githubusercontent.com/jmvillalobos/RNAseq_curso2025/main/Saccharomyces.txt",
                     header = TRUE, row.names = 1, sep = "\t")

head(counts)  # Preview first 6 rows
```
<img width="626" height="275" alt="Screenshot 2025-09-04 at 10 22 26 p m" src="https://github.com/user-attachments/assets/102725ca-5791-44af-9360-235271ed3868" />

```r
dim(counts)   # Check dimensions (genes × samples)
```
<img width="169" height="32" alt="Screenshot 2025-09-04 at 10 22 49 p m" src="https://github.com/user-attachments/assets/6d1aae73-a8dd-4ebc-8129-074b6dfe2330" />

- read.table imports tab-separated data with proper formatting
- header = TRUE uses first row as column names
- row.names = 1 uses first column as row names (gene IDs)
- sep = "\t" specifies tab-separated format

# 4. Data Filtering 
```r
# Filter low-expression genes using CPM threshold
counts = counts[rowSums(cpm(counts) >= 2) >=10,]
head(counts)
```
<img width="629" height="275" alt="Screenshot 2025-09-04 at 10 24 20 p m" src="https://github.com/user-attachments/assets/19121d70-5b04-4c9d-acd4-d72aa6eea701" />

- cpm() converts raw counts to counts per million (normalizes for library size)
- Filter keeps genes with at least 2 CPM in at least 10 samples
- Removes genes with low expression that could affect statistical analysis

# 5. Data Exploration and Quality Control
```r
colnames(counts) #nombres de las columnas (muestras) para verificar que no se modificaron 
```
<img width="669" height="51" alt="Screenshot 2025-09-04 at 10 25 35 p m" src="https://github.com/user-attachments/assets/0fd1a32b-b8dc-4aca-94cc-f09c9ad6841a" />

```r
dim(counts)
```
<img width="174" height="34" alt="Screenshot 2025-09-04 at 10 26 03 p m" src="https://github.com/user-attachments/assets/ab287a8c-eeb8-4a58-8917-7a0ec806ca7a" />

```r
# Scatter plot to assess technical variability
plot(log2(counts[,c("wt_sc_1", "st_sc_1")]), col="violet")
```
<img width="659" height="416" alt="Screenshot 2025-09-04 at 10 27 01 p m" src="https://github.com/user-attachments/assets/fee8703a-b1f3-4076-8bd2-deca1e2db282" />

```r
# Create experimental groups from sample names
grp = sub("..$", "", colnames(counts))
grp
```
<img width="664" height="53" alt="Screenshot 2025-09-04 at 10 27 40 p m" src="https://github.com/user-attachments/assets/a98154f8-889a-44fb-b170-abb09c6aaa18" />

- Log transformation stabilizes variance for visualization
- Group creation using sub() removes replicate identifiers (last two characters)
- This creates meaningful biological groups for comparison

# 6. edgeR Object Creation and Normalization
```r
# Create DGEList object for edgeR analysis
dge = DGEList(counts=counts, group=grp)
```
<img width="625" height="459" alt="Screenshot 2025-09-04 at 10 28 53 p m" src="https://github.com/user-attachments/assets/4ee08554-e139-42f0-a329-65a65c6ed725" />

```r
# Multidimensional scaling for quality control
plotMDS(dge, col="violet")
```
<img width="658" height="417" alt="Screenshot 2025-09-04 at 10 29 41 p m" src="https://github.com/user-attachments/assets/7a3b9a01-36d6-4cf8-a534-e44ad7fb4ae7" />

```r
# Normalize using TMM method
dgeNorm = calcNormFactors(dge)
dgeNorm$samples
```
<img width="346" height="257" alt="Screenshot 2025-09-04 at 10 30 55 p m" src="https://github.com/user-attachments/assets/9ff5cb4f-dc83-4c1f-939b-a7f3163daa71" />

- DGEList creates an object that stores counts and sample information
- plotMDS visualizes sample similarities and identifies outliers
- calcNormFactors performs TMM normalization to correct for library size differences

# 7. Dispersion Estimation
```r
# Estimate common dispersion across all genes
dgeNorm = estimateCommonDisp(dgeNorm)
dgeNorm$common.dispersion
```
<img width="171" height="32" alt="Screenshot 2025-09-04 at 10 32 49 p m" src="https://github.com/user-attachments/assets/bc5a1a38-38b1-4b93-bf87-38db0d79574c" />

- Dispersion estimation models biological variability
- Essential for accurate statistical testing in count-based data

# 8. Differential Expression Analysis
```r
# Perform exact test for pairwise comparisons
diff_exp = exactTest(dgeNorm, dispersion = dgeNorm$common.dispersion, pair = c("wt_sc", "wt_sl"))
diff_exp2 = exactTest(dgeNorm, dispersion = dgeNorm$common.dispersion, pair = c("st_sc", "st_sl"))
diff_exp
```
<img width="426" height="293" alt="Screenshot 2025-09-04 at 10 35 12 p m" src="https://github.com/user-attachments/assets/48717ecc-55f9-4c7e-a512-def6041470fa" />

```r
# Extract and examine results
?exactTest
dim(diff_exp)
```
<img width="187" height="36" alt="Screenshot 2025-09-04 at 10 35 51 p m" src="https://github.com/user-attachments/assets/3c62929a-2d36-4ce0-908e-d3108c7a6d2b" />

- exactTest performs pairwise differential expression testing
- Uses negative binomial distribution appropriate for count data

```r
topTags(diff_exp)
```
<img width="490" height="238" alt="Screenshot 2025-09-04 at 10 36 55 p m" src="https://github.com/user-attachments/assets/a4ff5eec-dde6-40f8-8d10-5d0a53d10d68" />

```r
dim(topTags(diff_exp))
```
<img width="139" height="34" alt="Screenshot 2025-09-04 at 10 38 01 p m" src="https://github.com/user-attachments/assets/fc382940-1eac-4fcb-bc54-bfa9c25d7901" />

```r
deTab = topTags(diff_exp, n=Inf)$table
deTab[c(15,30),]
```

<img width="488" height="70" alt="Screenshot 2025-09-04 at 10 38 50 p m" src="https://github.com/user-attachments/assets/49c8cb56-b404-419f-bd53-488851974a0f" />

```r
row.names(deTab)[deTab$logFC > 5]
```

<img width="627" height="106" alt="Screenshot 2025-09-04 at 10 40 25 p m" src="https://github.com/user-attachments/assets/f9d0d939-effc-4e08-af60-b288f2b094f1" />

```r
deTab["YNL284C-A",]
```

<img width="464" height="49" alt="Screenshot 2025-09-04 at 10 41 05 p m" src="https://github.com/user-attachments/assets/2e220156-fdc7-41c9-8eab-02999ac4e036" />

- topTags extracts and sorts results by statistical significance

# 9. Result Filtering and Interpretation
```r
# Filter results using FDR < 0.05 and |logFC| > 1 thresholds
deGenes = rownames(deTab)[deTab$FDR < 0.05 & abs(deTab$logFC) > 1]
down = row.names(deTab)[deTab$logFC < -1]  # Downregulated genes
over = row.names(deTab)[deTab$logFC > 1]   # Upregulated genes
```

```r
print(paste("total de diferenciales:", length(deGenes)))
```

<img width="329" height="35" alt="Screenshot 2025-09-04 at 10 42 32 p m" src="https://github.com/user-attachments/assets/0bf35513-f8bf-41ab-9c0c-9e17e64d7998" />

```r
print(paste("número de genes inducidos:", length(over)))
```

<img width="368" height="34" alt="Screenshot 2025-09-04 at 10 43 04 p m" src="https://github.com/user-attachments/assets/2b451d02-344c-48c9-abca-cfff47729d5c" />

```r
print(paste("número de genes reprimidos:", length(down)))
```

<img width="362" height="37" alt="Screenshot 2025-09-04 at 10 43 27 p m" src="https://github.com/user-attachments/assets/dff9ecc0-57f9-4989-ad2d-176391741ccc" />

- False Discovery Rate (FDR) < 0.05 controls for multiple testing
- |logFC| > 1 filters for biologically relevant changes (2-fold change)
- Separate identification of up- and down-regulated genes

# 10. Visualization

```r
# Create diagnostic plots
plotSmear(dge, de.tags=deGenes, ylab = "WT-sc vs WT-sl")
```
<img width="667" height="417" alt="Screenshot 2025-09-04 at 10 44 59 p m" src="https://github.com/user-attachments/assets/bc1a7f80-1ead-4255-b194-58cbb642710f" />

```r
# Heatmap of differentially expressed genes
normalizados = cpm(counts)
normalizados_diferenciales = normalizados[deGenes,]
head(normalizados_diferenciales)
```
<img width="658" height="275" alt="Screenshot 2025-09-04 at 10 46 12 p m" src="https://github.com/user-attachments/assets/e7ef8fda-6c13-449d-a1ce-5eec5b79db70" />

```r
heatmap(normalizados_diferenciales, col=my_colors, scale="row", margins=c(5,10))
```
<img width="418" height="482" alt="Screenshot 2025-09-04 at 10 47 22 p m" src="https://github.com/user-attachments/assets/4854f3b8-d4cc-4490-99c6-e055d3b008d6" />

- plotSmear shows log-fold changes versus average expression
- Heatmap visualizes expression patterns across samples
- scale="row" z-score normalizes for better visualization

# 11. Advanced Analysis with GLM Framework

```r
par(mfrow=c(1,2)) 
boxplot(log(counts),col=rainbow(6), main="antes de la normalización")
boxplot(log(normalizados), col=rainbow(6), main="después de la normalización")
```

<img width="625" height="433" alt="Screenshot 2025-09-04 at 10 48 51 p m" src="https://github.com/user-attachments/assets/ce5ed6ee-63e6-4012-abb7-2998658aa218" />

```r
barplot(apply(normalizados_diferenciales,2,sum),las=2, cex.names = 1, col = (1:6))
```
<img width="625" height="394" alt="Screenshot 2025-09-04 at 10 49 18 p m" src="https://github.com/user-attachments/assets/7c146f38-ce07-4432-b70b-9df696b4be52" />

```r
pca <- princomp(normalizados_diferenciales[,c(1:6)])
plot(pca$loadings, col=as.factor(colnames(normalizados_diferenciales[,c(1:6)])),  pch=19, cex=2, main="con nitrógeno")
text(pca$loadings, as.vector(colnames(normalizados_diferenciales[,c(1:6)])), pos=3, cex=0.8)
```
<img width="665" height="449" alt="Screenshot 2025-09-04 at 10 49 50 p m" src="https://github.com/user-attachments/assets/31357e8a-dfab-4ec3-b49b-25660d66719f" />

```r
pca <- princomp(normalizados_diferenciales[,c(7:12)])
plot(pca$loadings, col=as.factor(colnames(normalizados_diferenciales[,c(7:12)])),  pch=19, cex=2, main="si nitrógeno")
text(pca$loadings, as.vector(colnames(normalizados_diferenciales[,c(7:12)])), pos=3, cex=0.8)
```

<img width="676" height="472" alt="Screenshot 2025-09-04 at 10 50 22 p m" src="https://github.com/user-attachments/assets/bb9c561e-b5e9-4b3e-a640-d96f1aaead4f" />

```r
with(deTab, plot(logFC, -log10(FDR), pch=20, cex=0.8, col="purple", main="WT+N vs WT-N", xlim=c(-8, 8), ylim=c(0,300)))
text(deTab[1:20,]$logFC,-log(deTab[1:20,]$FDR,10),labels=rownames(deTab[1:20,]),cex=0.7,pos=1)
with(subset(deTab, FDR<.01 & abs(logFC)>2), points(logFC, -log10(FDR), pch=20, cex=0.5, col="pink"))
abline(v=2,lty=2, col="blue")
abline(v=-2,lty=2, col="blue")
legend("bottomright","Up_regulated",cex=1)
legend("bottomleft","Down_regulated",cex=1)
```

<img width="665" height="456" alt="Screenshot 2025-09-04 at 10 50 49 p m" src="https://github.com/user-attachments/assets/b7fa85bf-85df-479d-859f-9ac1dda66dc5" />

```r
deTab2 = topTags(diff_exp2, n=Inf)$table
deGenes2 = rownames(deTab2)[deTab2$FDR < 0.05 & abs(deTab2$logFC) > 1]

with(deTab2, plot(logFC, -log10(FDR), pch=20, cex=0.8, col="purple", main="ste12+N vs ste12-N", xlim=c(-10, 10), ylim=c(0,320)))
text(deTab2[1:20,]$logFC,-log(deTab2[1:20,]$FDR,10),labels=rownames(deTab2[1:20,]),cex=0.7,pos=1)
with(subset(deTab2, FDR<.01 & abs(logFC)>2), points(logFC, -log10(FDR), pch=20, cex=0.5, col="pink"))
abline(v=2,lty=2, col="blue")
abline(v=-2,lty=2, col="blue")
legend("bottomright","Up_regulated",cex=1)
legend("bottomleft","Down_regulated",cex=1)
```
<img width="662" height="450" alt="Screenshot 2025-09-04 at 10 51 14 p m" src="https://github.com/user-attachments/assets/3a1ef681-a30f-4c71-85e2-39e905908535" />

```r
WTover= head(rownames(deTab), 30)
ste12over= head(rownames(deTab2), 30)
setdiff(WTover, ste12over)
```
<img width="424" height="33" alt="Screenshot 2025-09-04 at 10 51 53 p m" src="https://github.com/user-attachments/assets/4db1d0cc-f1be-40d1-bc2a-019a9cacefe8" />

# RNA-Seq Analysis with Generalized Linear Models (GLM)

## Overview
This repository contains an advanced RNA-Seq differential expression analysis pipeline using edgeR's Generalized Linear Model (GLM) framework. The analysis examines Saccharomyces cerevisiae gene expression under different nitrogen conditions and genetic backgrounds (wild-type vs. ste12 mutant).

```r
# Install Bioconductor manager if not already present
if (!requireNamespace("BiocManager", quietly = TRUE)) 
    install.packages("BiocManager")

# Install required Bioconductor packages
for (pkg in c("edgeR","limma")) {
  if (!requireNamespace(pkg, quietly = TRUE)) 
      BiocManager::install(pkg, ask = FALSE, update = TRUE)
}

# Load necessary libraries
library(limma)
library(edgeR)
```

- Uses BiocManager for controlled installation of Bioconductor packages
- Installs edgeR for differential expression analysis with GLM framework
- Installs limma for linear models and contrast matrix functionality
- ask = FALSE and update = TRUE ensure smooth, non-interactive installation

# Project Structure Setup
```r
# Create dedicated directory for results
outpathcount <- "results_edger_yeast"
dir.create(outpathcount, showWarnings = FALSE, recursive = TRUE)
```

- recursive = TRUE allows creating nested directories if needed
- Organized output structure ensures reproducible results management

# Data Loading and Inspection
```r
# Read count data from remote repository
counts <- read.table(
  "https://raw.githubusercontent.com/jmvillalobos/RNAseq_curso2025/main/Saccharomyces.txt",
  header = TRUE, row.names = 1, sep = "\t", comment.char = ""
)

# Check data dimensions
dim(counts)
```

<img width="167" height="35" alt="Screenshot 2025-09-04 at 11 00 44 p m" src="https://github.com/user-attachments/assets/a2fe3b80-e705-4400-900e-b25d2d5e4400" />

```r
# Preview the data structure
head(counts[,1:4])
```

<img width="382" height="147" alt="Screenshot 2025-09-04 at 11 01 13 p m" src="https://github.com/user-attachments/assets/7faf259f-f42a-4d82-93b4-c14738f139da" />

- Reads tab-separated count data with proper formatting options
- comment.char = "" ensures no lines are treated as comments
- Dimension check verifies data integrity (genes × samples)
- Preview shows the structure of count data

# Experimental Design Setup

```r
# Create experimental groups from sample names
grp_names <- sub("..$", "", colnames(counts))
grp <- factor(grp_names, levels = c("wt_sc","wt_sl","st_sc","st_sl"))
table(grp)
```

<img width="246" height="69" alt="Screenshot 2025-09-04 at 11 02 15 p m" src="https://github.com/user-attachments/assets/8352058c-cf33-4ff7-b3b9-74151767fe23" />

- sub("..$", "", colnames(counts)) removes the last two characters (replicate identifiers)
- Converts to factor with explicit level ordering for controlled comparisons
- table(grp) provides sample distribution across conditions

# Data Preprocessing 
```r
# Create DGEList object for edgeR analysis
dge <- DGEList(counts = counts, group = grp)

# Filter low-expression genes using filterByExpr
keep <- filterByExpr(dge, group = grp)
dge <- dge[keep, , keep.lib.sizes = FALSE]
summary(keep); dim(dge)
```

<img width="278" height="48" alt="Screenshot 2025-09-04 at 11 03 08 p m" src="https://github.com/user-attachments/assets/4e0db709-36d6-4624-929d-068b083669c4" />
<img width="205" height="36" alt="Screenshot 2025-09-04 at 11 03 17 p m" src="https://github.com/user-attachments/assets/e5ad1c6f-ac78-4dff-beaf-5e6265d97bec" />

- DGEList object stores counts and sample metadata
- filterByExpr uses statistical methods to identify genes with sufficient expression
- More sophisticated than simple CPM-based filtering
- keep.lib.sizes = FALSE recalculates library sizes after filtering

# Normalization
```r
# Perform TMM normalization
dge <- calcNormFactors(dge, method = "TMM")
dge$samples
```

<img width="358" height="256" alt="Screenshot 2025-09-04 at 11 04 17 p m" src="https://github.com/user-attachments/assets/7ba12a95-b88a-45db-acf5-6efd1f41fd2e" />

```r
# Check dimensions post-filtering
dim(dge)
```

<img width="171" height="34" alt="Screenshot 2025-09-04 at 11 05 01 p m" src="https://github.com/user-attachments/assets/ec9f6b85-28f8-48d6-92b8-12b561b71e18" />

- TMM (Trimmed Mean of M-values) normalization corrects for:
    - Library size differences
    - RNA composition biases
- Normalization factors are added to the samples data frame

# Quality Control: MDS Plots
```r
# Define color scheme for experimental groups
my_colors <- c("wt_sc"="purple", "wt_sl"="pink", 
               "st_sc"="violet", "st_sl"="magenta")

# MDS plot using top 500 most variable genes
plotMDS(dge, top = 500, labels = colnames(dge),
        col = my_colors, main = "MDS: separación por grupo (top=500 genes)")
```

<img width="669" height="454" alt="Screenshot 2025-09-04 at 11 06 16 p m" src="https://github.com/user-attachments/assets/ef41e932-5923-45f3-87ba-e49d78377a1f" />

```r
# MDS plot using all genes
plotMDS(dge, top = nrow(dge), labels = colnames(dge), 
        col = my_colors, main = "MDS: separación por grupo (todos los genes)")
```

<img width="657" height="456" alt="Screenshot 2025-09-04 at 11 06 54 p m" src="https://github.com/user-attachments/assets/99a56c9d-688e-4235-879d-f4d06f719f55" />

- MDS (Multidimensional Scaling) plots visualize sample-to-sample distances
- Using top 500 genes focuses on most biologically relevant variation
- Using all genes provides comprehensive overview of technical and biological variation
- Color-coding by experimental condition enables visual assessment of group separation

```r
# Create design matrix for GLM analysis
design <- model.matrix(~ 0 + dge$samples$group)
colnames(design) <- levels(dge$samples$group)
design
```

<img width="381" height="350" alt="Screenshot 2025-09-04 at 11 07 37 p m" src="https://github.com/user-attachments/assets/62cf112b-9d55-4ea0-bc93-9f4f72ef7826" />

- model.matrix(~ 0 + group) creates a design matrix without intercept
- Each column represents one experimental condition
- Essential for specifying complex comparisons in the GLM framework

# Dispersion Estimation
```r
# Three-level dispersion estimation
dge <- estimateGLMCommonDisp(dge, design = design)  # Common dispersion
dge <- estimateGLMTrendedDisp(dge, design = design) # Trended dispersion
dge <- estimateGLMTagwiseDisp(dge, design = design) # Gene-specific dispersion

print(dge$common.dispersion)
```

<img width="213" height="35" alt="Screenshot 2025-09-04 at 11 08 33 p m" src="https://github.com/user-attachments/assets/ef921b70-74f5-4d41-bd0d-3aba498e4502" />

- Common dispersion: Average dispersion across all genes
- Trended dispersion: Dispersion as a function of gene expression level
- Tagwise dispersion: Gene-specific dispersion estimates
- This three-level approach provides robust variance estimation

# GLM Model Fitting
```r
# Fit GLM model to the data
fit <- glmFit(dge, design = design)

# Define biological comparisons of interest
contVector <- c(
  "WT_Nneg_vs_WT_Npos" = "wt_sl - wt_sc",      # WT: -N vs +N
  "ST_Nneg_vs_ST_Npos" = "st_sl - st_sc",      # ste12: -N vs +N
  "ST_vs_WT_at_Npos"   = "st_sc - wt_sc",      # ste12 vs WT at +N
  "ST_vs_WT_at_Nneg"   = "st_sl - wt_sl"       # ste12 vs WT at -N
)

# Convert to contrast matrix
contMatrix <- makeContrasts(contrasts = contVector, levels = design)
colnames(contMatrix) <- names(contVector)
contMatrix
```

<img width="691" height="126" alt="Screenshot 2025-09-04 at 11 09 43 p m" src="https://github.com/user-attachments/assets/12775bf2-8d4b-4a95-b19f-031d9c542c0e" />

- Fits negative binomial generalized linear model to count data
- Model incorporates design matrix and dispersion estimates
- Provides foundation for statistical testing of specific contrasts
- Defines four biologically meaningful comparisons:
    1. Nitrogen response in wild-type
    2. Nitrogen response in ste12 mutant
    3. Genetic difference under nitrogen sufficiency
    4. Genetic difference under nitrogen limitation
- makeContrasts converts verbal descriptions to mathematical contrasts

# Differential Expression Analysis
```r
# Set statistical thresholds
alpha <- 0.05   # FDR cutoff
minLFC <- 1     # Minimum log2 fold change

# Analyze each contrast
for (comp in colnames(contMatrix)) {
  cat("===>", comp, "\n")
  
  # Perform likelihood ratio test
  lrt <- glmLRT(fit, contrast = contMatrix[, comp])
  topTab <- topTags(lrt, n = Inf)$table
  
  # Identify significant genes
  deGenes <- rownames(topTab)[topTab$FDR < alpha & abs(topTab$logFC) > minLFC]
  deGenes_up <- rownames(topTab)[topTab$FDR < alpha & topTab$logFC > minLFC]
  deGenes_down <- rownames(topTab)[topTab$FDR < alpha & topTab$logFC < -minLFC]
  
  # Output results summary
  print(length(deGenes)); print(length(deGenes_up)); print(length(deGenes_down))
  
  # Create diagnostic plot
  pdf(file.path(outpathcount, paste0("Smear_", comp, ".pdf")))
  plotSmear(lrt, de.tags = deGenes, ylab = paste(comp, "logFC"))
  abline(h = c(1, -1), col = "violet", lty = 2)
  dev.off()
  
  # Save results to files
  write.table(topTab, file = file.path(outpathcount, paste0(comp, "_table.txt")),
              row.names = TRUE, col.names = NA, quote = FALSE, sep = "\t")
  
  write.table(deGenes_up, file = file.path(outpathcount, paste0(comp, "_up.txt")),
              row.names = TRUE, col.names = NA, quote = FALSE, sep = "\t")
  
  write.table(deGenes_down, file = file.path(outpathcount, paste0(comp, "_down.txt")),
              row.names = TRUE, col.names = NA, quote = FALSE, sep = "\t")
}
```

<img width="596" height="88" alt="Screenshot 2025-09-04 at 11 14 38 p m" src="https://github.com/user-attachments/assets/45e3c67d-4345-4e26-a7af-c3c047b532a7" />
<img width="725" height="52" alt="Screenshot 2025-09-04 at 11 14 55 p m" src="https://github.com/user-attachments/assets/55d7b6a6-fbcb-4dba-8d79-049bc44d8a7b" />
<img width="624" height="90" alt="Screenshot 2025-09-04 at 11 15 13 p m" src="https://github.com/user-attachments/assets/d01ee47c-23a7-4ea6-8f40-a818532cca3b" />
<img width="692" height="51" alt="Screenshot 2025-09-04 at 11 15 26 p m" src="https://github.com/user-attachments/assets/7b768d2a-d564-4d01-9c3b-b8a9aa6a2e8e" />
<img width="575" height="87" alt="Screenshot 2025-09-04 at 11 15 40 p m" src="https://github.com/user-attachments/assets/ae078bf8-fc72-4243-b621-5d452f1fa581" />
<img width="757" height="50" alt="Screenshot 2025-09-04 at 11 15 57 p m" src="https://github.com/user-attachments/assets/283c6e75-4c60-4dd3-82f1-00840c3e0584" />
<img width="635" height="87" alt="Screenshot 2025-09-04 at 11 16 54 p m" src="https://github.com/user-attachments/assets/606b4111-e885-4853-a30d-2638682d0062" />
<img width="749" height="55" alt="Screenshot 2025-09-04 at 11 17 06 p m" src="https://github.com/user-attachments/assets/69a68b16-8bdb-4e21-b329-80e2b2f88ea0" />

- Likelihood Ratio Test: Compares full vs reduced models for each contrast
- Statistical thresholds: FDR < 0.05 and |log2FC| > 1 (2-fold change)
- Smear plots: Visualize log-fold changes versus average expression
- Comprehensive output: Full results tables plus separate up/down gene lists

# Pathway Enrichment Analysis Results - ShinyGo

## Overview
This part presents the results of pathway enrichment analysis performed on differentially expressed genes from the RNA-Seq study comparing Saccharomyces cerevisiae ste12 mutant under nitrogen deficiency versus nitrogen sufficiency conditions (ST → N- vs N+).

## Analysis Summary
The enrichment analysis reveals distinct biological pathways significantly affected by nitrogen availability in the ste12 mutant background, showing both up-regulated and down-regulated functional categories.

### 1. Enrichment Summary Table
The enrichment summary tables contain the main statistics of pathway analysis. They provide a concise view of which pathways are significantly enriched among up-regulated and down-regulated genes for each condition comparison.

Up-regulated ST → N- vs N+
<img width="364" height="260" alt="Screenshot 2025-09-04 at 11 28 04 p m" src="https://github.com/user-attachments/assets/b67a1911-b1a6-4456-9a09-8b87a551f07f" />

Down-regulated ST → N- vs N+
<img width="364" height="140" alt="Screenshot 2025-09-04 at 11 28 32 p m" src="https://github.com/user-attachments/assets/bee494d7-05f9-47ab-9543-0453eacb9d52" />

### 2. Pathway Enrichment
This output displays the detailed enrichment results, including p-values, adjusted p-values (FDR), and enrichment scores. It highlights biological pathways that are overrepresented in the dataset and helps identify condition-specific molecular mechanisms.

Up-regulated ST → N- vs N+
<img width="526" height="263" alt="Screenshot 2025-09-04 at 11 29 01 p m" src="https://github.com/user-attachments/assets/d85705eb-e816-431f-875b-a7e957b0328f" />

Down-regulated ST → N- vs N+
<img width="528" height="262" alt="Screenshot 2025-09-04 at 11 29 17 p m" src="https://github.com/user-attachments/assets/511d680e-e416-44b4-b80a-60d1c7192778" />

### 3. Tree Pathway Diagram
The tree diagram illustrates the hierarchical structure of enriched pathways. It organizes related pathways into clusters, making it easier to visualize how specific functional categories are interconnected.

Up-regulated ST → N- vs N+
<img width="463" height="340" alt="Screenshot 2025-09-04 at 11 29 40 p m" src="https://github.com/user-attachments/assets/9462d9e1-96bd-46c5-9a8b-54ede95d6280" />

Down-regulated ST → N- vs N+
<img width="482" height="341" alt="Screenshot 2025-09-04 at 11 29 55 p m" src="https://github.com/user-attachments/assets/66a1d4a9-a4a9-47a3-995e-2332b22b285c" />

### 4. Pathway Network
The pathway network visualizes the relationships between significantly enriched pathways. Nodes represent pathways, and edges represent shared genes, allowing us to understand cross-talk between biological processes.

Up-regulated ST → N- vs N+
<img width="490" height="422" alt="Screenshot 2025-09-04 at 11 30 18 p m" src="https://github.com/user-attachments/assets/b9d2a4df-286d-4f85-87b0-bef33a7d207e" />

Down-regulated ST → N- vs N+
<img width="487" height="418" alt="Screenshot 2025-09-04 at 11 30 30 p m" src="https://github.com/user-attachments/assets/fe67303b-32c9-4b34-b358-4a21771eb543" />

### 5. Significant KEGG Pathways
This section lists only the most significant KEGG pathways, filtered by statistical thresholds. It provides a focused view of the core pathways impacted by the experimental conditions.

Up-regulated ST → N- vs N+
<img width="437" height="482" alt="Screenshot 2025-09-04 at 11 30 53 p m" src="https://github.com/user-attachments/assets/e45b5ac1-ac4a-41ca-9879-8d6937688e75" />

Down-regulated ST → N- vs N+
<img width="539" height="483" alt="Screenshot 2025-09-04 at 11 31 08 p m" src="https://github.com/user-attachments/assets/4f9c0128-f011-4a7b-a18f-e4cb0eef6491" />

### Condition Comparisons
- ST → N- vs N+
Enrichment results comparing the stressed treatment (ST) under nitrogen deprivation (N-) versus nitrogen supplementation (N+).
- WT → N- vs N+
Enrichment results comparing wild type (WT) under N- versus N+.
- N+ → ST vs WT
Enrichment analysis comparing stressed treatment (ST) against wild type (WT) under nitrogen supplementation.
- N- → ST vs WT
Enrichment analysis comparing stressed treatment (ST) against wild type (WT) under nitrogen deprivation.
