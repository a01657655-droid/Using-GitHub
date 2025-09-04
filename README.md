# Using GitHub

## RNAseq
### Joana Paulina Pineda Corella 
### A01657655

```{r}
#checa si el paquete BiocManager ya está instalado y lo carga desde CRAN
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("edgeR")
BiocManager::install("limma")

library(limma)
library(edgeR) #análisis de expresión diferencial en RNAseq
```

```{r}
#get working directory - te devuelve la carpeta actual en la que R está trabajando 
getwd()

#crear una carpeta dentro de tu directorio actual
dir.create("RNAseq_project", showWarnings = FALSE)

#set working directory - te mueves a la carpeta que indiques
setwd("RNAseq_project")

dir.create("results_diffexp", showWarnings = FALSE)

#muestra todos los archivos y carpetas que hay en el working directory actual 
list.files()

#guardar la ruta de la carpeta results en variable outpath
outpath = "results_diffexp"
dir.create(outpath, showWarnings=FALSE,recursive=TRUE)
```

```{r}
#leer archivos de texto en R
#header = la primera fila del archivo contiene los nombres de las columnas 
#row = la primera columna del archivo se usará como nombre de las filas 
#sep = indica que el archivo está separado por tabulaciones 
counts <- read.table("https://raw.githubusercontent.com/jmvillalobos/RNAseq_curso2025/main/Saccharomyces.txt",
                     header = TRUE, row.names = 1, sep = "\t")

#muestra las primeras 6 filas de la tabla 
head(counts)

#devuelve un vector con número de filas y columnas 
dim(counts)
```

```{r}
#cpm = counts per million 
#rowSums = suma los TRUE por fila 
#counts = se usa para seleccionar solo las filas (genes) que cumplen la condición 
counts = counts[rowSums(cpm(counts) >= 2) >=10,]

head(counts)
colnames(counts) #nombres de las columnas (muestras) para verificar que no se modificaron 
dim(counts)
```

```{r}
#gráficar
#log2 para aplanar diferencias grandes en conteos 
plot(log2(counts[,c("wt_sc_1", "st_sc_1")]), col="violet") 
```

```{r}
#crear los grupos
#colnames(counts) = devuelve los nombres de las columnas de la tabla
#grp = #vector que indica a qué grupo pertenece cada muestra 
grp = sub("..$", "", colnames(counts)) 

grp #imprime el vector para comprobar que quedo bien 
dge = DGEList(counts=counts, group=grp) #crea un objeto especial para análisis de expresión diferencial 
dge #revisar objeto 
```

```{r}
#MDS = multidimensional scaling - visualiza similitudes o diferencias entre muestras en un gráfico
plotMDS(dge, col="violet")
```

```{r}
#calcNormFactors() calcula factores de normalización para tus muestras
dgeNorm = calcNormFactors(dge)

#tabla de información de las muestras que ahora incluye los factores de normalización calculados 
dgeNorm$samples
```

```{r}
#el comando estima la dispersión común de todas las librerias - una medida promedio de qué tan variables son las muestras entre sí
dgeNorm = estimateCommonDisp(dgeNorm)

#te devuelve un número que representa la variabilidad biológica compartida entre todos los genes 
dgeNorm$common.dispersion
```

```{r}
#exactTest() = para saber si un gen tiene diferencias significativas en expresión entre dos grupos de muestras 
diff_exp = exactTest(dgeNorm, dispersion = dgeNorm$common.dispersion, pair = c("wt_sc", "wt_sl" ))
diff_exp2 = exactTest(dgeNorm, dispersion = dgeNorm$common.dispersion, pair = c("st_sc", "st_sl" ))

diff_exp
```

```{r}
#abre la documentación de la función y explica qué onda 
?exactTest

#da el número de genes analizados (filas) y las columnas con resultados 
dim(diff_exp)
```

```{r}
#toma el objeto DGEExact, ordena los genes según significancia estadística y muestra los genes con mayor cambio y los valores asociados 
topTags(diff_exp)
```

```{r}
#topTags(diff_exp) = tabla ordenada de genes según la significancia de la diferencia de expresión 
dim(topTags(diff_exp))
```

```{r}
#topTags() = top 10 genes más significativos 
#n=Inf = muestra todos los genes en la tabla 
#$table = extrae solo la tabla de resultados dentro del objeto que devuelve topTags()
deTab = topTags(diff_exp, n=Inf)$table

#selecciona filas 15 y 30 de deTab
deTab[c(15,30),]
```

```{r}
#da un vector con los nombres de los genes más sobreexpresados según el criterio
row.names(deTab)[deTab$logFC > 5]  #recuerden que el FC esta dado en log2.
```

```{r}
#da toda la información de ese gen
deTab["YNL284C-A",]
```

```{r}
#filtrar automáticamente por cualquier valor que predeterminemos
#valores estándar de filtrado pero esto se puede modificar
deGenes = rownames(deTab)[deTab$FDR < 0.05 & abs(deTab$logFC) > 1]
down=row.names(deTab)[deTab$logFC< -1]
over=row.names(deTab)[deTab$logFC> 1]

print(paste("total de diferenciales:", length(deGenes)))
print(paste("número de genes inducidos:", length(over)))
print(paste("número de genes reprimidos:", length(down)))
```

```{r}
plotSmear(dge, de.tags=deGenes, ylab = "WT-sc vs WT-sl")
#gráfico que permite ver cómo se distribuyen los genes según su expresión y cambio entre condiciones 
```

```{r}
#exporta un objeto de R a un archivo de texto 
write.table(deTab, file=paste(outpath, "diff_gene_wt.txt", sep=""), row.names=TRUE, col.names=NA, quote=FALSE, sep="\t")

#aquí exportamos solo los genes que se sobreexpresan o se reprimen
write.table(down, file=paste(outpath, "down.txt", sep=""), row.names=TRUE, col.names=NA, quote=FALSE, sep="\t")
write.table(over, file=paste(outpath, "up.txt", sep=""), row.names=TRUE, col.names=NA, quote=FALSE, sep="\t")
```

```{r}
install.packages(
  c("RColorBrewer",  "ggplot2"),
  repos = "https://cloud.r-project.org"
)

library("gplots")
library("RColorBrewer")
```

```{r}
#normalizamos nuestros datos de expresión por cuentas por millón
normalizados= cpm(counts)

#extraemos la expresión de los genes diferenciales
normalizados_diferenciales= normalizados[deGenes,]

#veamos cómo se ve esta tablita
head(normalizados_diferenciales)
```

```{r}
library(RColorBrewer)
my_colors <- colorRampPalette(c("purple", "magenta", "violet"))(100)

heatmap(normalizados_diferenciales, col=my_colors, scale="row", margins=c(5,10))
```

```{r}
#configura la disposición de los gráficos en la ventana de R
#mfrow = 1 fila 2 columnas - se muestran dos gráficos lado a lado
par(mfrow=c(1,2)) 

#diagrama de caja para cada columna de los datos
boxplot(log(counts),col=rainbow(6), main="antes de la normalización")
boxplot(log(normalizados), col=rainbow(6), main="después de la normalización")
```

```{r}
barplot(apply(normalizados_diferenciales,2,sum),las=2, cex.names = 1, col = (1:6))
```

```{r}
pca <- princomp(normalizados_diferenciales[,c(1:6)])
plot(pca$loadings, col=as.factor(colnames(normalizados_diferenciales[,c(1:6)])),  pch=19, cex=2, main="con nitrógeno")
text(pca$loadings, as.vector(colnames(normalizados_diferenciales[,c(1:6)])), pos=3, cex=0.8)
```

```{r}
pca <- princomp(normalizados_diferenciales[,c(7:12)])
plot(pca$loadings, col=as.factor(colnames(normalizados_diferenciales[,c(7:12)])),  pch=19, cex=2, main="si nitrógeno")
text(pca$loadings, as.vector(colnames(normalizados_diferenciales[,c(7:12)])), pos=3, cex=0.8)
```

```{r}
with(deTab, plot(logFC, -log10(FDR), pch=20, cex=0.8, col="purple", main="WT+N vs WT-N", xlim=c(-8, 8), ylim=c(0,300)))
text(deTab[1:20,]$logFC,-log(deTab[1:20,]$FDR,10),labels=rownames(deTab[1:20,]),cex=0.7,pos=1)
with(subset(deTab, FDR<.01 & abs(logFC)>2), points(logFC, -log10(FDR), pch=20, cex=0.5, col="pink"))
abline(v=2,lty=2, col="blue")
abline(v=-2,lty=2, col="blue")
legend("bottomright","Up_regulated",cex=1)
legend("bottomleft","Down_regulated",cex=1)
```

```{r}
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

```{r}
#los 30 genes principales sobreexpresados en la comparación WT 
WTover= head(rownames(deTab), 30)

#los 30 genes principales sobreexpresados en la comparación ST 
ste12over= head(rownames(deTab2), 30)

#da un vector con los nombres de esos genes exclusivos de WT
setdiff(WTover, ste12over)
```
