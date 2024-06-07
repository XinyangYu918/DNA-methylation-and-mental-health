#### Perform WGCNA model ####
library(WGCNA);
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

# Load data
datExpr0 <- readRDS("Resid_14_for_WGCNA.rds") # Residuals of methylation data, regressed out for recruitment sites and waves (hybridisation dates)

# Check bad samples
gsg <- goodSamplesGenes(datExpr0[,-c(1)], verbose = 3);
gsg$allOK

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

# Plot the sample tree
sampleTree <- hclust(dist(datExpr0), method = "average");
#save(sampleTree, file = "./sampleTree for sample at age 14.rda")
sizeGrWindow(12, 9)
#pdf(file = "sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0, 4, 2, 0))
plot(sampleTree, 
     main = "Sample clustering to detect outliers", 
     sub = "", 
     xlab = "", 
     cex.lab = 1.5,
     cex.axis = 1.5, 
     cex.main = 2)

# Plot a line to show the cut
abline(h = 15, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 15, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# Select the optimal soft threshold
# Choose a set of soft-thresholding powers
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
# Call the network topology analysis function
sft <- pickSoftThreshold(datExpr0[,-c(1)], powerVector = powers, verbose = 5, blockSize = 50000)

# Scale-free topology fit index as a function of the soft-thresholding power
sizeGrWindow(9, 5)
par(mfrow = c(1,2));

plot(sft$fitIndices[,1], 
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], 
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels = powers,
     cex = 0.9,
     col = "red")
abline(h = 0.80, col = "red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], 
     sft$fitIndices[,5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity", 
     type = "n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], 
     sft$fitIndices[,5], 
     labels = powers, 
     cex = 0.9,
     col = "red")

save(sft, file = "./sft for WGCNA at age 14.RData")

# Here, the best power is 5.
# Identify modules using WGCNA
bwnet <- blockwiseModules(datExpr0[,-c(1)], 
                          maxBlockSize = 50000,
                          power = 5, 
                          TOMType = "signed", 
                          minModuleSize = 30,
                          reassignThreshold = 0, 
                          mergeCutHeight = 0.25,
                          numericLabels = TRUE,
                          saveTOMs = TRUE,
                          saveTOMFileBase = "TOM-blockwise-age14-power5-min30-height25",
                          verbose = 3)

# Plot the dendrogram and the module colors underneath
sizeGrWindow(12, 9)
mergedColors <- labels2colors(bwnet$colors)
plotDendroAndColors(bwnet$dendrograms[[1]], 
                    mergedColors[bwnet$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, 
                    hang = 0.03,
                    addGuide = TRUE, 
                    guideHang = 0.05)

save(bwnet, file = "./bwnet for WGCNA at age 14.RData")
