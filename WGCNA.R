#### Perform WGCNA model ####
library(WGCNA);
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

# Load data
datExpr0 <- readRDS("Resid_14_for_WGCNA.rds") # Residuals of methylation data, regressed out for recruitment sites

# Identify bad samples
sampleTree = hclust(dist(datExpr0), method = "average");
#save(sampleTree, file = "./sampleTree for sample at age 14.rda")

# Plot the sample tree
sizeGrWindow(12,9)
pdf(file = "sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0, 4, 2, 0))
plot(sampleTree, 
     main = "Sample clustering to detect outliers", 
     sub="", 
     xlab="", 
     cex.lab = 1.5,
     cex.axis = 1.5, 
     cex.main = 2)
dev.off()

# Select the optimal soft threshold
powers <- c(c(1:10), seq(from = 12, to=20, by=2))
sft <- pickSoftThreshold(datExpr0[,-c(1)], powerVector = powers, verbose = 5, blockSize = 50000)
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
save(bwnet, file = "./bwnet for WGCNA at age 14.RData")
