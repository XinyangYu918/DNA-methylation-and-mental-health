#### Perform EWAS analysis for twins study ####
library(minfi)
library(minfiData)
library(gee)

# Provide cohort name
cohort = "YourCohortName"
load("./Quan-norm.rda")
load("./RGset.rda")

# Probe quality check
# Add SNP info to the data
objectWithSNPinfo <- addSnpInfo(object)

# Drop probes that contain either an SNP at the CpG interrogation or at the single nucleotide extension
objectSNPQCed <- dropLociWithSnps(objectWithSNPinfo, snps=c("SBE", "CpG", "Probe"), maf=0.05)
rm(objectWithSNPinfo)

detP <- detectionP(RGset)
Match1 <- match(colnames(objectSNPQCed),colnames(detP))
Match2 <- match(rownames(objectSNPQCed),rownames(detP))
detPSNPQCed <- detP[Match2[!is.na(Match2)],Match1[!is.na(Match1)]]
rm(Match1,Match2,detP)

failed <- detPSNPQCed >0.01
rm(detPSNPQCed)
beta <- getBeta(objectSNPQCed)

# Drop probes that failed quality control via the detection p-value in greater than 20% of samples
failedCG02 <- rowMeans(failed)>0.2

# Get the list non-variable CpG sites i.e. those where beta values for all samples are <= 20% or >= 80%
ProbeInvar <- (rowSums(beta<=0.2)==ncol(beta))|(rowSums(beta>=0.8)==ncol(beta))
ListInvarProbe <- rownames(beta)[which(ProbeInvar)]

# Remove sex chromosome probes
keepIndex=!seqnames(objectSNPQCed) %in% c("chrX", "chrY") 
rm(objectSNPQCed)

keepIndex <- keepIndex&(!failedCG02)
rm(failedCG02)

beta[failed] <- NA
rm(failed)

betaQC <- beta[which(keepIndex),]
rm(beta,keepIndex)

# Reformat beta value
Methy <- as.data.frame(t(betaQC))
Methy$Subject <- colnames(betaQC)
rm(betaQC);gc()
Methy <- Methy[,c(ncol(Methy),1:(ncol(Methy)-1))] 
MethyName <- colnames(Methy)[-c(1)]

# Load and format covariates
load("./fast_svd.rda")
PC_Beta <- as.data.frame(ss$v[,1:4])
PC_Beta$Subject <- Methy$Subject
load("./cellcount.rda") 
tmp <- prcomp(cellcount) 
pc_cell <- as.data.frame(tmp$x[,1:2]) 
pc_cell$Subject <- rownames(pc_cell) 

# Read and format subcortical volumes
Raw_Struc <- read.csv("./LandRvolumes.csv", header=T) 
Struc <- matrix(data=0, ncol=8, nrow=nrow(Raw_Struc)) 
colnames(Struc) <- c("Subject","Mthal", "Mcaud", "Mput", "Mpal", "Mhippo", "Mamyg", "Maccumb") 
Struc <- as.data.frame(Struc) 

Struc$Subject = Raw_Struc$SubjID 
Struc$Mthal <- rowMeans(Raw_Struc[,c("Lthal","Rthal")]);
Struc$Mcaud <- rowMeans(Raw_Struc[,c("Lcaud","Rcaud")]);
Struc$Mput <- rowMeans(Raw_Struc[,c("Lput","Rput")]); 
Struc$Mpal <- rowMeans(Raw_Struc[,c("Lpal","Rpal")]); 
Struc$Mhippo <- rowMeans(Raw_Struc[,c("Lhippo","Rhippo")]); 
Struc$Mamyg <- rowMeans(Raw_Struc[,c("Lamyg","Ramyg")]); 
Struc$Maccumb <- rowMeans(Raw_Struc[,c("Laccumb","Raccumb")]);
StrucName <- colnames(Struc)[-1]
ICV <- Raw_Struc[,c(1,ncol(Raw_Struc))]
colnames(ICV)[1] <- "Subject" 

# Combine covariates
Raw_Cov <- read.csv("./SubCortCovs.csv", header=T) 
Raw_Cov$Age_Square <- (Raw_Cov$Age)^2
colnames(Raw_Cov)[2] <- "Subject"
Cov <- merge(Raw_Cov, PC_Beta, by="Subject", all=F) 
Cov <- merge(Cov, pc_cell, by="Subject", all=F) 
Cov <- merge(Cov, ICV, by="Subject", all=F) 
CovName <- colnames(Cov)[-c(1,2)] 

# Generate files with complete data across all data
Data <- merge(Cov, Struc, by="Subject", all=F) 
Match <- match(Methy$Subject, Data$Subject)

# Generate new covariates, structure and methylation data, which have had their individuals matched
Cov <- Data[Match[!is.na(Match)], c("FID", "Subject", CovName)]    
Struc <- Data[Match[!is.na(Match)], c("Subject", StrucName)]  
Methy <- Methy[!is.na(Match),]

# For GEE, data MUST be ordered by family ID. Script below assumes that the covariates file contains a column called "FID" which contains the family ID/twin pair ID, which is a numeric variable.
rownames(Cov) <- Cov$Subject
rownames(Struc) <- Struc$Subject
rownames(Methy) <- Methy$Subject

Cov   <- Cov[order(Cov$FID),]       
Struc <- Struc[rownames(Cov),]      
Methy <- Methy[rownames(Cov),]      

# remove Subject ids
Methy <- Methy[,-c(1)]
Cov <- Cov[,-c(2)]
Struc  <- Struc[,-c(1)]

# Perform EWAS analysis: linear regression analysis between DNA methylation and subcortical volumes
Num_Methy <- ncol(Methy)
Num_Cov <- ncol(Cov)-1
Num_Struc <- ncol(Struc)

# Prepare the output files
Origin_Beta <- matrix( data= NA, nrow=Num_Methy, ncol= Num_Struc, byrow=F, dimnames=NULL)
colnames(Origin_Beta) <- colnames(Struc)
rownames(Origin_Beta) <- colnames(Methy)
Origin_SE <- matrix( data= NA, nrow=Num_Methy, ncol= Num_Struc, byrow=F, dimnames=NULL)   
colnames(Origin_SE) <- colnames(Struc)
rownames(Origin_SE) <- colnames(Methy)
Origin_P <- matrix( data= NA, nrow=Num_Methy, ncol= Num_Struc, byrow=F, dimnames=NULL)
colnames(Origin_P) <- colnames(Struc)
rownames(Origin_P) <- colnames(Methy)

# GEE analysis
Covar <- matrix(unlist(Cov[,-c(1)]), ncol=ncol(Cov)-1, byrow=F)
for (i in 1:Num_Methy) {
	for (j in 1:Num_Struc) {
		Out <- summary(gee(Methy[,i]~Covar+Struc[,j],id=Cov$FID, family=gaussian, corstr="exchangeable", maxiter=100, na.action=na.omit))
		Origin_Beta[i,j] <- Out$coefficients[nrow(Out$coefficients),1]
		Origin_SE[i,j] <- Out$coefficients[nrow(Out$coefficients),4]
		Origin_P[i,j] <- 2*pnorm(-abs(Out$coefficients[nrow(Out$coefficients),5]))
	}
}

# Save output file
save(Origin_Beta, Origin_SE, Origin_P, ListInvarProbe, file=paste("./Output_of_",cohort, "_Methylation_and_All_Subcortical_Structure.RData", sep=""), compress = T)

