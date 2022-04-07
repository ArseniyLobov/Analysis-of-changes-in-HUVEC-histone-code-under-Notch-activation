### Shotgun proteomics analysis of histone enriched fraction from 
###Activation of Notch signaling in endothelium cause upregulation of N-terminal acetylated histone 1

##Protein identification data
dat_s <- data.frame(read.csv("proteins_hist_pr.csv"))
dat_s[,3] <- gsub("\\|.*","",dat_s[,3])
rownames(dat_s) <- dat_s[,3]
dat_s <- dat_s[,18:29]

##Sample info
library(readxl)
fact_s <- data.frame(read_excel("fact_pr.xlsx"))

fact_s$Condition <- as.factor(fact_s$Condition)
fact_s$Condition

rownames(fact_s) <- fact_s[,1]
colnames(dat_s) <- rownames(fact_s)

head(dat_s)
head(fact_s)

##removing sample with high amount of NA
dat_s <- dat_s[,-c(5,11)]
fact_s <- fact_s[-c(5,11),]

##removing proteins with too much NA
dat_s1 <- dat_s[which(rowMeans(!is.na(dat_s)) >= 8/10), ]
mean(complete.cases(dat_s1))
colSums(is.na(dat_s1))

##removing protein contaminations
dat_s1 <- dat_s1[!rownames(dat_s1) %in% c("P08729|K2C7_HUMAN","P05783|K1C18_HUMAN" ,"P35527|K1C9_HUMAN","P04264|K2C1_HUMAN", "P13645|K1C10_HUMAN", "P08779|K1C16_HUMAN", "P02533|K1C14_HUMAN"),]

## Imputation of missed values
library(impute)
tdats <- t(dat_s1)
dat_knns <- impute.knn(tdats, k = 5)
dat_knns <- t(dat_knns$data)
mean(complete.cases(dat_knns))


##Data normalization
library(RColorBrewer)
pal <- brewer.pal(n = 9, name = "Set1")
cols <- pal[fact_s$Condition]
boxplot(dat_knns, outline = FALSE, col = cols, main = "Raw data")
legend("topright", levels(fact_s$Condition), fill = pal, bty = "n", xpd = T)

dat_logs <- log2(dat_knns+1)
head(dat_logs)
mean(complete.cases(dat_logs))
boxplot(dat_logs, outline = FALSE, col = cols, main = "Log-transformed data")
legend("topright", levels(fact_s$Condition), fill = pal, bty = "n", xpd = T)

library(limma)
dat_norms <- normalizeQuantiles(dat_logs)
head(dat_norms)
boxplot(dat_norms, col = cols, main = "Normalized data")
legend("topright", levels(fact_s$Condition), fill = pal, bty = "n", xpd = T)
mean(complete.cases(dat_norms))
colSums(is.na(dat_norms))

#MA-plot
maplot <- function(X1, X2, pch = 21, main = "MA-plot", xlab = "Average log-expression", ylab = "Expression log-ratio", lpars = list(col = "blue", lwd = 2), ...){
  X <- (rowMeans(X2) + rowMeans(X1)) / 2
  Y <- rowMeans(X2) - rowMeans(X1)
  scatter.smooth(x = X, y = Y,
                 main = main, pch = pch,
                 xlab = xlab, ylab = ylab,
                 lpars = lpars, ...)
  abline(h = c(-1, 0, 1), lty = c(2, 1, 2))
}

maplot(dat_norms[, rownames(subset(fact_s,Condition=="NICD"))], dat_norms[, rownames(subset(fact_s,Condition=="Control"))], main = "Normalized data")


##Limma analysis of differential expression
X_s <- model.matrix(~ fact_s$Condition)
X_s

fits <- lmFit(dat_norms, design = X_s, method = "robust", maxit = 10000)

efits <- eBayes(fits)

topTable(efits, coef = 2)
numGenes_s <- length(dat_norms[,2])
full_list_efits <- topTable(efits, number = numGenes_s)
#write.csv(full_list_efits,'shotgun_dif.csv')
head(full_list_efits)

library(EnhancedVolcano)
#tiff('Vulcano_shotgun.tiff', units="in", width=11, height=8, res=600, compression = 'lzw')
EnhancedVolcano(full_list_efits,
                lab = rownames(full_list_efits),
                x = 'logFC',
                y = 'adj.P.Val', # что такое adj. p.val? а почему мы просто p. val не используем??
                pCutoff = 0.05,
               # xlim = c(-5, 5),
              #  ylim = c(0, 9),
                FCcutoff = 1,
                title ="Vulcano",
                labSize = 4.0,
                boxedLabels = F,
                colAlpha = 1)
#dev.off()


##Manual check of top-proteins
head(full_list_efits)
dat_knns <- as.matrix(dat_knns)
boxplot(dat_knns[c("O94875"),] ~ Condition, data = fact_s,
        varwidth = TRUE, log = "y", las = 1)
boxplot(dat_knns[c("P51397"),] ~ Condition, data = fact_s,
        varwidth = TRUE, log = "y", las = 1)
boxplot(dat_knns[c("P02533"),] ~ Condition, data = fact_s,
        varwidth = TRUE, log = "y", las = 1)
boxplot(dat_knns[c("Q8N3V7"),] ~ Condition, data = fact_s,
        varwidth = TRUE, log = "y", las = 1)
boxplot(dat_knns[c("P28072"),] ~ Condition, data = fact_s,
        varwidth = TRUE, log = "y", las = 1)
boxplot(dat_knns[c("Q13823"),] ~ Condition, data = fact_s,
        varwidth = TRUE, log = "y", las = 1)

##PCA ordination
library(mixOmics)
dat_pca <- pca(t(dat_norms), ncomp = 8, center = TRUE)

#tiff('PCA_group.tiff', units="in", width=10, height=6, res=600, compression = 'lzw')
plotIndiv(dat_pca, comp = c(1, 2), ind.names = F, 
          group = fact_s$Condition, legend = TRUE, ellipse = T,
          title = 'PCA')
#dev.off()




### N-acetylated peptide analysis
dat <- data.frame(read.csv("protein-peptides.csv"))

head(dat)
##Extraction of N-acetylated peptides
dat$PTM <- as.factor(dat$PTM)
dat1 <- dat[dat$PTM=="Acetylation (Protein N-term)",]
new_names <- with(dat1, paste0(Protein.Accession, Peptide))
dat2 <- dat1[,c(14:25)]
rownames(dat2) <- new_names
head(dat2)
str(dat2)

##sample info
fact <- data.frame(read_excel("fact_pep.xlsx"))

rownames(fact) <- fact[,1]
fact <- fact[,-1]
fact$NICD <- as.factor(fact$NICD)
fact$NICD

fact$Donor <- as.factor(fact$Donor)
fact$Donor

colnames(dat2) <- rownames(fact)


##removing peptides with too much NA
dat3 <- dat2[which(rowMeans(!is.na(dat2)) >= 9/12), ]
mean(complete.cases(dat3))
colSums(is.na(dat2))

dat3[is.na(dat3)] <- 0
mean(complete.cases(dat3))

##removing samples with too much missed values
dat3 <- dat3[,-c(5,11)]
fact <- fact[-c(5,11),]

## Data normalization
pal <- brewer.pal(n = 9, name = "Set1")
cols <- pal[fact$NICD]
boxplot(dat3, outline = FALSE, col = cols, main = "Raw data")
legend("topright", levels(fact$NICD), fill = pal, bty = "n", xpd = T)
colSums(dat3)

dat_log <- log2(dat3+1)
head(dat_log)
mean(complete.cases(dat_log))
boxplot(dat_log, outline = FALSE, col = cols, main = "Log-transformed data")
legend("topright", levels(fact$NICD), fill = pal, bty = "n", xpd = T)

dat_norm <- normalizeQuantiles(dat_log)
head(dat_norm)
boxplot(dat_norm, col = cols, main = "Normalized data")
legend("topright", levels(fact$NICD), fill = pal, bty = "n", xpd = T)
mean(complete.cases(dat_norm))
colSums(is.na(dat_norm))


##Limma analysis of differential expression
X <- model.matrix(~ fact$NICD)
X

fit <- lmFit(dat_norm, design = X, method = "robust", maxit = 10000)

efit <- eBayes(fit)

topTable(efit, coef = 2)
numGenes <- length(dat3$N5_3)
full_list_efit <- topTable(efit, number = numGenes)
#write.csv(full_list_efit,'Dif_hist_pept.csv')
head(full_list_efit)
#tiff('Vulcano.tiff', units="in", width=11, height=8, res=600, compression = 'lzw')
EnhancedVolcano(full_list_efit,
                lab = rownames(full_list_efit),
                x = 'logFC',
                y = 'adj.P.Val', 
                pCutoff = 0.05,
               #xlim = c(-8, 10),
               #ylim = c(0, 5),
                FCcutoff = 1, 
                title ="Vulcano",
                labSize = 4.0,
                boxedLabels = F,
                colAlpha = 1)
#dev.off()