### Shotgun proteomics analysis of HUVEC secretomes from 
###Activation of Notch signaling in endothelium cause upregulation of N-terminal acetylated histone 1

##Protein identification data
dat <- data.frame(read.csv("proteins_secr.csv"))


head(dat)
str(dat)
dat1 <- dat
rownames(dat1) <- dat[,3]
dat1 <- dat1[,c(18:29)]
head(dat1)

##sample info
library(readxl)
fact <- data.frame(read_excel("fact_secr.xlsx"))

rownames(fact) <- fact[,1]
fact$Type <- as.factor(fact$Type)
fact$Type


colnames(dat1) <- rownames(fact)


#Removing rows with a lot of missing values
dat3 <- dat1[which(rowMeans(!is.na(dat1)) >= 9/12), ]
mean(complete.cases(dat3))
colSums(is.na(dat3))

##Data normalization
library(RColorBrewer)

#tiff('Raw_dat.tiff', units="in", width=16, height=8, res=600, compression = 'lzw')
pal <- brewer.pal(n = 9, name = "Set1")
cols <- pal[fact$Type]
boxplot(dat3, outline = FALSE, col = cols, main = "Raw data")
legend("topright", levels(fact$Type), fill = pal, bty = "n", xpd = T)

library(impute)
dat_knn <- impute.knn(t(dat3), k = 5)
dat_knn <- t(dat_knn$data)
mean(complete.cases(dat_knn))

#tiff('Raw_dat.tiff', units="in", width=16, height=8, res=600, compression = 'lzw')
boxplot(dat_knn, outline = FALSE, col = cols, main = "Data after missed values impuration")
legend("topright", levels(fact$Type), fill = pal, bty = "n", xpd = T)
#dev.off()

colSums(dat_knn)
head(dat_knn)
dat_log <- log2(dat_knn+1)

library(limma)
data_norm <- normalizeQuantiles(dat_log)
boxplot(data_norm, outline = FALSE, col = cols, main = "Normalized data")
legend("topright", levels(fact$Type), fill = pal, bty = "n", xpd = T)



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

maplot(data_norm[, rownames(subset(fact,Type=="NICD"))], data_norm[, rownames(subset(fact,Type=="Control"))], main = "Normalized data")



##PCA ordination
library(mixOmics)
dat_pca <- pca(t(data_norm), ncomp = 8, center = TRUE)

#tiff('PCA_group.tiff', units="in", width=10, height=8, res=600, compression = 'lzw')
plotIndiv(dat_pca, comp = c(1, 2), ind.names = T, 
          group = fact$Type, legend = TRUE, ellipse = T,
          title = 'PCA')
#dev.off()

##Limma analysis of differential expression
X <- model.matrix(~ fact$Type)
X

fit <- lmFit(data_norm, design = X, method = "robust", maxit = 10000)

efit <- eBayes(fit)

topTable(efit, coef = 2)
full_list <- topTable(efit, coef = 2, number = length(data_norm[,2]))
#write.csv(full_list,'Contr_vs_NICD.csv')
head(full_list)


library(EnhancedVolcano)

#tiff('Dif_expr.tiff', units="in", width=15, height=11, res=300, compression = 'lzw')
EnhancedVolcano(full_list,
                lab = rownames(full_list),
                x = 'logFC',
                y = 'P.Value',
                pCutoff = 0.05,
                FCcutoff = 1.5,
                #xlim = c(-3, 2),
                #ylim = c(0, 10),
                title ="Control vs NICD",
                labSize = 4.0,
                boxedLabels = F,
                colAlpha = 1)
#dev.off()


p_aboveS <- full_list$adj.P.Val <= 0.05
sum(p_aboveS)

head(full_list)

##Manual check of top-proteins
dat_knn1 <- as.matrix(dat_knn)
dat3_1 <- as.matrix(dat3)
boxplot(dat3_1[c("P26022"),] ~ Type, data = fact,
        varwidth = TRUE, log = "y", las = 1)
boxplot(dat3_1[c("P05121"),] ~ Type, data = fact,
        varwidth = TRUE, log = "y", las = 1)
boxplot(dat3_1[c("P07996"),] ~ Type, data = fact,
        varwidth = TRUE, log = "y", las = 1)