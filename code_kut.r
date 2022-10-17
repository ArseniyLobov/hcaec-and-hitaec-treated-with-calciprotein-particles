##Готовим данные_HCA
#protein_expression_data
setwd("C:/R/kut_2")

dat <- data.frame(read.csv("proteins_HCA.csv"))


dat$Accession <- sub("\\|.*", "", dat$Accession)                   # Extract first three characters

dat1 <- dat[,c(3,18:29)]
head(dat1)

str(dat1)
rownames(dat1) <- dat1[,1]
dat1 <- dat1[,-1]
head(dat1)



#делаем таблицу для количественного анализа
#по какому принципу я обрезаю датафрейм ниже?
dat2 <- dat1[which(rowMeans(!is.na(dat1)) >= 0.90), ]
mean(complete.cases(dat2))
colSums(is.na(dat2))

#Импутируем пропущенные значения по К-ближайшим соседям. Гуглим про это и читаем методичку :)
library(impute)
tdat <- t(dat2)
dat_knn1 <- impute.knn(tdat, k = 5)
dat_knn <- t(dat_knn1$data)
mean(complete.cases(dat_knn))

#Грузим нашу факторную матрицу (легенда к биологическим группам в образцах)
fact <- data.frame(c("sample", colnames(dat2)), c("type", "K", "K", "K", "M", "M", "M", "C", "C", "C", "I", "I", "I"), c("typeB", "KM", "KM", "KM", "KM", "KM", "KM", "CI", "CI", "CI", "CI", "CI", "CI"))
colnames(fact) <- fact[1,]
fact <- fact[-1,]
str(fact)
rownames(fact) <- fact[,1]

fact$type <- as.factor(fact$type)

fact$type

fact$typeB <- as.factor(fact$typeB)

fact$typeB

fact$typeB <- relevel(fact$typeB, ref = "KM")



#Контроль должен быть базовым уровнем !!! если нет, используем функцию relevel

#смотрим как выглядят данные и делаем их трансформацию - чтобы средние и распределения совпадали.
#сырые
library(RColorBrewer)
pal <- brewer.pal(n = 9, name = "Set1")
cols <- pal[fact$type]
boxplot(dat_knn, outline = FALSE, col = cols, main = "Raw data")
legend("topright", levels(fact$type), fill = pal, bty = "n", xpd = T)
colSums(dat_knn)
#log transformation (погуглите зачем это делать). В чем разница с тем, что было? почему так удобно? зачем прибавлять 1 к dat_knn?
dat_log <- log2(dat_knn+1)
head(dat_log)
mean(complete.cases(dat_log))
boxplot(dat_log, outline = FALSE, col = cols, main = "Log-transformed data")
legend("topright", levels(fact$type), fill = pal, bty = "n", xpd = T)

#Quantile normalization
library(limma)
dat_norm <- normalizeQuantiles(dat_log)
head(dat_norm)
boxplot(dat_norm, col = cols, main = "Normalized data")
legend("topright", levels(fact$type), fill = pal, bty = "n", xpd = T)
mean(complete.cases(dat_norm))
colSums(is.na(dat_norm))


## Ординации
t_dat1 <- t(dat_norm)
t_dat1 <- t_dat1[-c(4:6),]
fact <- fact[-c(4:6),]
fact$type <- as.factor(as.character(fact$type))

library(mixOmics)
#PLS-DA 
ordination.optimum.splsda <- splsda(t_dat1, fact$type, ncomp = 3, keepX = c(15,15,15))
selectVar(ordination.optimum.splsda, comp=1)
selectVar(ordination.optimum.splsda, comp=2)
selectVar(ordination.optimum.splsda, comp=3)

tiff('PLSDA_HCA_three.tiff', units="in", width=16, height=8, res=600, compression = 'lzw')
layout(matrix(c(1, 2, 3, 3, 3, 3), 2, 3))
plotLoadings(ordination.optimum.splsda, comp = 1, size.name = 1, size.title = 1.2, title = "Loadings\n on 1st component", contrib = "max", legend = FALSE, col.ties="black", ndisplay = 15)
plotLoadings(ordination.optimum.splsda, comp = 2, size.name = 1, size.title = 1.2, title = "Loadings\n on 2nd component", contrib = "max",ndisplay = 15,  legend = FALSE, col.ties="black")
plotIndiv(ordination.optimum.splsda, ind.names = F, ellipse = T, style = "graphics", abline = TRUE, cex = 1.5, size.axis = 1.2, size.xlabel = 1.5, size.ylabel = 1.5, title = "PLS-DA ordination", size.title = 1.5, legend=TRUE)
dev.off()
layout(1,1)



## Ординации. Пока только PCA


dat_pca <- pca(t_dat1, ncomp = 8, center = TRUE)
dat_pca
plot(dat_pca)
dat_pca <- pca(t_dat1, ncomp = 3, center = TRUE)

plotIndiv(dat_pca, comp = c(1, 2), ind.names = F, 
          group = fact$type, legend = TRUE, ellipse = TRUE,
          title = 'PCA')

## Сохраняем. Что нового появилось на рисунке?
tiff('PCA_HCA_three.tiff', units="in", width=16, height=8, res=600, compression = 'lzw')
layout(matrix(c(1, 2, 3, 3, 3, 3), 2, 3))
plotLoadings(dat_pca, comp = 1, size.name = 1, size.title = 1.2, title = "Loadings\n on 1st component", contrib = "max", legend = FALSE, col.ties="black", ndisplay = 15)
plotLoadings(dat_pca, comp = 2, size.name = 1, size.title = 1.2, title = "Loadings\n on 2nd component", contrib = "max",ndisplay = 15,  legend = FALSE, col.ties="black")
plotIndiv(dat_pca, comp = c(1, 2), ind.names = F, 
          group = fact$type, legend = TRUE, ellipse = TRUE, style = "graphics", 
          title = 'PCA')
dev.off()
layout(1,1)

#К с С
datKC <- data.frame(dat_norm[,c(1:3, 7:9)])
facKC <- fact[c(1:3, 7:9),]
str(datKC)
facKC$type <- as.factor(as.character(facKC$type))
facKC$type <- relevel(facKC$type, ref = "K")

Xkc <- model.matrix(~ facKC$type)
Xkc

fitkc <- lmFit(datKC, design = Xkc, method = "robust", maxit = 10000)

# Empirical Bayes statistics
efit_kc <- eBayes(fitkc)

# Dif_expr_table
topTable(efit_kc, coef = 2)
full_list_kc <- topTable(efit_kc, number = length(dat_norm))
write.csv(full_list_kc,'Dif_expr_K_vs_C.csv')
head(full_list_kc)
#Vulcano
library(EnhancedVolcano)

tiff('Vulcano_K_vs_C.tiff', units="in", width=11, height=8, res=600, compression = 'lzw')
EnhancedVolcano(full_list_kc,
                lab = rownames(full_list_kc),
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                xlim = c(-10, 11),
                ylim = c(0, 7),
                FCcutoff = 2,
                title ="K_vs_C",
                labSize = 4.0,
                boxedLabels = F,
                colAlpha = 1)
dev.off()


#К с I
datKI <- data.frame(dat_norm[,c(1:3, 10:12)])
facKI <- fact[c(1:3, 10:12),]
str(datKI)
facKI$type <- as.factor(as.character(facKI$type))
facKI$type <- relevel(facKI$type, ref = "K")
facKI

Xki <- model.matrix(~ facKI$type)
Xki

fitki <- lmFit(datKI, design = Xki, method = "robust", maxit = 10000)

# Empirical Bayes statistics
efit_ki <- eBayes(fitki)

# Dif_expr_table
topTable(efit_ki, coef = 2)
full_list_ki <- topTable(efit_ki, number = length(dat_norm))
write.csv(full_list_ki,'Dif_expr_K_vs_I.csv')
head(full_list_ki)
#Vulcano
tiff('Vulcano_K_vs_I.tiff', units="in", width=11, height=8, res=600, compression = 'lzw')
EnhancedVolcano(full_list_ki,
                lab = rownames(full_list_ki),
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                xlim = c(-10, 11),
                ylim = c(0, 5),
                FCcutoff = 2,
                title ="K_vs_I",
                labSize = 4.0,
                boxedLabels = F,
                colAlpha = 1)
dev.off()

#С с I
datСI <- data.frame(dat_norm[,c(7:12)])
facСI <- fact[c(7:12),]
str(datСI)
facСI$type <- as.factor(as.character(facСI$type))
facСI$type

Xci <- model.matrix(~ facСI$type)
Xci

fitci <- lmFit(datСI, design = Xci, method = "robust", maxit = 10000)

# Empirical Bayes statistics
efit_ci <- eBayes(fitci)

# Dif_expr_table
topTable(efit_ci, coef = 2)
full_list_ci <- topTable(efit_ci, number = length(dat_norm))
write.csv(full_list_ci,'Dif_expr_c_vs_I_HCA.csv')
head(full_list_ci)
#Vulcano
tiff('Vulcano_c_vs_I_HCA.tiff', units="in", width=11, height=8, res=600, compression = 'lzw')
EnhancedVolcano(full_list_ci,
                lab = rownames(full_list_ci),
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                xlim = c(-10, 8),
                ylim = c(0, 4),
                FCcutoff = 2,
                title ="C_vs_I",
                labSize = 4.0,
                boxedLabels = F,
                colAlpha = 1)
dev.off()


#КM с CI

X_HCA <- model.matrix(~ fact$typeB)
X_HCA

fit_HCA <- lmFit(dat_norm, design = X_HCA, method = "robust", maxit = 10000)

# Empirical Bayes statistics
efit_HCA <- eBayes(fit_HCA)

# Dif_expr_table
topTable(efit_HCA, coef = 2)
full_list_HCA <- topTable(efit_HCA, number = length(dat_norm))
write.csv(full_list_HCA,'Dif_expr_KM_vs_CI.csv')
head(full_list_HCA)
#Vulcano
tiff('Vulcano_KM_vs_CI_HCA.tiff', units="in", width=11, height=8, res=600, compression = 'lzw')
EnhancedVolcano(full_list_HCA,
                lab = rownames(full_list_HCA),
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                xlim = c(-8, 9),
                ylim = c(0, 10),
                FCcutoff = 2,
                title ="K_vs_I",
                labSize = 4.0,
                boxedLabels = F,
                colAlpha = 1)
dev.off()




###############################################


##Готовим данные_HIT
#protein_expression_data

dat_HIT <- data.frame(read.csv("proteins_HIT.csv"))


dat_HIT$Accession <- sub("\\|.*", "", dat_HIT$Accession)                   # Extract first three characters

datdat_HIT1 <- dat_HIT[,c(3,18:29)]
head(datdat_HIT1)

str(datdat_HIT1)
rownames(datdat_HIT1) <- datdat_HIT1[,1]
datdat_HIT1 <- datdat_HIT1[,-1]
head(datdat_HIT1)



#делаем таблицу для количественного анализа
#по какому принципу я обрезаю датафрейм ниже?
dat_HIT2 <- datdat_HIT1[which(rowMeans(!is.na(datdat_HIT1)) >= 0.90), ]
mean(complete.cases(dat_HIT2))
colSums(is.na(dat_HIT2))

#Импутируем пропущенные значения по К-ближайшим соседям. Гуглим про это и читаем методичку :)
library(impute)
tdat_HIT <- t(dat_HIT2)
dat_knn_HIT1 <- impute.knn(tdat_HIT, k = 5)
dat_knn_HIT <- t(dat_knn_HIT1$data)
mean(complete.cases(dat_knn_HIT))

#Грузим нашу факторную матрицу (легенда к биологическим группам в образцах)
fact <- data.frame(c("sample", colnames(dat_HIT2)), c("type", "K", "K", "K", "M", "M", "M", "C", "C", "C", "I", "I", "I"), c("typecomb", "KM", "KM", "KM", "KM", "KM", "KM", "CI", "CI", "CI", "CI", "CI", "CI"))
colnames(fact) <- fact[1,]
fact <- fact[-1,]
str(fact)
rownames(fact) <- fact[,1]

fact$type <- as.factor(fact$type)

fact$type

fact$typecomb <- as.factor(fact$typecomb)
fact$typecomb
fact$typecomb <- relevel(fact$typecomb, ref = "KM")



#Контроль должен быть базовым уровнем !!! если нет, используем функцию relevel

#смотрим как выглядят данные и делаем их трансформацию - чтобы средние и распределения совпадали.
#сырые
pal <- brewer.pal(n = 9, name = "Set1")
cols <- pal[fact$type]
boxplot(dat_knn_HIT, outline = FALSE, col = cols, main = "Raw data")
legend("topright", levels(fact$type), fill = pal, bty = "n", xpd = T)
colSums(dat_knn_HIT)
#log transformation (погуглите зачем это делать). В чем разница с тем, что было? почему так удобно? зачем прибавлять 1 к dat_knn?
dat_log_HIT <- log2(dat_knn_HIT+1)
head(dat_log_HIT)
mean(complete.cases(dat_log_HIT))
boxplot(dat_log_HIT, outline = FALSE, col = cols, main = "Log-transformed data")
legend("topright", levels(fact$type), fill = pal, bty = "n", xpd = T)

#Quantile normalization
dat_norm_HIT <- normalizeQuantiles(dat_log_HIT)
head(dat_norm_HIT)
boxplot(dat_norm_HIT, col = cols, main = "Normalized data")
legend("topright", levels(fact$type), fill = pal, bty = "n", xpd = T)
mean(complete.cases(dat_norm_HIT))
colSums(is.na(dat_norm_HIT))


## Ординации
t_dat_HIT1 <- t(dat_norm_HIT[,-c(4:6)])
fact <- fact[-c(4:6),]
fact$type <- as.factor(as.character(fact$type))

library(mixOmics)
#PLS-DA 
ordination.optimum.splsda <- splsda(t_dat_HIT1, fact$type, ncomp = 3, keepX = c(15,15,15))
selectVar(ordination.optimum.splsda, comp=1)
selectVar(ordination.optimum.splsda, comp=2)
selectVar(ordination.optimum.splsda, comp=3)

tiff('PLSDA_HIT_three.tiff', units="in", width=16, height=8, res=600, compression = 'lzw')
layout(matrix(c(1, 2, 3, 3, 3, 3), 2, 3))
plotLoadings(ordination.optimum.splsda, comp = 1, size.name = 1, size.title = 1.2, title = "Loadings\n on 1st component", contrib = "max", legend = FALSE, col.ties="black", ndisplay = 15)
plotLoadings(ordination.optimum.splsda, comp = 2, size.name = 1, size.title = 1.2, title = "Loadings\n on 2nd component", contrib = "max",ndisplay = 15,  legend = FALSE, col.ties="black")
plotIndiv(ordination.optimum.splsda, ind.names = F, ellipse = T, style = "graphics", abline = TRUE, cex = 1.5, size.axis = 1.2, size.xlabel = 1.5, size.ylabel = 1.5, title = "PLS-DA ordination", size.title = 1.5, legend=TRUE)
dev.off()
layout(1,1)



## Ординации. Пока только PCA


dat_pca <- pca(t_dat_HIT1, ncomp = 8, center = TRUE)
dat_pca
plot(dat_pca)
dat_pca <- pca(t_dat_HIT1, ncomp = 3, center = TRUE)

plotIndiv(dat_pca, comp = c(1, 2), ind.names = F, 
          group = fact$type, legend = TRUE, ellipse = TRUE,
          title = 'PCA')

## Сохраняем. Что нового появилось на рисунке?
tiff('PCA_HIT_three.tiff', units="in", width=16, height=8, res=600, compression = 'lzw')
layout(matrix(c(1, 2, 3, 3, 3, 3), 2, 3))
plotLoadings(dat_pca, comp = 1, size.name = 1, size.title = 1.2, title = "Loadings\n on 1st component", contrib = "max", legend = FALSE, col.ties="black", ndisplay = 15)
plotLoadings(dat_pca, comp = 2, size.name = 1, size.title = 1.2, title = "Loadings\n on 2nd component", contrib = "max",ndisplay = 15,  legend = FALSE, col.ties="black")
plotIndiv(dat_pca, comp = c(1, 2), ind.names = F, 
          group = fact$type, legend = TRUE, ellipse = TRUE, style = "graphics", 
          title = 'PCA')
dev.off()
layout(1,1)

#КM с СI
X <- model.matrix(~ fact$typecomb)
X

fit_HIT <- lmFit(dat_norm_HIT, design = X, method = "robust", maxit = 10000)

# Empirical Bayes statistics
efit_HIT <- eBayes(fit_HIT)

# Dif_expr_table
topTable(efit_HIT, coef = 2)
full_list_HIT <- topTable(efit_HIT, number = length(dat_norm_HIT))
write.csv(full_list_HIT,'Dif_expr_KM_vs_CI_HIT.csv')
head(full_list_HIT)
#Vulcano
tiff('Vulcano_KM_vs_CI_HIT.tiff', units="in", width=11, height=8, res=600, compression = 'lzw')
EnhancedVolcano(full_list_HIT,
                lab = rownames(full_list_HIT),
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                xlim = c(-6, 8),
                ylim = c(0, 8),
                FCcutoff = 2,
                title ="K_vs_C",
                labSize = 4.0,
                boxedLabels = F,
                colAlpha = 1)
dev.off()



#К с I
dat_HIT_KI <- data.frame(dat_norm_HIT[,c(1:3, 10:12)])
facKI <- fact[c(1:3, 10:12),]
str(dat_HIT_KI)
facKI$type <- as.factor(as.character(facKI$type))
facKI$type <- relevel(facKI$type, ref = "K")
facKI$type

X_HIT_KI <- model.matrix(~ facKI$typecomb)
X_HIT_KI

fit_HIT_KI <- lmFit(dat_HIT_KI, design = X_HIT_KI, method = "robust", maxit = 10000)

# Empirical Bayes statistics
efit_HIT_KI <- eBayes(fit_HIT_KI)

# Dif_expr_table
topTable(efit_HIT_KI, coef = 2)
full_list_HIT_KI <- topTable(efit_HIT_KI, number = length(dat_norm_HIT))
write.csv(full_list_HIT_KI,'Dif_expr_K_vs_I_HIT.csv')
head(full_list_HIT_KI)
#Vulcano
tiff('Vulcano_K_vs_I_HIT.tiff', units="in", width=11, height=8, res=600, compression = 'lzw')
EnhancedVolcano(full_list_HIT_KI,
                lab = rownames(full_list_HIT_KI),
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                xlim = c(-12, 8),
                ylim = c(0, 5),
                FCcutoff = 2,
                title ="K_vs_I_HIT",
                labSize = 4.0,
                boxedLabels = F,
                colAlpha = 1)
dev.off()


#С vs I
dat_HIT_CI <- data.frame(dat_norm_HIT[,c(7:12)])
facCI <- fact[c(4:9),]
str(dat_HIT_CI)
facCI$type <- as.factor(as.character(facCI$type))
facCI$type <- relevel(facCI$type, ref = "C")
facCI$type


X_HIT_CI <- model.matrix(~ facCI$type)
X_HIT_CI

fit_HIT_CI <- lmFit(dat_HIT_CI, design = X_HIT_CI, method = "robust", maxit = 10000)

# Empirical Bayes statistics
efit_HIT_CI <- eBayes(fit_HIT_CI)

# Dif_expr_table
topTable(efit_HIT_CI, coef = 2)
full_list_HIT_CI <- topTable(efit_HIT_CI, number = length(dat_norm_HIT))
write.csv(full_list_HIT_CI,'Dif_expr_C_vs_I_HIT.csv')
head(full_list_HIT_CI)
#Vulcano
tiff('Vulcano_C_vs_I_HIT.tiff', units="in", width=11, height=8, res=600, compression = 'lzw')
EnhancedVolcano(full_list_HIT_CI,
                lab = rownames(full_list_HIT_CI),
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                xlim = c(-7, 8),
                ylim = c(0, 4),
                FCcutoff = 2,
                title ="C_vs_I",
                labSize = 4.0,
                boxedLabels = F,
                colAlpha = 1)
dev.off()

#К с С
dat_HIT_KC <- data.frame(dat_norm_HIT[,c(1:3, 7:9)])
facKC <- fact[c(1:3, 7:9),]
str(dat_HIT_KC)
facKC$type <- as.factor(as.character(facKC$type))
facKC$type <- relevel(facKC$type, ref = "K")
facKC$type


X_HIT_KC <- model.matrix(~ facKC$type)
X_HIT_KC

fit_HIT_KC <- lmFit(dat_HIT_KC, design = X_HIT_KC, method = "robust", maxit = 10000)

# Empirical Bayes statistics
efit_HIT_KC <- eBayes(fit_HIT_KC)

# Dif_expr_table
topTable(efit_HIT_KC, coef = 2)
full_list_HIT_KC <- topTable(efit_HIT_KC, number = length(dat_norm_HIT))
write.csv(full_list_HIT_KC,'Dif_expr_K_vs_C_HIT.csv')
head(full_list_HIT_KC)
#Vulcano
tiff('Vulcano_I_vs_C_HIT.tiff', units="in", width=11, height=8, res=600, compression = 'lzw')
EnhancedVolcano(full_list_HIT_KC,
                lab = rownames(full_list_HIT_KC),
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                xlim = c(-10, 10),
                ylim = c(0, 5),
                FCcutoff = 2,
                title ="K_vs_C",
                labSize = 4.0,
                boxedLabels = F,
                colAlpha = 1)
dev.off()

