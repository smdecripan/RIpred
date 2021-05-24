
#

load("gmdRi.rda")

db <- data.frame(matrix(vector(), nrow = length(gmdRi), ncol = 3, dimnames = list(c(), c("Name", "RI.VAR5.ALK", "MW"))))
for (i in 1:length(gmdRi@database)) {
  db[i,1] <- gmdRi@database[[i]]$Name
  db[i,2] <- as.numeric(gmdRi@database[[i]]$RI.VAR5.ALK)
  db[i,3] <- as.numeric(gmdRi@database[[i]]$MW)
}

load("fgSet/fgGmd_dragonTMS.rda")
dataTMS <- as.data.frame(matrix(unlist(fgGmd), 
                                ncol = ncol(fgGmd), nrow = nrow(fgGmd)),
                         stringsAsFactors = FALSE)
colnames(dataTMS) <- colnames(fgGmd)

db <- db[dataTMS$Name,]
db$Nrowdb <- 1:nrow(db)
db$Nameid <- rownames(db)

itN <- 1:20
randomCI <- NULL
for (i in itN) {
  set.seed(100 + i)
  randomCI[[i]] <- sort(sample(seq(1,nrow(db))[-which(db$MW == 0)], size = 70, replace = FALSE))
}
getPPM <- function(mass.error, mass){
  return((mass.error*10^6)/round(mass))
}

CI_ids <- NULL
for (i in itN) {
  x <- db[,"MW"]
  y = db[randomCI[[i]], "MW"]
  CI_ids[[i]] <- list()
  for (m in 1:length(y)) {
    CI_ids[[i]][[m]] <- which(sapply(x, function(mz) getPPM(abs(mz-y[m]),mz) < 10))
  }
}

dbSpectra <- NULL
for (i in as.numeric(db$Nameid)) {
  dbSpectra[[length(dbSpectra)+1]] <- matrix(as.numeric(unlist(strsplit(unlist(strsplit(x = gmdRi@database[[i]]$Spectra, split = " ")), ":"))), ncol = 2, byrow = TRUE, dimnames = list(NULL, c("mass", "intensity")))
}
randomEI <- NULL
for (i in itN) {
  set.seed(200 + i)
  randomEI[[i]] <- sort(sample(1:nrow(db), size = 20, replace = FALSE))
}

SpectMat <- matrix(70:600, ncol = 1, dimnames = list(NULL, "mass"))
for (i in 1:length(dbSpectra)) {
  namevect <- colnames(SpectMat)
  SpectMat <- merge(SpectMat, dbSpectra[[i]], by = "mass", all = TRUE)
  colnames(SpectMat) <- c(namevect, db$Nameid[i])
}
rownames(SpectMat) <- SpectMat$mass
SpectMat <- SpectMat[,-1]
SpectMat <- t(SpectMat)
SpectMat[,c(73:75,147:149)-69] <- 0
SpectMat[is.na(SpectMat)] <- 0

dot.product <- function(x,y){
  (x%*%t(y))^2/(rowSums(x^2)%*%t(rowSums(y^2)))
}
dotMat <- NULL
for (i in itN) {
  x <- SpectMat
  y <- SpectMat[randomEI[[i]],]
  dotMat[[i]] <- dot.product(x,y)
}

EI_ids <- NULL
for (i in itN) {
  EI_ids[[i]] <- apply(dotMat[[i]], 2, function(x){which(x >= 0.8)})
  
}

CIres <- NULL
EIres <- NULL
for (i in itN) {
  CIres[[i]] <- db[randomCI[[i]], -3]
  CIres[[i]]$Ncand <- unlist(lapply(CI_ids[[i]], length))
  CIres[[i]][,c(paste("Top", 1:5, sep = ""),paste("Top",1:5,"error", sep = ""))] <- NA
  EIres[[i]] <- db[randomEI[[i]], -3]
  EIres[[i]]$Ncand <- unname(unlist(lapply(EI_ids[[i]], length)))
  EIres[[i]][,c(paste("Top", 1:5, sep = ""),paste("Top",1:5,"error", sep = ""))] <- NA
}


#SVM linear model

library(e1071)
resSVM_CI <- list()
ModelSum_CI <- data.frame(matrix(vector(), ncol = 5,
                                 dimnames = list(NULL, c("Model",
                                                         "MAE", "MAPE",
                                                         "MdAE", "MdAPE"
                                 ))))
for (i in itN) {
  train_data <- dataTMS[-unique(unlist(CI_ids[[i]])), 3:ncol(dataTMS)]
  test_data <- dataTMS[unique(unlist(CI_ids[[i]])), 3:ncol(dataTMS)]
  #To facilitate models predicting arbitrary values lets transform labels, 
  #expressed in seconds (?), multiplying them by an order of 10:
  train_labels <- dataTMS[-unique(unlist(CI_ids[[i]])), 2]*10
  test_labels <- dataTMS[unique(unlist(CI_ids[[i]])), 2]*10
  
  ##SVM
  #linear kernel
  set.seed(123)
  linear.SVM <- NULL
  linear.SVM <- svm(x = train_data, y = train_labels,
                    kernel = "linear",
                    type = "eps-regression",
                    cost = 100,
                    epsilon = 0.1)
  testPred <- predict(linear.SVM, newdata = test_data)
  resSVM_CI[[i]] <- data.frame("EmpiricalRI" = test_labels,
                               "PredictedRI" = testPred, 
                               "AE" = c(abs(testPred - test_labels)),
                               "APE" = c(abs((testPred - test_labels) / test_labels)))
  
  ModelSum_CI[i,] = c("SVMlinear", mean(resSVM_CI[[i]]$AE), mean(resSVM_CI[[i]]$APE),
                      median(resSVM_CI[[i]]$AE), median(resSVM_CI[[i]]$APE)
  )
  
  
}

resSVM_EI <- list()
ModelSum_EI <- data.frame(matrix(vector(), ncol = 5,
                                 dimnames = list(NULL, c("Model",
                                                         "MAE", "MAPE",
                                                         "MdAE", "MdAPE"
                                 ))))
for (i in itN) {
  train_data <- dataTMS[-unique(unlist(EI_ids[[i]])), 3:ncol(dataTMS)]
  test_data <- dataTMS[unique(unlist(EI_ids[[i]])), 3:ncol(dataTMS)]
  #To facilitate models predicting arbitrary values lets transform labels, 
  #expressed in seconds (?), multiplying them by an order of 10:
  train_labels <- dataTMS[-unique(unlist(EI_ids[[i]])), 2]*10
  test_labels <- dataTMS[unique(unlist(EI_ids[[i]])), 2]*10
  
  ##SVM
  #linear kernel
  set.seed(123)
  linear.SVM <- NULL
  linear.SVM <- svm(x = train_data, y = train_labels,
                    kernel = "linear",
                    type = "eps-regression",
                    cost = 100,
                    epsilon = 0.1)
  testPred <- predict(linear.SVM, newdata = test_data)
  resSVM_EI[[i]] <- data.frame("EmpiricalRI" = test_labels,
                               "PredictedRI" = testPred,
                               "AE" = c(abs(testPred - test_labels)),
                               "APE" = c(abs((testPred - test_labels) / test_labels)))
  
  ModelSum_EI[i,] = c("SVMlinear", mean(resSVM_EI[[i]]$AE), mean(resSVM_EI[[i]]$APE),
                      median(resSVM_EI[[i]]$AE), median(resSVM_EI[[i]]$APE)
  )
  
  
}

resSVM_CItransf <- resSVM_CI
resSVM_EItransf <- resSVM_EI
for (i in itN) {
  resSVM_CItransf[[i]][,1:3] <- resSVM_CItransf[[i]][,1:3]/10
  resSVM_EItransf[[i]][,1:3] <- resSVM_EItransf[[i]][,1:3]/10
}

for (i in itN) {
  for (m in 1:nrow(CIres[[i]])) {
    ids <- CI_ids[[i]][[m]]
    x <- resSVM_CItransf[[i]][as.character(ids),]
    errVect <- abs(CIres[[i]][m,2]-x$PredictedRI)/CIres[[i]][m,2]
    
    ids <- ids[1:5]
    CIres[[i]][m,6:(5+length(ids))] <- rownames(x)[order(errVect)][1:5]
    CIres[[i]][m, 11:(10+length(ids))] <- sort(errVect)[1:5]
  }
}

for (i in itN) {
  for (m in 1:nrow(EIres[[i]])) {
    ids <- EI_ids[[i]][[m]]
    x <- resSVM_EItransf[[i]][as.character(ids),]
    errVect <- abs(EIres[[i]][m,2]-x$PredictedRI)/EIres[[i]][m,2]
    
    ids <- ids[1:5]
    EIres[[i]][m,6:(5+length(ids))] <- rownames(x)[order(errVect)][1:5]
    EIres[[i]][m, 11:(10+length(ids))] <- sort(errVect)[1:5]
  }
}


EI.pie <- data.frame(Ncand = c(rep(2,4), rep(3,4), rep(4,5)),
                     Top = as.factor(c(rep(c(1,2,3,"T"),3), "Other")),
                     Value = 0, Percent = NA)
for (i in itN) {
  EI.pie[1,3] <- sum(EI.pie[1,3], nrow(subset(EIres[[i]], Ncand <= 2 & Top1 == Nrowdb)))
  EI.pie[2,3] <- sum(EI.pie[2,3], nrow(subset(EIres[[i]], Ncand <= 2 & Top2 == Nrowdb)))
  EI.pie[3,3] <- sum(EI.pie[3,3], nrow(subset(EIres[[i]], Ncand <= 2 & Top3 == Nrowdb)))
  EI.pie[4,3] <- sum(EI.pie[4,3], nrow(subset(EIres[[i]], Ncand <= 2)))
  EI.pie[5,3] <- sum(EI.pie[5,3], nrow(subset(EIres[[i]], Ncand == 3 & Top1 == Nrowdb)))
  EI.pie[6,3] <- sum(EI.pie[6,3], nrow(subset(EIres[[i]], Ncand == 3 & Top2 == Nrowdb)))
  EI.pie[7,3] <- sum(EI.pie[7,3], nrow(subset(EIres[[i]], Ncand == 3 & Top3 == Nrowdb)))
  EI.pie[8,3] <- sum(EI.pie[8,3], nrow(subset(EIres[[i]], Ncand == 3)))
  EI.pie[9,3] <- sum(EI.pie[9,3], nrow(subset(EIres[[i]], Ncand >= 4 & Top1 == Nrowdb)))
  EI.pie[10,3] <- sum(EI.pie[10,3], nrow(subset(EIres[[i]], Ncand >= 4 & Top2 == Nrowdb)))
  EI.pie[11,3] <- sum(EI.pie[11,3], nrow(subset(EIres[[i]], Ncand >= 4 & Top3 == Nrowdb)))
  EI.pie[12,3] <- sum(EI.pie[12,3], nrow(subset(EIres[[i]], Ncand >= 4)))
  EI.pie[13,3] <- EI.pie[12,3]-sum(EI.pie[9:11,3])
}
EI.pie
for (i in c(1,5,9)) {
  EI.pie[i,4] <- EI.pie[i,3]/EI.pie[i+3,3]*100
  EI.pie[i+1,4] <- EI.pie[i+1,3]/EI.pie[i+3,3]*100
  EI.pie[i+2,4] <- EI.pie[i+2,3]/EI.pie[i+3,3]*100
  EI.pie[i+3,4] <- EI.pie[i+3,3]/EI.pie[i+3,3]*100
}
EI.pie[13,4] <- EI.pie[13,3]/EI.pie[12,3]*100
EI.pie$Percent <- round(EI.pie$Percent, digits = 2)
EI.pie

CI.pie <- data.frame(Ncand = c(rep(2,4), rep(3,4), rep(4,5)),
                     Top = as.factor(c(rep(c(1,2,3,"T"),3), "Other")),
                     Value = 0, Percent = NA)
for (i in itN) {
  CI.pie[1,3] <- sum(CI.pie[1,3], nrow(subset(CIres[[i]], Ncand <= 2 & Top1 == Nrowdb)))
  CI.pie[2,3] <- sum(CI.pie[2,3], nrow(subset(CIres[[i]], Ncand <= 2 & Top2 == Nrowdb)))
  CI.pie[3,3] <- sum(CI.pie[3,3], nrow(subset(CIres[[i]], Ncand <= 2 & Top3 == Nrowdb)))
  CI.pie[4,3] <- sum(CI.pie[4,3], nrow(subset(CIres[[i]], Ncand <= 2)))
  CI.pie[5,3] <- sum(CI.pie[5,3], nrow(subset(CIres[[i]], Ncand == 3 & Top1 == Nrowdb)))
  CI.pie[6,3] <- sum(CI.pie[6,3], nrow(subset(CIres[[i]], Ncand == 3 & Top2 == Nrowdb)))
  CI.pie[7,3] <- sum(CI.pie[7,3], nrow(subset(CIres[[i]], Ncand == 3 & Top3 == Nrowdb)))
  CI.pie[8,3] <- sum(CI.pie[8,3], nrow(subset(CIres[[i]], Ncand == 3)))
  CI.pie[9,3] <- sum(CI.pie[9,3], nrow(subset(CIres[[i]], Ncand >= 4 & Top1 == Nrowdb)))
  CI.pie[10,3] <- sum(CI.pie[10,3], nrow(subset(CIres[[i]], Ncand >= 4 & Top2 == Nrowdb)))
  CI.pie[11,3] <- sum(CI.pie[11,3], nrow(subset(CIres[[i]], Ncand >= 4 & Top3 == Nrowdb)))
  CI.pie[12,3] <- sum(CI.pie[12,3], nrow(subset(CIres[[i]], Ncand >= 4)))
  CI.pie[13,3] <- CI.pie[12,3]-sum(CI.pie[9:11,3])
}

CI.pie
for (i in c(1,5,9)) {
  CI.pie[i,4] <- CI.pie[i,3]/CI.pie[i+3,3]*100
  CI.pie[i+1,4] <- CI.pie[i+1,3]/CI.pie[i+3,3]*100
  CI.pie[i+2,4] <- CI.pie[i+2,3]/CI.pie[i+3,3]*100
  CI.pie[i+3,4] <- CI.pie[i+3,3]/CI.pie[i+3,3]*100
}
CI.pie[13,4] <- CI.pie[13,3]/CI.pie[12,3]*100
CI.pie$Percent <- round(CI.pie$Percent, digits = 2)
CI.pie

PieC2 <- ggplot(CI.pie[1:3,2:3], aes(x="", y=Value, fill=Top)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="none")+
  #geom_text(aes(y = Value, label = Value), color = "white", size=6) 
  scale_fill_manual(values=c("#9ebcda", "#8c96c6", "#8c6bb1", "#88419d"))#+
#annotate(geom="text", x=2, y=0, label='bold("CI")', parse = TRUE, size = 5)
PieC3 <- ggplot(CI.pie[5:7,2:3], aes(x="", y=Value, fill=Top)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="none")+
  #geom_text(aes(y = Value, label = Value), color = "white", size=6) 
  scale_fill_manual(values=c("#9ebcda", "#8c96c6", "#8c6bb1", "#88419d"))
PieC4 <- ggplot(CI.pie[c(9:11,13),2:3], aes(x="", y=Value, fill=Top)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  #theme(legend.position="bottom") +
  theme(legend.position="none")+
  #geom_text(aes(y = Value, label = Value), color = "white", size=6) 
  scale_fill_manual(values=c("#9ebcda", "#8c96c6", "#8c6bb1", "#88419d"))

PieE2 <- ggplot(EI.pie[1:3,2:3], aes(x="", y=Value, fill=Top)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="none")+
  #geom_text(aes(y = Value, label = Value), color = "white", size=6) 
  scale_fill_manual(values=c("#9ebcda", "#8c96c6", "#8c6bb1", "#88419d"))
PieE3 <- ggplot(EI.pie[5:7,2:3], aes(x="", y=Value, fill=Top)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="none")+
  #geom_text(aes(y = Value, label = Value), color = "white", size=6) 
  scale_fill_manual(values=c("#9ebcda", "#8c96c6", "#8c6bb1", "#88419d"))
PieE4 <- ggplot(EI.pie[c(9:11,13),2:3], aes(x="", y=Value, fill=Top)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="none")+
  #theme(legend.position="bottom", legend.key.size = unit(2.5, 'cm'),
  #legend.title = element_text(size=22, face = "bold"), legend.text = element_text(size = 20)) +
  #geom_text(aes(y = Value, label = Value), color = "white", size=6) 
  scale_fill_manual(values=c("#9ebcda", "#8c96c6", "#8c6bb1", "#88419d"))

PiePlot <- ggarrange(PieC2, PieE2, PieC3, PieE3, PieC4, PieE4,
                     labels = c("A", "B"),
                     font.label = list(size = 22, color = "black", face = "bold", family = NULL),
                     ncol = 2, nrow = 3, common.legend = TRUE
)
PiePlot


CI_all <- CIres[[1]][0,]
EI_all <- EIres[[1]][0,]
for (i in itN) {
  CI_all <- rbind(CI_all, CIres[[i]])
  EI_all <- rbind(EI_all, EIres[[i]])
}

roc <- function(data = EI_all, ThresVect = seq(0.01,1,0.01)){
  rocData <- data.frame(ThresVect, "TP"=NA, "TN"=NA, "FP"=NA, "FN"=NA, 
                        "Sens"=NA, "Spec"=NA, "FDR"=NA)
  a <- data[,c("Nrowdb","Top1", "Top2", "Top3", "Top1error", "Top2error", "Top3error")]
  
  for (t in 1:nrow(rocData)) {
    a$Top1class <- NA
    a$Top2class <- NA
    a$Top3class <- NA
    Thres <- rocData[t,1]
    for (n in 1:nrow(a)) {
      for (i in 2:4) {
        b <- as.numeric(a[n,i])
        if (is.na(b)) {next} else {#nrow of the candidate (Top1/2/3)
          e <- a[n,(i+3)] <= Thres #TRUE/FALSE
          c <- b == a[n,1] #TRUE/FALSE (match)
          if (e == TRUE) { #POSITIVES
            a[n,i+6] <- ifelse(c, "TP", "FP")
          } else {#NEGATIVES
            a[n,i+6] <- ifelse(c, "FN", "TN")
          }
        }
      }
    }
    rocData[t,'TP'] <- table(c(a[,8],a[,9],a[,10]))['TP']
    rocData[t,'TN'] <- table(c(a[,8],a[,9],a[,10]))['TN']
    rocData[t,'FP'] <- table(c(a[,8],a[,9],a[,10]))['FP']
    rocData[t,'FN'] <- table(c(a[,8],a[,9],a[,10]))['FN']
    rocData[t,"Sens"] <- rocData[t,"TP"]/sum(rocData[t,c("TP","FN")])
    rocData[t,"Spec"] <- rocData[t,"TN"]/sum(rocData[t,c("TN","FP")])
    rocData[t, "FDR"] <- rocData[t,"FP"]/sum(rocData[t,c("TP","FP")])
  }
  return(rocData)
}

rocDataEI <- roc(data = EI_all, ThresVect = seq(0,1,0.01))
rocDataEI[1,c(6,7)] <- c(0,1)

rocEI <- rocDataEI[-which(is.na(rocDataEI$Sens)),]
inForAUC <- approx(c(1-rocEI$Spec, 1), y = c(rocEI$Sens,1), xout = seq(0,1, by = 0.01))$y
AUC.EI <- sum(inForAUC)/length(inForAUC)
AUC.EI

plot(1-rocDataEI[,'Spec'],rocDataEI[,'Sens'], type = 'l', 
     xlab = "FPR", ylab = "TPR", 
     xlim = c(0,1), ylim = c(0,1),
     xaxs="i",
     yaxs="i",
     main = "AUC: 0.56", font.main = 1)
abline(0,1, lty = 2)

rocDataCI <- roc(data = CI_all, ThresVect = seq(0,1,0.01))
rocDataCI[1,c(6,7)] <- c(0,1)

rocCI <- rocDataCI[-which(is.na(rocDataCI$Sens)),]
inForAUC <- approx(c(1, 1-rocCI$Spec,0),c(1, rocCI$Sens,0), xout = seq(0,1, by = 0.01))$y
AUC.CI <- sum(inForAUC)/length(inForAUC)
AUC.CI

plot(1-rocDataCI[,'Spec'],rocDataCI[,'Sens'], type = 'l', 
     xlab = "FPR", ylab = "TPR", 
     xlim = c(0,1), ylim = c(0,1),
     xaxs="i",
     yaxs="i",
     main = "AUC: 0.53", font.main = 1)
abline(0,1, lty = 2)



barDataEI <- data.frame(Error = c(EI_all$Top1error, EI_all$Top2error, EI_all$Top3error),
                        Total = as.factor(rep("T", nrow(EI_all)*3))
)
barDataEI$Filtered <- barDataEI$Error <= 0.03
barDataEI <- data.frame(Name = c("T", "F", rep("C", 4)),
                        Variable = c("T", "F", "TP", "FP", "TN", "FN"),
                        Value = c(length(which(c(EI_all$Top1error, EI_all$Top2error, EI_all$Top3error) != "NA")),
                                  length(which(c(EI_all$Top1error, EI_all$Top2error, EI_all$Top3error) <= 0.03)),
                                  rocDataEI[4,'TP'],
                                  rocDataEI[4,'FP'],
                                  rocDataEI[4,'TN'],
                                  rocDataEI[4,'FN']
                        )
)
barDataEI$Name <- factor(barDataEI$Name, levels = c("T", "F", "C"))

barDataCI <- data.frame(Name = c("T", "F", rep("C", 4)),
                        Variable = c("T", "F", "TP", "FP", "TN", "FN"),
                        Value = c(length(which(c(CI_all$Top1error, CI_all$Top2error, CI_all$Top3error) != "NA")),
                                  length(which(c(CI_all$Top1error, CI_all$Top2error, CI_all$Top3error) <= 0.03)),
                                  rocDataCI[4,'TP'],
                                  rocDataCI[4,'FP'],
                                  rocDataCI[4,'TN'],
                                  rocDataCI[4,'FN']
                        )
)
barDataCI$Name <- factor(barDataCI$Name, levels = c("T", "F", "C"))

library(ggplot2)
library(ggpubr)
barEI <- ggplot(data = barDataEI, aes(y = Value, x = Name, 
                                      fill = factor(Variable, levels = c("T","F" ,"TP", "FP","TN", "FN")))) +
  geom_bar(stat = "identity", colour = "black") +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        plot.title = element_text(color="black", size=14, face="bold.italic", hjust = 0.5)
  )+ ggtitle("EI")+
  xlab("") + ylab("") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 20), axis.text.y = element_text(size = 15)) +
  scale_fill_manual(name = "",
                    breaks = c("TP", "FP","TN", "FN"),
                    values = c("#2166ac", "#92c5de", "#fee060", "#fee090", "#f4a582", "#f46d43"))
barCI <- ggplot(data = barDataCI, aes(y = Value, x = Name, 
                                      fill = factor(Variable, levels = c("T","F" ,"TP", "FP","TN", "FN")))) +
  geom_bar(stat = "identity", colour = "black") +
  theme_bw() + ggtitle("CI")+
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        plot.title = element_text(color="black", size=14, face="bold.italic", hjust = 0.5)
  )+ 
  xlab("") + ylab("") + #theme(legend.position="none") 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 20), axis.text.y = element_text(size = 15)) +
  scale_fill_manual(name = "",
                    breaks = c("TP", "FP","TN", "FN"),
                    values = c("#2166ac", "#92c5de", "#fee060", "#fee090", "#f4a582", "#f46d43"))

BarSim <- ggarrange(barCI, barEI, 
                    labels = c("C", "D"),
                    font.label = list(size = 32, color = "black", face = "bold", family = NULL),
                    ncol = 2, nrow = 1,
                    common.legend = TRUE)
BarSim








