###MACHINE LEARNING MODELS TRAINING AND TESTING WITH DIFFERENT FINGERPRINT CLASSES
#Import files containing diferent fingerprint molecular representations
file.dir <- "fgSet/"
file.list <- paste(file.dir , list.files(file.dir), sep = "")


df <- list()
for (i in 1:length(file.list)) {
  load(file.list[i])
  df[[i]] <- matrix(unlist(fgGmd), ncol = ncol(fgGmd), nrow = nrow(fgGmd))
}
names(df) <- c("TMSdragon", "TMSecfp2", "TMSecfp4", "TMSlayered",
               "TMSmaccs", "TMSobfp2", "TMSpb")

rowsN <- as.numeric(unique(unlist(lapply(df, nrow))))


##Preparing data sets
###Data will be randomly split as: train set (75% = 869), and test set (25% = 290) 
###but due to the limited number of compounds each model will be tested iteratively 20 times.
itN <- 1:20
random_ids <- list()
for (i in itN) {
  set.seed(120+i)
  random_ids[[i]] <- sample(rowsN, size = 290, replace = FALSE)
}

###With 20 iterations for data split we ensure that almost every compound is tested in 
###both training and test sets. It will be useful to assess whether the degree of similarity 
###between molecular structures has any effect on the predicted values.

idNames <- lapply(random_ids, function(x) df[[1]][x,1])

#Training the models
library(keras)
library(e1071)
library(randomForest)

##DNN model
model <- NULL
set.seed(123)
build_model <- NULL
build_model <- function(){
  model <- NULL
  model <- keras_model_sequential() %>%
    layer_dense(units = unitsLayer1, activation = "relu",
                input_shape = ncol(train_data)) %>%
    layer_dense(units = unitsLayer2, activation = "relu") %>%
    layer_dense(units = 1)
  model %>% compile(
    optimizer = "rmsprop",
    loss = "mse",
    metrics = c("mae")
  )
}
##DNN parameters
unitsLayer1 <- 200
unitsLayer2 <- 64
epochsDNN <- 30
batchSizeDNN <- 1

##CNN model
set.seed(1234)
build_model_cnn <- NULL
build_model_cnn <- function(){
  model <- NULL
  model <- keras_model_sequential() %>%
    layer_conv_1d(filters = filtersLayer1, kernel_size = kernelSize,
                  activation = "relu",
                  input_shape = list(ncol(train_data), 1)) %>%
    layer_average_pooling_1d(pool_size = poolSize) %>%
    layer_conv_1d(filters = filtersLayer2, kernel_size = kernelSize,
                  activation = "relu") %>%
    layer_average_pooling_1d(pool_size = poolSize) %>%
    layer_conv_1d(filters = filtersLayer3, kernel_size = kernelSize, activation = "relu") %>%
    layer_average_pooling_1d(pool_size = poolSize) %>%
    layer_flatten()%>%
    layer_dense(units = unitsLayer4, activation = "relu", ) %>%
    layer_dense(units = 1, activation = "relu")
  
  model %>% compile(
    optimizer = "rmsprop",
    loss = "mse",
    metrics = c("mae")
  )
}

##Parameterns CNN
kernelSize <- 6
poolSize <- 4
filtersLayer1 <- 600
filtersLayer2 <- 320
filtersLayer3 <- 164
unitsLayer4 <- 64
epochsCNN <- 130
batchSizeCNN <- 1


resDNN <- list()
resSVMinear <- list()
resSVMpoly <- list()
resRF <- list()
resCNN <- list()
ModelSum <- list()
for (f in 1) { #Each model will be trained and evaluated with dragon fingerprints
  ##DNN
  resDNN[[f]] <- list()
  resSVMinear[[f]] <- list()
  resSVMpoly[[f]] <- list()
  resRF[[f]] <- list()
  resCNN[[f]] <- list()
  ModelSum[[f]] <- list()
  ModelSum[[f]]$DNN <- data.frame(matrix(vector(), ncol = 5,
                                         dimnames = list(NULL, c("Model",
                                                                 "MAE", "MAPE",
                                                                 "MdAE", "MdAPE"
                                         ))))
  ModelSum[[f]]$SVMlinear <- data.frame(matrix(vector(), ncol = 5,
                                               dimnames = list(NULL, c("Model",
                                                                       "MAE", "MAPE",
                                                                       "MdAE", "MdAPE"
                                               ))))
  ModelSum[[f]]$SVMpoly <- data.frame(matrix(vector(), ncol = 5,
                                             dimnames = list(NULL, c("Model",
                                                                     "MAE", "MAPE",
                                                                     "MdAE", "MdAPE"
                                             ))))
  ModelSum[[f]]$RandomForest <- data.frame(matrix(vector(), ncol = 5,
                                                  dimnames = list(NULL, c("Model",
                                                                          "MAE", "MAPE",
                                                                          "MdAE", "MdAPE"
                                                  ))))
  ModelSum[[f]]$CNN <- data.frame(matrix(vector(), ncol = 5,
                                         dimnames = list(NULL, c("Model",
                                                                 "MAE", "MAPE",
                                                                 "MdAE", "MdAPE"
                                         ))))
  for (i in itN) {
    train_data <- df[[f]][-random_ids[[i]], 3:ncol(df[[f]])]
    test_data <- df[[f]][random_ids[[i]], 3:ncol(df[[f]])]
    #To facilitate models predicting arbitrary values lets transform labels,
    #multiplying them by an order of 10:
    train_labels <- df[[f]][-random_ids[[i]], 2]*10
    test_labels <- df[[f]][random_ids[[i]], 2]*10
    
    
    ##DNN
    
    modeldnn = NULL
    modeldnn <- build_model()
    #historydnn <- NULL
    historydnn <- modeldnn %>% fit(train_data, train_labels, 
                                   epochs = epochsDNN, batch_size = batchSizeDNN, verbose = 0)
    modEval <- modeldnn %>% evaluate(test_data, test_labels)
    testPred <- modeldnn %>% predict(test_data)
    resDNN[[f]][[i]] <- data.frame("EmpiricalRI" = test_labels,
                                   "PredictedRI" = testPred, 
                                   "AE" = c(abs(testPred - test_labels)),
                                   "APE" = c(abs((testPred - test_labels) / test_labels)))
    ModelSum[[f]]$DNN[i,] = c("DNN", mean(resDNN[[f]][[i]]$AE), mean(resDNN[[f]][[i]]$APE),
                              median(resDNN[[f]][[i]]$AE), median(resDNN[[f]][[i]]$APE)
    )
    
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
    resSVMinear[[f]][[i]] <- data.frame("EmpiricalRI" = test_labels,
                                        "PredictedRI" = testPred, 
                                        "AE" = c(abs(testPred - test_labels)),
                                        "APE" = c(abs((testPred - test_labels) / test_labels)))
    
    ModelSum[[f]]$SVMlinear[i,] = c("SVMlinear", mean(resSVMinear[[f]][[i]]$AE), mean(resSVMinear[[f]][[i]]$APE),
                                    median(resSVMinear[[f]][[i]]$AE), median(resSVMinear[[f]][[i]]$APE)
    )
    #polynomial kernel
    set.seed(123)
    poly.SVM <- svm(x = train_data, y = train_labels,
                    kernel = "polynomial",
                    type = "eps-regression",
                    cost = 1,
                    degree = 6,
                    coef0 = 6,
                    epsilon = 0.1)
    testPred <- predict(poly.SVM, newdata = test_data)
    resSVMpoly[[f]][[i]] <- data.frame("EmpiricalRI" = test_labels,
                                       "PredictedRI" = testPred,
                                       "AE" = c(abs(testPred - test_labels)),
                                       "APE" = c(abs((testPred - test_labels) / test_labels)))
    ModelSum[[f]]$SVMpoly[i,] = c("SVMpoly", mean(resSVMpoly[[f]][[i]]$AE), mean(resSVMpoly[[f]][[i]]$APE),
                                  median(resSVMpoly[[f]][[i]]$AE), median(resSVMpoly[[f]][[i]]$APE)
    )
    
    ##Random Forest
    set.seed(123)
    RFmodel <- randomForest(x = train_data, y = train_labels, 
                            ntree = 600, importance = FALSE)
    testPred <- predict(RFmodel, test_data)
    resRF[[f]][[i]] <- data.frame("EmpiricalRI" = test_labels,
                                  "PredictedRI" = testPred,
                                  "AE" = c(abs(testPred - test_labels)),
                                  "APE" = c(abs((testPred - test_labels) / test_labels)))
    ModelSum[[f]]$RandomForest[i,] = c("SVMlinear", mean(resRF[[f]][[i]]$AE), mean(resRF[[f]][[i]]$APE),
                                       median(resRF[[f]][[i]]$AE), median(resRF[[f]][[i]]$APE)
    )
    
    ##CNN
    x_train <- array_reshape(train_data, dim= c(nrow(train_data), ncol(train_data), 1))
    x_test <- array_reshape(test_data, dim = c(nrow(test_data), ncol(test_data), 1))
    
    modelcnn = NULL
    modelcnn <- build_model_cnn()
    historycnn <- NULL
    historycnn <- modelcnn %>%  fit(x_train, train_labels,
                                    epochs = epochsCNN,
                                    batch_size = batchSizeCNN)
    modEval <- modelcnn %>% evaluate(x_test, test_labels)
    testPred <- modelcnn %>% predict(x_test)
    resCNN[[f]][[i]] <- data.frame("EmpiricalRI" = test_labels,
                                   "PredictedRI" = testPred, 
                                   "AE" = c(abs(testPred - test_labels)),
                                   "APE" = c(abs((testPred - test_labels) / test_labels)))
    ModelSum[[f]]$CNN[i,] = c("CNN", mean(resCNN[[f]][[i]]$AE), mean(resCNN[[f]][[i]]$APE),
                              median(resCNN[[f]][[i]]$AE), median(resCNN[[f]][[i]]$APE)
    )
    
  }
  
}





for (f in 2:length(df)) { #SVM model with linear kernel was evaluated with remaining fingerprint classes
  resSVMinear[[f]] <- list()
  ModelSum[[f]] <- list()
  ModelSum[[f]]$SVMlinear <- data.frame(matrix(vector(), ncol = 5,
                                               dimnames = list(NULL, c("Model",
                                                                       "MAE", "MAPE",
                                                                       "MdAE", "MdAPE"
                                               ))))
  for (i in itN) {
    train_data <- df[[f]][-random_ids[[i]], 3:ncol(df[[f]])]
    test_data <- df[[f]][random_ids[[i]], 3:ncol(df[[f]])]
    train_labels <- df[[f]][-random_ids[[i]], 2]*10
    test_labels <- df[[f]][random_ids[[i]], 2]*10
    
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
    resSVMinear[[f]][[i]] <- data.frame("EmpiricalRI" = test_labels,
                                        "PredictedRI" = testPred, 
                                        "AE" = c(abs(testPred - test_labels)),
                                        "APE" = c(abs((testPred - test_labels) / test_labels)))
    
    ModelSum[[f]]$SVMlinear[i,] = c("SVMlinear", mean(resSVMinear[[f]][[i]]$AE), mean(resSVMinear[[f]][[i]]$APE),
                                    median(resSVMinear[[f]][[i]]$AE), median(resSVMinear[[f]][[i]]$APE)
    )
  }
}

#Prediction error for each ML model with dragon fingerprints as input data
boxplot(unlist(lapply(resSVMinear[[1]], function(x){x$APE*100})),
        unlist(lapply(resSVMpoly[[1]], function(x){x$APE*100})), 
        unlist(lapply(resDNN[[1]], function(x){x$APE*100})),
        unlist(lapply(resCNN[[1]], function(x){x$APE*100})),
        unlist(lapply(resRF[[1]], function(x){x$APE*100})),
        
        outline = FALSE, col = c("#6c5d53", "#917c6f", "#ac9d93", "#c8beb7", "#e3dedb"), 
        axes = FALSE, main = "",
        boxlty = 0, whisklty = 3, whisklwd = 2,
        staplelwd = 2, ylim = c(0,20), ylab = "Relative error (%)"
        )

axis(2, seq(0,20, 5), seq(0, 20, 5), lwd = 3, las = 1)
text(1:5, 
     par("usr")[3], 
     srt = 45, adj = c(1, 3), xpd = TRUE,
     labels = c("SVM-Lin", "SVM-Poly", "DNN", "CNN", "RF"), 
     cex = 1.5)

#Prediction error for each FGPT class. 
boxplot(SVMlinearAE[[1]]$APE*100,
        SVMlinearAE[[7]]$APE*100,
        SVMlinearAE[[11]]$APE*100,
        SVMlinearAE[[3]]$APE*100,
        SVMlinearAE[[5]]$APE*100,
        SVMlinearAE[[9]]$APE*100,
        SVMlinearAE[[13]]$APE*100,
        outline = FALSE, col = c("#49ad84", "#b09e4b", "#cb994e", "#d7944e", 
                                 "#a7b1ec", "#e692c3", "#6daee2"), 
        axes = FALSE, main = "",
        boxlty = 0, whisklty = 3, whisklwd = 2,
        staplelwd = 2,
        ylim = c(0,50), ylab = "Relative error (%)"
        )

axis(2, seq(0,50, 5), seq(0, 50, 5), lwd = 3, las = 1)
text(1:7, 
     par("usr")[3], 
     srt = 45, adj = c(1, 3), xpd = TRUE,
     labels = c("Dagon7", "Layered", "OBFP2", "ECFP2", "ECFP4", "MACCS", "PB"), 
     cex = 1.5)


MLdata <- data.frame(groups = rep(c("SVMlinear", "SVMpoly", "DNN", "CNN", "RF"), each = 5800),
                     APE = c(unlist(lapply(resSVMinear[[1]], function(x){x$APE*100})),
                             unlist(lapply(resSVMpoly[[1]], function(x){x$APE*100})),
                             unlist(lapply(resDNN[[1]], function(x){x$APE*100})),
                             unlist(lapply(resCNN[[1]], function(x){x$APE*100})),
                             unlist(lapply(resRF[[1]], function(x){x$APE*100}))))

##Paired Wilcox test across ML models
wilcox.test(APE ~ groups, 
            data = MLdata[c(which(MLdata$groups == "SVMlinear"),which(MLdata$groups == "SVMpoly")),], 
            paired = TRUE)

wilcox.test(APE ~ groups, 
            data = MLdata[c(which(MLdata$groups == "SVMlinear"),which(MLdata$groups == "DNN")),], 
            paired = TRUE)

wilcox.test(APE ~ groups, 
            data = MLdata[c(which(MLdata$groups == "SVMlinear"),which(MLdata$groups == "CNN")),], 
            paired = TRUE)

wilcox.test(APE ~ groups, 
            data = MLdata[c(which(MLdata$groups == "SVMlinear"),which(MLdata$groups == "RF")),], 
            paired = TRUE)

FPdata <- data.frame(groups = rep(names(SVMlinearAE), each = 5800),
                     APE = c(SVMlinearAE[[1]]$APE*100,
                             SVMlinearAE[[2]]$APE*100,
                             SVMlinearAE[[3]]$APE*100,
                             SVMlinearAE[[4]]$APE*100,
                             SVMlinearAE[[5]]$APE*100,
                             SVMlinearAE[[6]]$APE*100,
                             SVMlinearAE[[7]]$APE*100))
##Paired Wilcox test across Fingerprint classes
wilcox.test(APE ~ groups, 
            data = FPdata[c(1:5800,which(FPdata$groups == "TMSecfp2")),], paired = TRUE)

wilcox.test(APE ~ groups, 
            data = FPdata[c(1:5800,which(FPdata$groups == "TMSecfp4")),], paired = TRUE)

wilcox.test(APE ~ groups, 
            data = FPdata[c(1:5800,which(FPdata$groups == "TMSlayered")),], paired = TRUE)

wilcox.test(APE ~ groups, 
            data = FPdata[c(1:5800,which(FPdata$groups == "TMSmaccs")),], paired = TRUE)

wilcox.test(APE ~ groups, 
            data = FPdata[c(1:5800,which(FPdata$groups == "TMSobfp2")),], paired = TRUE)

wilcox.test(APE ~ groups, 
            data = FPdata[c(1:5800,which(FPdata$groups == "TMSpb")),], paired = TRUE)

#RI of derivatized metabolites with 1TMS vs 2TMS


#Import Golm Metabolome Database containing metabolites with TMS derivetives metabolites
load(gmdRI.rda)
##Extracting and filtering database information
db <- data.frame(matrix(vector(), nrow = length(gmdRi), ncol = 2, dimnames = list(c(), c("Name", "RI.VAR5.ALK"))))
for (i in 1:length(gmdRi@database)) {
  db[i,1] <- gmdRi@database[[i]]$Name
  db[i,2] <- as.numeric(gmdRi@database[[i]]$RI.VAR5.ALK)
}
db$TMS <- as.numeric(unlist(lapply(db$Name, function(x){
  s <- unlist(strsplit(x, split = "TMS"))[1]
  s <- substring(s, nchar(s))
  return(s)
})))

db$CompName <- unlist(lapply(db$Name, function(x){
  s <- unlist(strsplit(x, split = "TMS"))[1]
  #s <- substring(s, nchar(s)-3)
  s <- substr(s, start = 1, stop = nchar(s)-3)
  return(s)
}))
for (i in which(is.na(db$TMS))) {
  db[i,"CompName"] <- db[i,"Name"]
}
db$id <- as.numeric(row.names(db))
##Select only the metabolites we calculated different fingerpring classes 
##(Dragon fingerprints were took as reference but metabolites are the same regarding all the sets)
dbFilt <- db[df[[1]][,1],]
dbFilt[is.na(dbFilt)] <- 0

barplot(table(dbFilt$TMS), col = "#91c3e9",
        las = 1,
        axes = TRUE, border = NA, xlab = "No. TMS groups", ylab = "No. metabolites")


TMS1_2 <- merge(dbFilt[which(dbFilt$TMS == 1),], dbFilt[which(dbFilt$TMS == 2),], by = "CompName")
TMS2_3 <- merge(dbFilt[which(dbFilt$TMS == 2),], dbFilt[which(dbFilt$TMS == 3),], by = "CompName")
TMS1_2Filt <- TMS1_2[-which(TMS1_2$CompName == "Phenethylamine"),] #Phenethylamine is an outlier so we discard it to calculate the linear model

TMS1_2$Pred2TMS <- predict(lm(RI.VAR5.ALK.y~RI.VAR5.ALK.x, data = TMS1_2Filt), 
                           newdata = data.frame("RI.VAR5.ALK.x" = TMS1_2$RI.VAR5.ALK.x))
TMS2_3$Pred3TMS <- predict(lm(RI.VAR5.ALK.y~RI.VAR5.ALK.x, data = TMS2_3), 
                           newdata = data.frame("RI.VAR5.ALK.x" = TMS2_3$RI.VAR5.ALK.x))

TMS1_2$PredError <- (abs(TMS1_2$RI.VAR5.ALK.y - TMS1_2$Pred2TMS)/TMS1_2$RI.VAR5.ALK.y)*100
TMS2_3$PredError <- (abs(TMS2_3$RI.VAR5.ALK.y - TMS2_3$Pred3TMS)/TMS2_3$RI.VAR5.ALK.y)*100

#RI of derivatized metabolites with 1TMS vs 2TMS
plot(x = TMS1_2Filt$RI.VAR5.ALK.x, y = TMS1_2Filt$RI.VAR5.ALK.y,
     xlim=c(1000,3100), ylim=c(1000,3100), pch=20, cex=2,
     xlab = "RI (1 TMS)", ylab = "RI (2 TMS)",
     cex.axis = 1.2, cex.lab = 1.5
)
abline(lm(RI.VAR5.ALK.y~RI.VAR5.ALK.x, data = TMS1_2Filt), col = "red", 
       lty = 2, lwd = 2)
abline(coef = c(0,1), lty = 4, lwd = 2)
#box plots showing the corresponding relative prediction errors as determined by linear regression
boxplot(TMS1_2$PredError, ylim = c(0,8),
        horizontal = TRUE, col = "red", outline = FALSE, 
        frame.plot = FALSE, sub = "RI error (%)",
        cex.sub = 1.5, cex.axis = 1.5, lty = 1, boxlty = 0)

#RI of derivatized metabolites with 2TMS vs 3TMS
plot(x = TMS2_3$RI.VAR5.ALK.x, y = TMS2_3$RI.VAR5.ALK.y,
     xlim=c(1000,3100), ylim=c(1000,3100), pch=20, cex=3,
     xlab = "RI (2 TMS)", ylab = "RI (3 TMS)",
     cex.axis = 1.2, cex.lab = 1.5
)
abline(lm(RI.VAR5.ALK.y~RI.VAR5.ALK.x, data = TMS2_3), col = "red", 
       lty = 2, lwd = 2)
abline(coef = c(0,1), lty = 4, lwd = 2)
#box plots showing the corresponding relative prediction errors as determined by linear regression
boxplot(TMS2_3$PredError, ylim = c(0,10),
        horizontal = TRUE, col = "red", outline = FALSE, 
        frame.plot = FALSE, sub = "RI error (%)",
        cex.sub = 1.5, cex.axis = 1.5, lty = 1, boxlty = 0)

#Training set structural similarity influence on prediction performance
getTanMat <- function(fingMat, includeDiag=FALSE, includeDiagAndUpperTri=FALSE){
  Z <- as.matrix(fingMat) %*% t(as.matrix(fingMat))
  O_AB <- matrix(rep(diag(Z), ncol(Z)), ncol=ncol(Z))
  AB <- O_AB + t(O_AB)
  tanimotoMatrix <- Z/(AB-Z)
  if(!includeDiag) diag(tanimotoMatrix) <- 0
  if(!includeDiagAndUpperTri) tanimotoMatrix[lower.tri(tanimotoMatrix, diag=TRUE)] <- 0
  tanimotoMatrix
}

TanMatTMS <- as.data.frame(getTanMat(df[[1]][,-c(1,2)], includeDiag = FALSE, includeDiagAndUpperTri = TRUE))
dim(TanMatTMS)

rownames(TanMatTMS) <- df[[1]][,1]
colnames(TanMatTMS) <- df[[1]][,1]

idInfo <- list()
for (i in itN) {
  idInfo[[i]] <- data.frame(matrix(c(1:290, random_ids[[i]]), ncol = 2,
                                   dimnames = list(NULL, c("nrow", "random_ids"))))
}

idInfo2 <- idInfo
for (i in itN) {
  idInfo2[[i]] <- merge(idInfo2[[i]], resSVMinear[[1]][[i]], by = 0, all = TRUE)
}

idInfo.df <- NULL
for (i in itN) {
  idInfo.df <- rbind(idInfo.df, idInfo2[[i]])
}

simm1_0.9 <- which(idInfo.df$Simm <= 1 & idInfo.df$Simm >= 0.9)
simm0.9_0.8 <- which(idInfo.df$Simm < 0.9 & idInfo.df$Simm >= 0.8)
simm0.8_0.7 <- which(idInfo.df$Simm < 0.8 & idInfo.df$Simm >= 0.7)
simm0.7_0.6 <- which(idInfo.df$Simm< 0.7 & idInfo.df$Simm >= 0.6)
simm0.6_0 <- which(idInfo.df$Simm < 0.6 & idInfo.df$Simm >= 0)

##
minMatch <- length(simm0.6_0)

randomSimm1_0.9 <- sample(simm1_0.9, size = minMatch, replace = FALSE)
randomSimm0.9_0.8 <- sample(simm0.9_0.8, size = minMatch, replace = FALSE)
randomSimm0.8_0.7 <- sample(simm0.8_0.7, size = minMatch, replace = FALSE)
randomSimm0.7_0.6 <- sample(simm0.7_0.6, size = minMatch, replace = FALSE)
randomSimm0.6_0 <- sample(simm0.6_0, size = minMatch, replace = FALSE)

SimmSVM1_0.9 <- idInfo.df[randomSimm1_0.9,]
SimmSVM0.9_0.8 <- idInfo.df[randomSimm0.9_0.8,]
SimmSVM0.8_0.7 <- idInfo.df[randomSimm0.8_0.7,]
SimmSVM0.7_0.6 <- idInfo.df[randomSimm0.7_0.6,]
SimmSVM0.6_0 <- idInfo.df[randomSimm0.6_0,]

library(RColorBrewer)

hist(idInfo.df$Simm,
     col = rgb(red = 0.03,green = 0.32,blue = 0.61, alpha = 0.6, maxColorValue = 1), breaks = 50, border = NA,
     xlim = c(0,1), ylim = c(0,1500), cex.axis = 1.3,
     xlab = "Tanimoto Value", cex.lab = 1.5, lwd = 3,
     main = "")
hist(c(SimmSVM1_0.9$Simm, SimmSVM0.9_0.8$Simm, SimmSVM0.8_0.7$Simm, 
       SimmSVM0.7_0.6$Simm, SimmSVM0.6_0$Simm),
     col = rgb(red = 0.62,green = 0.79,blue = 0.88, alpha = 0.6, maxColorValue = 1), 
     breaks = 50, border = NA,
     main = "", add = TRUE)
legend("topleft", inset = .02, 
       legend=c("N=5800", "N=423"),
       cex = 1.5,
       fill=brewer.pal(n = 6, name = "Blues")[c(6,4)],
       border = NA, box.lty = 0
)

Simm.df <- SimmSVM1_0.9
Simm.df <- rbind(Simm.df, SimmSVM0.9_0.8)
Simm.df <- rbind(Simm.df, SimmSVM0.8_0.7)
Simm.df <- rbind(Simm.df, SimmSVM0.7_0.6)
Simm.df <- rbind(Simm.df, SimmSVM0.6_0)
Simm.df$Groups <- as.factor(c(rep("1-0.9", nrow(SimmSVM1_0.9)),
                              rep("0.9-0.8", nrow(SimmSVM0.9_0.8)),
                              rep("0.8-0.7", nrow(SimmSVM0.8_0.7)),
                              rep("0.7-0.6", nrow(SimmSVM0.7_0.6)),
                              rep("0.6-0", nrow(SimmSVM0.6_0))
))

boxplot(SimmSVM1_0.9[,8]*100,
        SimmSVM0.9_0.8[,8]*100,
        SimmSVM0.8_0.7[,8]*100,
        SimmSVM0.7_0.6[,8]*100,
        SimmSVM0.6_0[,8]*100,
        outline = FALSE,
        col = c("#4D004B", "#9569c6", "#2171B5", "#08519C", "#08306B"),
        xlab = "Simmilarity range", ylab = "Relative Error( (%)",
        axes = FALSE, cex.lab = 1.5,
        boxlty = 0, whisklty = 3, whisklwd = 3,
        staplelwd = 2)
axis(side = 2, at = seq(0,20, 2), labels = seq(0,20,2), 
     cex.axis = 1.2, las = 1, lwd = 3)
text(1:5,
     par("usr")[3], 
     srt = 60, adj = 0.9, xpd = TRUE, 
     #adj = c(NA,2), 
     xpd = TRUE,
     labels = c("1-0.9",
                "0.9-0.8",
                "0.8-0.7",
                "0.7-0.6",
                "0.6-0"),
     cex = 1.5)

pairwise.wilcox.test(x = Simm.df$APE, g = Simm.df$Groups, p.adjust.method = "holm")


#Empirical cumulative distribution function
d_fun1_0.9 <- ecdf(SimmSVM1_0.9$APE*100)
d_fun0.9_0.8   <- ecdf(SimmSVM0.9_0.8$APE*100)
d_fun0.8_0.7 <- ecdf(SimmSVM0.8_0.7$APE*100)
d_fun0.7_0.6 <- ecdf(SimmSVM0.7_0.6$APE*100)
d_fun0.6_0 <- ecdf(SimmSVM0.6_0$APE*100)

d_fun1_0.9(1)
d_fun0.9_0.8(1)
d_fun0.8_0.7(1)
d_fun0.7_0.6(1)
d_fun0.6_0(1)

plot(d_fun0.6_0, do.points =FALSE, col.01line = NULL,
     main = "",
     xlab = "Relative Error (%)", ylab = "Cumulative probability",
     axes = FALSE, xlim = c(0,30),
     col = "#08306B", lwd = 1.5)
axis(side = 1, at = seq(0,60,10), labels = seq(0,60,10),
     cex.axis = 1.5, las = 1, lwd = 3)
axis(side = 2, at = seq(0,1, 0.1), labels = seq(0,1, 0.1), 
     cex.axis = 1.5, las = 1, lwd = 3)
plot(d_fun0.7_0.6, do.points = FALSE, col.01line = NULL,
     add = TRUE, 
     col = "#08519C", lwd = 1.5)
plot(d_fun0.8_0.7, do.points = FALSE, col.01line = NULL,
     add = TRUE, col = "#2171B5", lwd = 1.5)
plot(d_fun0.9_0.8, do.points = FALSE, col.01line = NULL,
     add = TRUE, col = "#9569c6", lwd = 1.5)
plot(d_fun1_0.9, do.points = FALSE, col.01line = NULL,
     add = TRUE, col = "#4D004B", lwd = 1.5)
legend("bottomright",
       legend = c("1-0.9", "0.9-0.8", "0.8-0.7", "0.7-0.6", "0.6-0"), 
       col = c("#4D004B", "#9569c6", "#2171B5", "#08519C", "#08306B"), 
       pch = 19, 
       bty = "n", 
       pt.cex = 2, 
       cex = 1.2,
       text.col = "black", 
       horiz = F , 
       inset = c(0.1, 0.1))

























