library(erah)
load("id_eRah.rda")
load("gmdRi.rda")

db <- data.frame(matrix(vector(), nrow = length(gmdRi), ncol = 2, 
                        dimnames = list(c(), c("Name", "RI.VAR5.ALK"))))
for (i in 1:length(gmdRi@database)) {
  db[i,1] <- gmdRi@database[[i]]$Name
  db[i,2] <- as.numeric(gmdRi@database[[i]]$RI.VAR5.ALK)
}

idF_info <- data.frame(Nrow = 1:nrow(idF),
                       idF$Name.1,
                       idF$Name.2,
                       idF$Name.3,
                       idF$Exp.RI)


testNames <- data.frame(Names = unique(sort(c(idF_info$idF.Name.1, 
                                              idF_info$idF.Name.2, 
                                              idF_info$idF.Name.3))))


load("fgSet/fgGmd_dragonTMS.rda")
dataTMS <- as.data.frame(matrix(unlist(fgGmd), 
                                ncol = ncol(fgGmd), nrow = nrow(fgGmd)), 
                         stringsAsFactors = FALSE)
colnames(dataTMS) <- colnames(fgGmd)

db <- db[dataTMS$Name,]
db$Nrowdb <- 1:nrow(db)
db$Nameid <- rownames(db)

dataF <- merge(db, testNames, 
               by.x = "Name", by.y = "Names", 
               all.x = FALSE, all.y = FALSE)

idF_info <- merge(idF_info, db[,1:2], by.x = "idF.Name.1", by.y = "Name", 
                  all.x = TRUE, all.y = FALSE, sort = FALSE, suffixes = 1)
idF_info <- merge(idF_info, db[,1:2], by.x = "idF.Name.2", by.y = "Name", 
                  all.x = TRUE, all.y = FALSE, sort = FALSE, suffixes = 2)
idF_info <- merge(idF_info, db[,1:2], by.x = "idF.Name.3", by.y = "Name", 
                  all.x = TRUE, all.y = FALSE, sort = FALSE, suffixes = 3)
colnames(idF_info) <- c("idF.Name.3", "idF.Name.2","idF.Name.1", "Nrow",
                        "Exp.RI",  "RI.1", "RI.2", "RI.3")
idF_info <- idF_info[order(idF_info$Nrow),c(4,3,2,1,6,7,8,5)]

test_data <- dataTMS[dataF$Nrowdb,-(1:2)]
test_labels <- dataTMS[dataF$Nrowdb,2]*10
train_data <- dataTMS[-dataF$Nrowdb,-(1:2)]
train_labels <- dataTMS[-dataF$Nrowdb,2]*10

library(e1071)

set.seed(123)
linear.SVM <- NULL
linear.SVM <- svm(x = train_data, y = train_labels,
                  kernel = "linear",
                  type = "eps-regression",
                  cost = 100,
                  epsilon = 0.1)
testPred <- predict(linear.SVM, newdata = test_data)
resSVMlinear <- data.frame("EmpiricalRI" = test_labels,
                           "PredictedRI" = testPred,
                           "AE" = c(abs(testPred - test_labels)),
                           "APE" = c(abs((testPred - test_labels) / test_labels)))

resSVMtransf <- resSVMlinear
resSVMtransf[,1:3] <- resSVMtransf[,1:3]/10

resMet <- cbind(dataF, resSVMtransf)

idF_info <- merge(idF_info, resMet[,c(1,6)], by.x = "idF.Name.1", by.y = "Name", 
                  all.x = TRUE, all.y = FALSE)
idF_info <- merge(idF_info, resMet[,c(1,6)], by.x = "idF.Name.2", by.y = "Name", 
                  all.x = TRUE, all.y = FALSE)
idF_info <- merge(idF_info, resMet[,c(1,6)], by.x = "idF.Name.3", by.y = "Name", 
                  all.x = TRUE, all.y = FALSE)
colnames(idF_info) <- c("idF.Name.3", "idF.Name.2", "idF.Name.1", "Nrow", 
                        "RI.1", "RI.2", "RI.3", "Exp.RI", "RI.Pred.1", 
                        "RI.Pred.2", "RI.Pred.3")
idF_info <- idF_info[order(idF_info$Nrow),c(4,3,2,1,5,6,7,9:11,8)]


idF_info <- cbind(idF_info, idF[,c(19:21,23)])

idF_info$RI.Pred.error.1 <- round(abs((idF_info$RI.Pred.1 - idF_info$Exp.RI)/idF_info$Exp.RI)*100, 
                                  digits = 2)
idF_info$RI.Pred.error.2 <- round(abs((idF_info$RI.Pred.2 - idF_info$Exp.RI)/idF_info$Exp.RI)*100, 
                                  digits = 2)
idF_info$RI.Pred.error.3 <- round(abs((idF_info$RI.Pred.3 - idF_info$Exp.RI)/idF_info$Exp.RI)*100, 
                                  digits = 2)
idF_infoFilt <- idF_info[-which(apply(idF_info[,15:17], 1, function(x){length(which(is.na(x) == TRUE))})>= 2),]

idF_infoFilt$Rank.Pred <- unname(apply(idF_infoFilt[,16:18], 1, which.min))

s <- summary(lm(unname(apply(idF_infoFilt[,c("RI.1", "RI.2", "RI.3", "Rank")], 1, 
                             function(x){x[x[4]]}))~idF_infoFilt$Exp.RI))
s.Pred <- summary(lm(unname(apply(idF_infoFilt[,c("RI.Pred.1", "RI.Pred.2", "RI.Pred.3", "Rank.Pred")], 
                                  1, function(x){x[x[4]]}))~idF_infoFilt$Exp.RI))


layout(matrix(c(1,1,2,3,3,4), 3, 2, byrow = FALSE))
plot(x = unname(apply(idF_infoFilt[,c("RI.1", "RI.2", "RI.3", "Rank")], 1, function(x){x[x[4]]})),
     y = idF_infoFilt$Exp.RI,
     pch = 19, cex = 1.5,
     #col = "#67001f",
     xlab = "Reference RI", ylab = "Experimental RI", 
     cex.axis = 1.5, cex.lab = 1.5,
     main = "")
abline(coef = c(0,1), lty = 4, lwd = 1.3)
abline(lm(unname(apply(idF_infoFilt[,c("RI.1", "RI.2", "RI.3", "Rank")], 1, function(x){x[x[4]]}))~idF_infoFilt$Exp.RI), 
       col = "red", lwd = 1.3
)
text(x = 1500, y = 2000, 
     paste("R2 =", round(s$adj.r.squared, 3)),
     cex = 2)
#legend("topleft", inset = .02, legend=c("Match 1", "Match 2", "Match 3"),
#fill=c("#67001f", "#053061", "#d6604d")
#)
boxplot(unname(apply(idF_infoFilt[,c("RI.error.1", "RI.error.2", "RI.error.3", "Rank")], 1, function(x){x[x[4]]})),
        horizontal = TRUE, col = "red", 
        outline = TRUE, frame.plot = FALSE, 
        ylim = c(0,1.5),
        sub = "RI error (%)",
        boxlty = 0, whisklty = 3, whisklwd = 3,
        staplelwd = 2,cex.axis = 1.2)
plot(x = unname(apply(idF_infoFilt[,c("RI.Pred.1", "RI.Pred.2", "RI.Pred.3", "Rank.Pred")], 1, function(x){x[x[4]]})),
     y = idF_infoFilt$Exp.RI,
     pch = 19, cex = 1.5,
     #col = "#67001f",
     xlab = "Predicted RI", ylab = "Experimental RI", 
     cex.axis = 1.5, cex.lab = 1.5,
     main = "")
abline(coef = c(0,1), lty = 4, lwd = 1.3)
abline(lm(unname(apply(idF_infoFilt[,c("RI.Pred.1", "RI.Pred.2", "RI.Pred.3", "Rank.Pred")], 1, function(x){x[x[4]]}))~idF_infoFilt$Exp.RI), 
       col = "red", lwd = 1.3)
text(x = 1500, y = 2000, 
     paste("R2 =", round(s.Pred$adj.r.squared, 3)),
  cex = 2)
#legend("topleft", inset = .02, legend=c("Match 1", "Match 2", "Match 3"),
#fill=c("#67001f", "#053061", "#d6604d")
#)
boxplot(unname(apply(idF_infoFilt[,c("RI.Pred.error.1", "RI.Pred.error.2", "RI.Pred.error.3", "Rank.Pred")], 
                     1, function(x){x[x[4]]})),
        horizontal = TRUE, col = "red", 
        outline = TRUE, frame.plot = FALSE, 
        ylim = c(0,10),
        sub = "RI error (%)",
        boxlty = 0, whisklty = 3, whisklwd = 3,
        staplelwd = 2, cex.axis = 1.2)









