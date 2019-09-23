# Xiaojing Zhu 
# Random Survival Forest to Analyze High-Dimensional Breast Cancer Data

########## Preparation ##########
# Load the libraries
library(randomForestSRC)
library(ggRandomForests)
library(ggplot2)
library(Hmisc)
library(akima)

# Import the NKI breast canser data
data <- read.csv("~/Desktop/646 Project/NKI_cleaned.csv")
# This high-dimensional data set contains 272 rows and 1570 columns
dim(data)
# Check for missing values - there is no missing value
sum(is.na(data))

# The response time variable is "survival" (continuous) and the response status 
# variable is "eventdeath" (integers of 0 and 1)
summary(data$survival)
summary(data$eventdeath)

# Explore all 10 clinical predictors
class(data$age)
summary(data$age)

class(data$chemo)
table(data$chemo)

class(data$hormonal)
table(data$hormonal) 

class(data$amputation)
table(data$amputation)

class(data$histtype)
table(data$histtype)

class(data$diam) 
summary(data$diam) 

class(data$posnodes) 
summary(data$posnodes) 

class(data$grade) 
table(data$grade) 

class(data$angioinv) 
table(data$angioinv) 

class(data$lymphinfil) 
table(data$lymphinfil) 

# Explore the first four and the last one genomic predictors
summary(data$esr1)
summary(data$G3PDH_570)
summary(data$Contig45645_RC)
summary(data$Contig44916_RC)
summary(data$AF067420)

# Convert integers with binary or nominal interpretations into factors to feed 
# in the rfsrc function
data$chemo <- factor(data$chemo)
data$hormonal <- factor(data$hormonal)
data$amputation <- factor(data$amputation)
data$histtype <- factor(data$histtype)

# Drop irrelevant columns: "Patient", "timerecurrence", "ID", "barcode"
data <- data[-c(1, 2, 6, 16)]

# Training data consists of 70% of the original data
train_size <- floor(0.7 * nrow(data))

# Set a seed to make reproducible sampling
set.seed(2019)
train_ind <- sample(seq_len(nrow(data)), size = train_size)
train <- data[train_ind, ]
test <- data[-train_ind, ]
dim(train)
dim(test)

############ Random Survival Forest Using Clinical Data Only ############
clinical.train <- train[ ,1:12]
clinical.test <- test[ ,1:12]

# Find the best mtry and nodesize
set.seed(2019)
tuneC <- tune(Surv(survival, eventdeath) ~., clinical.train, samptype = "swr")
print(tuneC$rf)

# The following plot.tune function is adapted from Ishwaran & Kogalur (2019).
plot.tune <- function(model , linear = TRUE) {
  x <- model$results[,1]
  y <- model$results[,2]
  z <- model$results[,3]
  so <- interp(x=x, y=y, z=z, linear = linear)
  idx <- which.min(z)
  x0 <- x[idx]
  y0 <- y[idx]
  filled.contour(x = so$x,
                 y = so$y,
                 z = so$z,
                 xlim = range(so$x, finite = TRUE) + c(-2, 2),
                 ylim = range(so$y, finite = TRUE) + c(-2, 2),
                 color.palette =
                   colorRampPalette(c("yellow", "red")),
                 xlab = "nodesize",
                 ylab = "mtry",
                 main = "OOB error for nodesize and mtry",
                 key.title = title(main = "OOB error", cex.main = 1),
                 plot.axes = {axis(1);axis(2);points(x0,y0,pch="x",cex=1,font=2);
                   points(x,y,pch=16,cex=.25)})
}

# Visualize Figure 2
plot.tune(tuneC)

# Find the best ntree
C <- rfsrc(Surv(survival, eventdeath) ~., 
           data = clinical.train, ntree = 2000, 
           nodesize = 45, mtry = 3, sampsize = 120, nsplit = 10, 
           importance = TRUE, tree.err = TRUE, seed = 2019)
print(C)

err.clinical <- gg_error(C)
err.clinical <- na.omit(err.clinical)

# The following data frame gives the optimal number of trees (ntree) and the 
# corresponding error rate
(Cntree.data <- as.data.frame(err.clinical[err.clinical$error == min(err.clinical$error), ]))
names(Cntree.data) <- NULL
rownames(Cntree.data) <- NULL

# Now Cntree is the optimal number of trees
Cntree <- as.numeric(Cntree.data[2])
Cerr <- as.numeric(Cntree.data[1])

# Plot Figure 1
plot(err.clinical) + theme_bw() + geom_point(data = Cntree.data, 
                                             aes(x = Cntree , y = Cerr), 
                                             color = "red") + labs(title = "Finding Optimal ntree for Clinical Model")

# Cmodel refers to Clinical Model, constructed based on optimal ntree, mtry and nodesize
Cmodel <- rfsrc(Surv(survival, eventdeath) ~., 
                data = clinical.train, ntree = Cntree, 
                nodesize = 45, mtry = 3, sampsize = 120, 
                nsplit = 10, importance = TRUE, tree.err = TRUE, seed = 2019)

# Results for the tuned model shown by Table 3
print(Cmodel)

# Not shown in the report. VIMP with CI is reported instead.
plot(gg_vimp(Cmodel)) + theme(legend.position = c(0.8, 0.2)) + labs(fill = "VIMP > 0")

# Plot forest estimated OOB survival probabilities for each training subject over time
# blue for censoring and red for event (Figure 3)
gg_dtaC <- gg_rfsrc(Cmodel)
plot(gg_dtaC, alpha = 0.3)  + scale_color_manual(values = c("blue", "red")) +  theme_bw() + labs(y = "Survival Probability", 
                                                                                                 title = "OOB Prediction from Clinical Model")

# Obtain plots of OOB survival, OOB mortality and OOB Brier score (Figure 3)
plot.survival(Cmodel, plots.one.page = FALSE)

# Obtain plot of OOB cumulative hazard and OOB hazard (Figure 3)
# subset includes all subject in the training data
plot.survival(Cmodel, plots.one.page = FALSE, subset = 1:nrow(train))

# 95% jackknife CI for VIMP (Figure 5)
plot.subsample(subsample(Cmodel, B = 50), cex = 0.7)

# Fit the training model to the testing data set (Table 3)
pred.test.clinical = predict(Cmodel, newdata = clinical.test )
C.clinical <- rcorr.cens(-pred.test.clinical$predicted, 
                         Surv(clinical.test$survival, 
                              clinical.test$eventdeath))["C Index"]

# Testing C Index
C.clinical
# Testing Error = 1 - C Index
1 - C.clinical


########## Random Survival Forest using Both Clinical and Genomic Data ##########
# Find the best mtry and nodesize
set.seed(2019)
tuneB <- tune(Surv(survival, eventdeath) ~., train, samptype = "swr")
print(tuneB$rf)

# Visualize Figure 2
plot.tune(tuneB)

# Find the best ntree 
B <- rfsrc(Surv(survival, eventdeath) ~., data = train, 
           ntree = 2000, nodesize = 30, mtry = 166, sampsize = 120,
           nsplit = 10, importance = TRUE, tree.err = TRUE, seed = 2019)
print(B)

err <- gg_error(B)
err <- na.omit(err)

# The following data frame gives the optimal number of trees (ntree) and the 
# corresponding error rate
(Bntree.data <- as.data.frame(err[err$error == min(err$error), ]))
names(Bntree.data) <- NULL
rownames(Bntree.data) <- NULL

# Now Bntree is the optimal number of trees
Bntree <- as.numeric(Bntree.data[2]) 
Berr <- as.numeric(Bntree.data[1]) 

# Figure 1
plot(err) + theme_bw() + geom_point(data = Bntree.data, 
                                    aes(x = Bntree , y = Berr), 
                                    color = "red") + labs(title = "Finding Optimal ntree for Clinical & Genomic Model")

Bmodel <- rfsrc(Surv(survival, eventdeath) ~., data = train, 
                ntree = Bntree, nodesize = 30, mtry = 166, sampsize = 120, 
                nsplit = 10, importance = TRUE, tree.err = TRUE, seed = 2019)
# Table 2
print(Bmodel)

# Top 15 VIMP
plot(gg_vimp(Bmodel, nvar = 15)) + theme(legend.position = c(0.8, 0.2)) + labs(fill = "VIMP > 0")
imp15 <- sort(Bmodel$importance, decreasing = TRUE)[1:15]
imp15

# The best model selects 594 out of 1564 predictors for us
sum(Bmodel$importance > 0)
length(Bmodel$importance > 0)

# Figure 4, top right panel
gg_dta <- gg_rfsrc(Bmodel)
plot(gg_dta, alpha = 0.3) + scale_color_manual(values = c("blue", "red")) + theme(legend.position = c(0.2, 0.2)) + theme_bw() + labs(y = "Survival Probability", 
                                                                                                                                     title = "OOB Prediction from Clinical & Genomic Model")

# Figure 4, OOB survival, OOB Brier score and OOB mortality
plot.survival(Bmodel, plots.one.page = FALSE)
# Figure 4, OOB cumulative hazard, OOB hazard
plot.survival(Bmodel, plots.one.page = FALSE, subset = 1:nrow(train))

# 95% jackknife CI for VIMP (Figure 5)
plot.subsample(subsample(Bmodel, B = 50), pmax = 15, cex = 0.4)

# Fit the training model to the testing data set (Table 3)
pred.test = predict(Bmodel, newdata = test )

# Testing C Index
C <- rcorr.cens(-pred.test$predicted, 
                Surv(test$survival, test$eventdeath))["C Index"]
# Testing Error = 1 - C Indexx
1 - C

# Partial plot for top 15 VIMP in the best model (Figure 6)
plot.variable(Bmodel, surv.type = "surv", nvar = 15,  
              plots.per.page = 3, time = 5, partial = TRUE)




