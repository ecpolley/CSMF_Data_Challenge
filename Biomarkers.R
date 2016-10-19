## discussion of biomarkers

## load data
if(file.exists("TRAIN.CSV")) {
  TRAIN <- read.csv("TRAIN.CSV")
} else {
  urlfile <- "https://raw.githubusercontent.com/ecpolley/
    CSMF_Data_Challenge/master/TRAIN.CSV"
  download.file(urlfile, destfile = "TRAIN.CSV")
  TRAIN <- read.csv("TRAIN.CSV")
}

if(file.exists("VALID.CSV")) {
  VALID <- read.csv("VALID.CSV")
} else {
  urlfile <- "https://raw.githubusercontent.com/ecpolley/
    CSMF_Data_Challenge/master/VALID.CSV"
  download.file(urlfile, destfile = "VALID.CSV")
  VALID <- read.csv("VALID.CSV")
}

## load packages
library(ggplot2)
library(reshape2)

## extract biomarker variables and melt (convert wide to long)
TRAIN_bio <- melt(TRAIN[, grep("Biomarker", colnames(TRAIN))])
VALID_bio <- melt(VALID[, grep("Biomarker", colnames(VALID))])

# add dataset label and combine
TRAIN_bio$Dataset <- "TRAIN"
VALID_bio$Dataset <- "VALID"
ALL_bio <- rbind(TRAIN_bio, VALID_bio)

ggplot(ALL_bio, aes(x = value, fill = Dataset)) + geom_histogram(alpha = .4, position = "identity") + facet_wrap(~variable, ncol = 5, scales = "free")
# note default "position" is to stack the histograms. in the facet_wrap, could remove scales = "free" so all figures on same scale

# plot data
library(GGally)
ggpairs(TRAIN[, grep("Biomarker|Therapeutic_Dose", colnames(TRAIN))])


## look at principal components
## Should we use the Biomarkers?
train.pca <- prcomp(TRAIN[, grep("Biomarker", colnames(TRAIN))])
screeplot(train.pca)
biplot(train.pca)