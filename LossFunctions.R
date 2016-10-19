## discussion of loss functions

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

# create indicator variables and fix variable names
GeneA_ind_t <- model.matrix(~-1 + as.factor(GeneA), data = TRAIN)
colnames(GeneA_ind_t) <- sub("as.factor(", "", colnames(GeneA_ind_t), fixed = TRUE)
colnames(GeneA_ind_t) <- sub(")", "_", colnames(GeneA_ind_t), fixed = TRUE)
colnames(GeneA_ind_t) <- gsub("*", "", colnames(GeneA_ind_t), fixed = TRUE)
colnames(GeneA_ind_t) <- sub("/", "_", colnames(GeneA_ind_t), fixed = TRUE)

GeneA_ind_v <- model.matrix(~-1 + as.factor(GeneA), data = VALID)
colnames(GeneA_ind_v) <- sub("as.factor(", "", colnames(GeneA_ind_v), fixed = TRUE)
colnames(GeneA_ind_v) <- sub(")", "_", colnames(GeneA_ind_v), fixed = TRUE)
colnames(GeneA_ind_v) <- gsub("*", "", colnames(GeneA_ind_v), fixed = TRUE)
colnames(GeneA_ind_v) <- sub("/", "_", colnames(GeneA_ind_v), fixed = TRUE)

setdiff(colnames(GeneA_ind_t), colnames(GeneA_ind_v))  # 3 levels not present in validation
setdiff(colnames(GeneA_ind_v), colnames(GeneA_ind_t))  # all in validation in training
GeneA_in_common <- intersect(colnames(GeneA_ind_t), colnames(GeneA_ind_v))

GeneB_ind_t <- model.matrix(~-1 + as.factor(GeneB), data = TRAIN)
colnames(GeneB_ind_t) <- sub("as.factor(", "", colnames(GeneB_ind_t), fixed = TRUE)
colnames(GeneB_ind_t) <- sub(")", "_", colnames(GeneB_ind_t), fixed = TRUE)
colnames(GeneB_ind_t) <- gsub("*", "", colnames(GeneB_ind_t), fixed = TRUE)
colnames(GeneB_ind_t) <- sub("/", "_", colnames(GeneB_ind_t), fixed = TRUE)

GeneB_ind_v <- model.matrix(~-1 + as.factor(GeneB), data = VALID)
colnames(GeneB_ind_v) <- sub("as.factor(", "", colnames(GeneB_ind_v), fixed = TRUE)
colnames(GeneB_ind_v) <- sub(")", "_", colnames(GeneB_ind_v), fixed = TRUE)
colnames(GeneB_ind_v) <- gsub("*", "", colnames(GeneB_ind_v), fixed = TRUE)
colnames(GeneB_ind_v) <- sub("/", "_", colnames(GeneB_ind_v), fixed = TRUE)

setdiff(colnames(GeneB_ind_t), colnames(GeneB_ind_v))  # All in training in validation
setdiff(colnames(GeneB_ind_v), colnames(GeneB_ind_t))  # all in validation in training

# put datasets together
train_sub <- as.data.frame(cbind(GeneA_ind_t, GeneB_ind_t, Therapeutic_Dose = TRAIN$Therapeutic_Dose))
valid_sub <- as.data.frame(cbind(GeneA_ind_v, GeneB_ind_v))


## cross-validated loss functions

# define loss function for squared error
sq_err <- function(pred, obs) (obs - pred)^2
# other loss function for large deviations
large_diff <- function(pred, obs) as.numeric(abs(obs - pred) > 16)

V <- 10
N <- nrow(TRAIN)
MSE_cv <- rep(NA, V)  # placeholder for CV MSE estimates
Diff_cv <- rep(NA, V)  # placeholder for CV MSE estimates
# list of row ids by V validation splits
validRows <- split(sample(1:N), rep(1:V, length=N)) 
for(v in seq(V)) {
  tempTRAIN <- train_sub[-validRows[[v]], ]
  tempVALID <- train_sub[validRows[[v]], ]
  fit_cv <- lm(Therapeutic_Dose ~ ., data = tempTRAIN)  # "." in formula means all variables as additive in data
  pred_cv <- predict(fit_cv, newdata = tempVALID)
  ## compute average loss within validation fold
  MSE_cv[v] <- mean(sq_err(pred = pred_cv, obs = tempVALID$Therapeutic_Dose))
  Diff_cv[v] <- mean(large_diff(pred = pred_cv, obs = tempVALID$Therapeutic_Dose))
}

# summarizes across the folds
summary(MSE_cv)
summary(Diff_cv)

