# Step 1: Install Required Packages

# Install packages if not already installed
if (!requireNamespace("GEOquery", quietly = TRUE)) BiocManager::install("GEOquery")
library(GEOquery)

if (!requireNamespace("limma", quietly = TRUE)) BiocManager::install("limma")
library(limma)

# For splitting data and evaluating the model
if (!requireNamespace("caret", quietly = TRUE)) install.packages("caret")
library(caret)

if (!requireNamespace("glmnet", quietly = TRUE)) install.packages("glmnet")
library(glmnet)

#pROC: For evaluating model performance using ROC curves.
if (!requireNamespace("pROC", quietly = TRUE)) install.packages("pROC")
library(pROC)

if (!requireNamespace("randomForest", quietly = TRUE)) install.packages("randomForest")
library(randomForest)

library(GEOquery)
# Download the dataset from GEO
gse <- getGEO("GSE15471", GSEMatrix = TRUE)

# Extract expression data (gene expression levels)
expression_data <- exprs(gse[[1]])

# Extract phenotype data (clinical information such as disease status)
phenotype_data <- pData(phenoData(gse[[1]]))

#Preprocessing data
# Log2 transformation to stabilize variance
expression_data <- log2(expression_data + 1)

# Normalize the data to adjust for technical variation between samples
library(limma)
expression_data <- normalizeBetweenArrays(expression_data)

# Design matrix for differential expression analysis
ls()
colnames(phenotype_data)
design <- model.matrix(~ 0 + phenotype_data$characteristics_ch1)

dim(expression_data)
dim(design)

condition <- factor(phenotype_data$status)  # If 'status' indicates the groups

unique(phenotype_data$treatment_protocol_ch1)
unique(phenotype_data$characteristics_ch1)
unique(phenotype_data$characteristics_ch1.1)

condition <- factor(phenotype_data$characteristics_ch1.1)
design <- model.matrix(~ condition)




# Fit a linear model and compute statistics for differential expression
fit <- lmFit(expression_data, design)
fit <- eBayes(fit)

# Select top 100 differentially expressed genes
top_genes <- topTable(fit, number = 100)

library(caret)
set.seed(123)  # For reproducibility

# Split data into training (80%) and testing (20%) sets
trainIndex <- createDataPartition(phenotype_data$characteristics_ch1, p = 0.8, list = FALSE)
trainData <- expression_data[trainIndex, ]
testData <- expression_data[-trainIndex, ]
trainLabels <- phenotype_data$characteristics_ch1[trainIndex]
testLabels <- phenotype_data$characteristics_ch1[-trainIndex]

library(randomForest)
# Train a Random Forest model with 500 trees
rf_model <- randomForest(x = trainData, y = as.factor(trainLabels), ntree = 500)

# Predict on the test data
predictions <- predict(rf_model, testData)

#Step 7: Evaluate the Model
library(pROC)

levels(predictions)  # Check predicted labels
levels(as.factor(testLabels))  # Check true labels

# Double-check that testData only has 78 rows
dim(testData)  # This should return something like [78, n_features]

# Define a random split (80% training, 20% test)
set.seed(123)  # For reproducibility
total_samples <- nrow(expression_data)
test_size <- 78  # Since you have 78 test labels

# Randomly sample indices for the test set
test_indices <- sample(1:total_samples, size = test_size)

# Optionally, create the training set indices (rest of the data)
train_indices <- setdiff(1:total_samples, test_indices)
# Subset the test data using the test indices
testData <- expression_data[test_indices, ]

# Generate predictions on the test data
predictions <- predict(rf_model, testData)

# Check the length of predictions
length(predictions)  # This should be 78, matching the number of test labels
# Check the length of testLabels
print(predictions)

length(testLabels)  # Should be 78


# Example: Extract the rows corresponding to the 78 test samples
testData <- expression_data[test_indices, ]  # test_indices should be the indices of your test set


# Confusion Matrix: Shows accuracy, sensitivity, specificity, etc.
confusionMatrix(predictions, as.factor(testLabels))

testLabels <- factor(testLabels, levels = levels(predictions))

length(predictions)
length(testLabels) 

predictions <- predict(rf_model, testData)

head(predictions)
levels(testLabels)  # It should output something like: "sample: normal", "sample: tumor"

# Check if there are missing values in the target variable (testLabels)
anyNA(testLabels)

# Check if there are missing values in the feature set (trainData)
anyNA(trainData)


# Remove rows with missing values in testLabels
trainData_clean <- trainData[!is.na(testLabels), ]
testLabels_clean <- testLabels[!is.na(testLabels)]

# Impute missing values with the most frequent category (mode)
mode_value <- names(sort(table(testLabels), decreasing = TRUE))[1]
testLabels_imputed <- ifelse(is.na(testLabels), mode_value, testLabels)

# Check the unique values (classes) in testLabels_clean or testLabels_imputed
unique(testLabels_clean)
# or if you imputed
unique(testLabels_imputed)

table(testLabels)

# Correctly assign testLabels, assuming you have the labels in `phenotype_data`
testLabels <- factor(phenotype_data$status)  # Or however the labels are stored

# Assuming the phenotype data has the correct class labels (normal/tumor)
# Correct the testLabels to reflect the proper class labels
testLabels <- factor(phenotype_data$status)

# Check the unique levels again
levels(testLabels)

# If necessary, handle missing values in `testLabels` appropriately
testLabels <- factor(testLabels, levels = c("sample: normal", "sample: tumor"))

# Check the distribution again
table(testLabels)

str(phenotype_data)

unique(phenotype_data$`sample:ch1`)

# Assign the correct labels to testLabels
testLabels <- factor(phenotype_data$`sample:ch1`)

# Check the levels of the testLabels
levels(testLabels)

# Train a random forest model using the clean training data and the test labels
rf_model <- randomForest(testLabels ~ ., data = trainData)

# Check the model output
print(rf_model)

predictions <- predict(rf_model, newdata = testData)  # Predict on the test data
confusionMatrix(predictions, testLabels) #UP TO HERE CODE WORKS 

testLabels <- factor(phenotype_data$status)  # Convert to factor
testLabels <- factor(phenotype_data$characteristics_ch1.1)  # Convert to factor
levels(testLabels)

dim(testData)  # This should show something like [78, n] where n is the number of features
# Generate predictions for the test set
predictions <- predict(rf_model, testData)
length(predictions)


testLabels <- phenotype_data$characteristics_ch1.1
testLabels <- phenotype_data$status
levels(testLabels)  # Check the unique categories in testLabels




# ROC curve: Evaluates the tradeoff between true positive and false positive rates
roc_curve <- roc(testLabels, as.numeric(predictions))
plot(roc_curve)
auc(roc_curve)  # Area under the ROC curve (AUC)

# Save the trained model to disk for future use
saveRDS(rf_model, "pancreatic_cancer_predictive_model.rds")




















