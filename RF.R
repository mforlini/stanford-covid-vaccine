##
## Select hyperparameters and train Random Forest models for RNA reactivity predictions
##

# Helper packages
library(dplyr)    # for data wrangling
library(ggplot2)  # for awesome graphics

# Modeling packages
library(ranger)   # a c++ implementation of random forest 


############################################################## 

# Read in data
all.data = read.csv("Processed_Data\\10rna.csv")

dim(all.data)
#Dimensions should be: [1] 140488     44



##############################################################
# Use reactivity values to select hyperparameters. Ideally, optimal
#   hyper parameters would be selected for each of the 5 columns, or
#   at least the 3 columns that will be scored, however, computing 
#   time was limited and I experienced a lot of crashing so I only 
#   used reactivity.


# Get one RNA id from each cluster that will be used for reactivity
#   NB: We could alternatively use a weighting system, however this 
#   method reduced the size of the training data set and therefore
#   improved run time. For default hyperparameter values, it did not
#   affect OOB prediction error much.
keep.ids <- all.data %>%
  as_tibble %>%
  filter(reactivity_weights>0) %>%
  distinct(cluster, .keep_all = TRUE)

# Select only RNA ids from keep.ids (one from each cluster with positive
#   weight for the reactivity model), and select columns for features
#   and target (remove other targets, weights, and id/index/cluster 
#   columns)
data.reactivity = all.data %>%
  as_tibble %>%
  filter(id %in% keep.ids$id) %>%
  select(-X,-id,-index,-cluster,-deg_50C,-deg_Mg_pH10,-deg_Mg_50C,-deg_pH10,-reactivity_weights,-deg_50C_weights,-deg_Mg_50C_weights,-deg_Mg_pH10_weights,-deg_pH10_weights) 

n_features = length(setdiff(names(data.reactivity),"reactivity"))
#colnames(data.reactivity)

# Train a random forest with default hyperparameters in order
#   to get baseline OOB RMSE for our hyperparameter search
reactivity.rf1 <- ranger( 
  reactivity ~ .,
  data = data.reactivity,
  respect.unordered.factors="partition",
  num.trees = n_features*10,
  mtry = floor(n_features/3),
  save.memory = TRUE,
  importance = "impurity",
  seed=123
)

# Store default RMSE as baseline for hyperparameter search
(default_rmse <- sqrt(reactivity.rf1$prediction.error))
#[1] 0.219561


# Set up grid search for hyperparameters
hyper_grid <- expand.grid(
  mtry = floor(n_features * c(.15, .25, .333, .4, .6)),
  min.node.size = c(3, 5, 10), 
  replace = c(TRUE, FALSE),                               
  sample.fraction = c(.5, .63, .8),                       
  rmse = NA                                               
)  

# Execute full cartesian grid search
for(i in seq_len(nrow(hyper_grid))) {
  if(is.na(hyper_grid$rmse[i])){
    # fit model for ith hyperparameter combination
    print(paste0("Fitting model ",i, " with parameters.... "))
    print(paste("\t mtry =", hyper_grid$mtry[i], ", min.node.size =", hyper_grid$min.node.size[i], ", replace =",hyper_grid$replace[i], ", sample.fraction =", hyper_grid$sample.fraction[i]))
    fit <- ranger(
      formula         = reactivity ~ ., 
      data            = data.reactivity, 
      num.trees       = n_features * 10,
      mtry            = hyper_grid$mtry[i],
      min.node.size   = hyper_grid$min.node.size[i],
      replace         = hyper_grid$replace[i],
      sample.fraction = hyper_grid$sample.fraction[i],
      verbose         = FALSE,
      seed            = 123,
      save.memory     = TRUE,
      respect.unordered.factors = 'order',
    )
    # export OOB error 
    hyper_grid$rmse[i] <- sqrt(fit$prediction.error)
  } 
}

# Assess top 10 models, select optimal hyperparameters
hyper_grid %>%
  arrange(rmse) %>%
  mutate(perc_gain = (default_rmse - rmse) / default_rmse * 100) %>%
  head(10)



# Check for optimal number of trees

# Set up grid for recording rmse
ntrees_grid = expand.grid(
  ntrees = c(1:15)*n_features,
  rmse = NA                                               
)  

for(i in seq_len(nrow(ntrees_grid))){
  fit <- ranger(
    reactivity ~ .,
    data = data.reactivity,
    num.trees = ntrees_grid$ntrees[i],
    mtry = floor(n_features/3),
    seed=123
  )
  
  ntrees_grid$rmse[i] = fit$prediction.error
}

# Plot number of trees vs OOB error
plot <- ggplot(data = ntrees_grid, mapping=aes(ntrees, rmse))  +
  geom_line(mapping=aes(y=rmse)) +
  geom_point(mapping=aes(y=rmse)) +
  ylim(0,NA) 



# ------------------------
# TRAIN FINAL MODEL
# ------------------------

hyper_selection = c("ntrees" = 10*n_features, "mtry" = 11, "min.node.size" = 3, "replace" = FALSE, "sample.fraction"=.63)

target_names = c("reactivity","deg_50C","deg_Mg_pH10","deg_Mg_50C","deg_pH10")
weight_names = (paste0(target_names,"_weights"))

# Set up data and weight variables for each column to be predicted
for(i in seq_len(length(target_names))){
  print(i)
# Get one id from each cluster for reactivity
# keep.ids <- all.data %>%
#   as_tibble %>%
#   filter((weight_names[i])>0) #%>%
#   #distinct(cluster, .keep_all = TRUE)

  temp <- all.data %>%
    as_tibble %>%
    filter(!!as.name(weight_names[i])>0) %>%
    select(-X,-id,-index,-cluster,-all_of(target_names[-i]),-all_of(weight_names[-i])) 

assign(paste0("weights.",target_names[i]),
       temp %>%
         select(!!as.symbol(weight_names[i]))
       )

assign(paste0("data.",target_names[i]),
       temp %>%
         select(-!!as.symbol(weight_names[i]))
)

}


# Fit 5 models
fit.deg_pH10 <- ranger(
  formula         = deg_pH10 ~ ., 
  data            = data.deg_pH10, #get(paste0("data.",target_names[i])), 
  num.trees       = hyper_selection["ntrees"],
  mtry            = hyper_selection["mtry"],
  min.node.size   = hyper_selection["min.node.size"],
  replace         = hyper_selection["replace"],
  sample.fraction = hyper_selection["sample.fraction"],
  case.weights    = get(paste0("weights.",target_names[i])),
  verbose         = TRUE,
  seed            = 123,
  save.memory     = TRUE,
  respect.unordered.factors = 'partition'
)

saveRDS(fit.deg_pH10,file="RFmodel_pH10.rds")
remove(fit.deg_pH10)
gc()


fit.reactivity <- ranger(
  formula         = reactivity ~ ., 
  data            = data.reactivity, #get(paste0("data.",target_names[i])), 
  num.trees       = hyper_selection["ntrees"],
  mtry            = hyper_selection["mtry"],
  min.node.size   = hyper_selection["min.node.size"],
  replace         = hyper_selection["replace"],
  sample.fraction = hyper_selection["sample.fraction"],
  case.weights    = get(paste0("weights.",target_names[i])),
  verbose         = TRUE,
  seed            = 123,
  save.memory     = TRUE,
  respect.unordered.factors = 'partition'
)

saveRDS(fit.reactivity,file="RFmodel_reactivity.rds")
remove(fit.deg_reactivity)
gc()


fit.deg_50C <- ranger(
  formula         = deg_50C ~ ., 
  data            = data.deg_50C, #get(paste0("data.",target_names[i])), 
  num.trees       = hyper_selection["ntrees"],
  mtry            = hyper_selection["mtry"],
  min.node.size   = hyper_selection["min.node.size"],
  replace         = hyper_selection["replace"],
  sample.fraction = hyper_selection["sample.fraction"],
  case.weights    = get(paste0("weights.",target_names[i])),
  verbose         = TRUE,
  seed            = 123,
  save.memory     = TRUE,
  respect.unordered.factors = 'partition'
)

saveRDS(fit.deg_50C,file="RFmodel_50C.rds")
remove(fit.deg_50C)
gc()


fit.deg_Mg_50C <- ranger(
  formula         = deg_Mg_50C ~ ., 
  data            = data.deg_Mg_50C, #get(paste0("data.",target_names[i])), 
  num.trees       = hyper_selection["ntrees"],
  mtry            = hyper_selection["mtry"],
  min.node.size   = hyper_selection["min.node.size"],
  replace         = hyper_selection["replace"],
  sample.fraction = hyper_selection["sample.fraction"],
  case.weights    = get(paste0("weights.",target_names[i])),
  verbose         = TRUE,
  seed            = 123,
  save.memory     = TRUE,
  respect.unordered.factors = 'partition'
)

saveRDS(fit.deg_Mg_50C,file="Mg_50C_model.rds")
remove(fit.deg_Mg_50C)
gc()

fit.deg_Mg_pH10 <- ranger(
  formula         = deg_Mg_pH10 ~ ., 
  data            = data.deg_Mg_pH10, #get(paste0("data.",target_names[i])), 
  num.trees       = hyper_selection["ntrees"],
  mtry            = hyper_selection["mtry"],
  min.node.size   = hyper_selection["min.node.size"],
  replace         = hyper_selection["replace"],
  sample.fraction = hyper_selection["sample.fraction"],
  case.weights    = get(paste0("weights.",target_names[i])),
  verbose         = TRUE,
  seed            = 123,
  save.memory     = TRUE,
  respect.unordered.factors = 'partition'
)

saveRDS(fit.deg_Mg_pH10,file="Mg_pH10_model.rds")
remove(fit.deg_Mg_pH10)
gc()


