##
## Run Random Forest Predictions for mRNA from models built in RF.R
## 

# Load required packages

library(dplyr)    # for data wrangling
library(ggplot2)  # for awesome graphics

# Modeling packages
library(ranger)   # a c++ implementation of random forest 



# Load up test data

test.data = read.csv("Processed_Data\\10test_rna.csv")


# 
output <- test.data %>%
  as_tibble() %>%
  select(seq_pos)

# Load models from .RDS files, generate predictions, remove models
#   (I had to do this to avoid crashing)
names = c("reactivity","deg_50C","deg_Mg_50C","deg_pH10","deg_Mg_pH10")
i=1
for( i in 1:length(names) ){
  model <- readRDS(paste0(names[i],"_model.rds"))
  predictions <- predict(model, test.data)$predictions
  assign(paste0(names[i],".predictions"), predictions)
  remove( model, predictions)
  gc()
}

output$reactivity <- reactivity.predictions
output$deg_Mg_pH10 <- deg_Mg_pH10.predictions
output$deg_pH10 <- deg_pH10.predictions
output$deg_Mg_50C <- deg_Mg_50C.predictions
output$deg_50C <- deg_50C.predictions

colnames(output)[1] = "id_seqpos"

write.csv(output, file="Predictions\\randomforest.csv")
