library(cdcfluview)

train_data <- get_flu_data("national", "ilinet", years=1997:2011)
write.csv(train_data, "data/train_data.csv")
test_data <- get_flu_data("national", "ilinet", years=2012:2015)
write.csv(test_data, "data/test_data.csv")
