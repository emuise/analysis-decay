# Load the boot package
library(boot)

# Sample dataset
data <- data.frame(
  x1 = rnorm(100, mean = 20),
  x2 = rnorm(100, mean = 40),
  x3 = rnorm(100)
)

mean_fun <- function(data, indices) {
  col_means <- sapply(data[indices, ], mean)
  return(col_means)
}

# Number of bootstrap iterations
R <- 1000

# Perform bootstrapping
boot_results <- boot(data = data, statistic = mean_fun, R = R)

get_ci <- function(col_index) {
  df <- boot.ci(boot_results, type = "bca", index = col_index)$bca |>
    as.data.frame()
  
  names(df) <- c("conf", "bias", "skewness", "lower", "upper")
  
  df["col_index"] = col_index
  
  df
}

df_list <- lapply(1:ncol(data), get_ci)

do.call(rbind, df_list)
