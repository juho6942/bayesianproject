## ============================================================
## 0) Packages
## ============================================================

required_pkgs <- c(
  "tidybayes", "brms", "ggplot2", "metadat",
  "dplyr", "lme4", "posterior",
  "haven", "naniar", "corrplot", "VIM", "R.utils", "tidyr"
)

for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

lapply(required_pkgs, library, character.only = TRUE)

## Optional: cmdstanr backend for faster sampling (if you later fit models)
if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  install.packages(
    "cmdstanr",
    repos = c("https://mc-stan.org/r-packages/", getOption("repos"))
  )
}

## ============================================================
## 1) Data import and basic checks
## ============================================================

# Unzip and read SPSS file (assumes "bacteremia.sav" is gzipped)
gunzip("bacteremia.sav", destname = "bacteremia_unzipped.sav", overwrite = FALSE)
bacteremia <- haven::read_sav("bacteremia_unzipped.sav")

str(bacteremia)
glimpse(bacteremia)
summary(bacteremia)

## ============================================================
## 2) Basic preprocessing
## ============================================================

data <- as.data.frame(bacteremia)

data <- data %>%
  mutate(
    sex = ifelse(tolower(sex) == "male", 1, 0)
  ) %>%
  select(-id)

outcome_var <- "bacteremia"
data[[outcome_var]] <- as.factor(data[[outcome_var]])

n_patients <- nrow(data)
n_features <- ncol(data) - 1

cat("Patients:", n_patients, "\n")
cat("Features (including outcome):", ncol(data), "\n")

table(data[[outcome_var]])
prop.table(table(data[[outcome_var]]))

## ============================================================
## 3) Missingness overview
## ============================================================

missing_pct <- sapply(data, function(x) mean(is.na(x))) * 100
missing_df <- data.frame(
  variable    = names(missing_pct),
  missing_pct = as.numeric(missing_pct)
)

# Sort by most missing
missing_df <- missing_df[order(-missing_df$missing_pct), ]

# View top/bottom variables by missingness
head(missing_df, 15)   # most missing
tail(missing_df, 15)   # least missing

# Barplot of missingness
ggplot(missing_df, aes(x = reorder(variable, missing_pct), y = missing_pct)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Percentage of missing values per variable",
    x = "Variable",
    y = "Missing (%)"
  )

# Heatmap of missing vs observed (rows = patients, cols = variables)
vis_miss(data, sort_miss = TRUE) +
  labs(title = "Missingness pattern (sorted by % missing)")

## ============================================================
## 4) Numeric variables and quick plots
## ============================================================

numeric_vars <- names(data)[sapply(data, is.numeric)]

# Choose up to 20 numeric variables with lowest missingness
miss_num    <- missing_pct[numeric_vars]
top_numeric <- head(names(sort(miss_num)), 20)

top_numeric

# Histograms
for (v in top_numeric) {
  print(
    ggplot(data, aes(x = .data[[v]])) +
      geom_histogram(bins = 30) +
      labs(
        title = paste("Distribution of", v),
        x = v,
        y = "Count"
      )
  )
}

# Boxplots by outcome
for (v in top_numeric) {
  print(
    ggplot(data, aes(
      x = .data[[outcome_var]],
      y = .data[[v]],
      fill = .data[[outcome_var]]
    )) +
      geom_boxplot() +
      labs(
        title = paste(v, "by", outcome_var),
        x = outcome_var,
        y = v
      )
  )
}

## ============================================================
## 5) Correlation structure of numeric variables
## ============================================================

numeric_data <- data %>% dplyr::select(all_of(numeric_vars))

# Correlation matrix using pairwise complete observations
cor_mat <- cor(numeric_data, use = "pairwise.complete.obs")

# Visual correlation matrix
corrplot(
  cor_mat,
  method = "color",
  type   = "upper",
  tl.cex = 0.7,
  tl.col = "black"
)

## Correlation of missingness indicators themselves
miss_ind <- as.data.frame(sapply(data, function(x) as.numeric(is.na(x))))
colnames(miss_ind) <- paste0("miss_", names(data))

miss_cor <- cor(miss_ind)
corrplot(
  miss_cor,
  method = "color",
  type   = "upper",
  tl.cex = 0.5,
  tl.col = "black"
)

## ============================================================
## 6) Missingness vs outcome (MAR-ish diagnostics)
## ============================================================

# Attach outcome to missingness indicators
miss_ind2 <- miss_ind
miss_ind2[[outcome_var]] <- data[[outcome_var]]

# Missing rate by outcome class for each variable
miss_by_outcome <- lapply(names(data), function(var) {
  df <- miss_ind2 %>%
    group_by(.data[[outcome_var]]) %>%
    summarise(
      missing_rate = mean(.data[[paste0("miss_", var)]], na.rm = TRUE),
      .groups      = "drop"
    )
  df$variable <- var
  df
}) %>%
  bind_rows()

# Filter to variables with any missingness
miss_by_outcome_filtered <- miss_by_outcome %>%
  group_by(variable) %>%
  filter(any(missing_rate > 0)) %>%
  ungroup()

# Plot: missing rate by outcome for variables with some missingness
ggplot(
  miss_by_outcome_filtered,
  aes(x = .data[[outcome_var]], y = missing_rate, group = variable)
) +
  geom_col() +
  facet_wrap(~ variable, scales = "free_y") +
  labs(
    title = "Missingness rate by outcome class for each variable",
    x = outcome_var,
    y = "Missing rate"
  )

## ============================================================
## 7) High / low missingness variables
## ============================================================

high_missing <- missing_df %>% dplyr::filter(missing_pct > 20)
low_missing  <- missing_df %>% dplyr::filter(missing_pct < 10)

high_missing   # candidates to drop or treat carefully
low_missing    # good candidates for initial modelling

## ============================================================
## 8) Correlation with outcome (numeric) & missingness difference
## ============================================================

# Numeric outcome (0/1)
data <- data %>%
  mutate(
    bact_num = ifelse(
      .data[[outcome_var]] == levels(.data[[outcome_var]])[2],
      1, 0
    )
  )

# Numeric features
num_vars <- setdiff(names(data)[sapply(data, is.numeric)], "bact_num")

# Correlation with bacteremia
cor_df <- data.frame(
  variable    = num_vars,
  correlation = sapply(
    num_vars,
    function(v) cor(data[[v]], data$bact_num, use = "pairwise.complete.obs")
  )
) %>%
  arrange(correlation)

ggplot(cor_df, aes(x = reorder(variable, correlation), y = correlation)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  theme_minimal()

# Missingness difference between outcome classes
miss_ind_classes <- as.data.frame(sapply(data, function(x) as.numeric(is.na(x))))
colnames(miss_ind_classes) <- names(data)
miss_ind_classes[[outcome_var]] <- data[[outcome_var]]

out_lv <- levels(data[[outcome_var]])

miss_diff <- miss_ind_classes %>%
  pivot_longer(-all_of(outcome_var), names_to = "variable", values_to = "is_miss") %>%
  group_by(variable, .data[[outcome_var]]) %>%
  summarise(rate = mean(is_miss), .groups = "drop") %>%
  pivot_wider(names_from = .data[[outcome_var]], values_from = rate) %>%
  mutate(diff = .data[[out_lv[2]]] - .data[[out_lv[1]]]) %>%
  arrange(diff)

ggplot(miss_diff, aes(x = reorder(variable, diff), y = diff)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  theme_minimal()

