Script plots base case
================
S. Pinto
December 2022

## Load required packages

``` r
library( tidyverse )
library( car )
library( ggplot2 )
library( ggsignif )
library( ggpubr )
library( phyloseq )
```

## Load the data

# Investigate resident species

``` r
combined_df <- data.frame()

# Iterate over each element in the list
for (i in 1:length(listofresults)) {
  # Get the current data frame from the list
  current_df <- listofresults[[i]]
  
  # Add an extra column with the list element name
  current_df$ListElementName <- names(listofresults)[i]
  
  # Bind the current data frame to the combined data frame
  combined_df <- rbind(combined_df, current_df)
}

# Print the combined data frame
#print(combined_df)

# Filter the combined data frame based on the conditions
filtered_df <- combined_df %>%
  filter(if_else(!is.na(Week14), Week14 == "Resident",
                 if_else(!is.na(Week10), Week10 == "Resident",
                         if_else(!is.na(Week8), Week8 == "Resident", FALSE))))

# Print the filtered data frame
#print(filtered_df)

# Separate the "patient.characteristics" column into separate columns
separated_df <- filtered_df %>%
  separate(patient.characteristics, into = c("subject_id", "State", "Donor", "Pretreatment", "Age", "Sex"), sep = " - ")

# Print the separated data frame
#print(separated_df)

# Get the abundance values per mOTU and subject_id
physeq_mOTU <- readRDS( "/Volumes/My Passport for Mac/FECBUD studie/FECBUD GitHub/Data/FECBUD_physeq_filt_relativeabundances.rds" ) 

# Remove patients 109 and 117 because they have too many missing timepoints
physeq_mOTU.new = subset_samples( physeq_mOTU, subject_id != "109" &  subject_id != "117" )

abundance_data <- phyloseq::otu_table(physeq_mOTU.new) %>%
  t() %>%
    as.data.frame() %>% 
  rownames_to_column(var = "id") %>%
  separate(id, into = c("subject_id", "timepoint"), sep = "-") %>% 
  filter(timepoint == "FMT1" | timepoint == "Week8" | timepoint ==  "Week10" | timepoint ==  "Week14")

abundance_data$timepoint[abundance_data$timepoint == 'FMT1'] <- 'Pre-FMT'

#abundance_data %>% colnames()
column_names <- colnames(abundance_data)
# Get the index of the last column
last_column_index <- ncol(abundance_data)
# Get the name of the last column
last_column_name <- column_names[last_column_index]

abundance_data_long <- gather(abundance_data, mOTUs, abundance, ref_mOTU_v3_00077:ext_mOTU_v3_33798, factor_key=TRUE)
abundance_data_wide <- reshape(abundance_data_long, idvar = c("mOTUs","subject_id"), timevar = "timepoint", direction = "wide")

# Join the abundance data to the separated_df
final_df <- left_join(separated_df, abundance_data_wide, by = c("mOTUs","subject_id"))

# Create the abundances.final column based on the conditions
final_df$abundance.final <- ifelse(!is.na(final_df$abundance.Week14),
                                          final_df$abundance.Week14,
                                          ifelse(!is.na(final_df$abundance.Week10),
                                                 final_df$abundance.Week10,
                                                 ifelse(!is.na(final_df$abundance.Week8),
                                                        final_df$abundance.Week8,
                                                        NA)))

# Remove rows with NA values in the abundances.final column
final_df <- final_df[!is.na(final_df$abundance.final), ]

final_df$abundance.difference <- final_df$"abundance.Pre-FMT" - final_df$abundance.final
final_df$abundance.difference.trans <- sign(final_df$abundance.difference) * sqrt(abs(final_df$abundance.difference))

# Print the updated dataframe
#head(final_df)

# Create the histogram plot
histogram_plot <- ggplot(final_df, aes(x = abundance.difference.trans, fill=State)) +
    #geom_density(color = "black", 
     #            size = 1) +
  geom_histogram(bins = 50, color = "black") +
  facet_wrap(~ factor(subject_id, levels=c('101', '104', '106', '107', '108', '113', '115', '120', '125',
                                           '102', '103', '110', '111', '112', '119', '121', '122', '124')), ncol = 3, scales = "free") + 
  scale_fill_manual(values = c(Good = "#00BFC4", None = "#EC7063")) +
  labs(x = "Difference pre- and post FMT relative abundance (square root)", y = "Count") +
    scale_x_continuous(limits = c(min(final_df$abundance.difference.trans), max(final_df$abundance.difference.trans))) +  # Set the x-axis limits

  #scale_x_log10() +
  theme_bw()

# Add a black vertical line at a specific x-value
x_value <- 0  # Adjust this value as needed
histogram_plot <- histogram_plot +
  geom_vline(xintercept = x_value, color = "black", linetype = "dashed") +
    theme(text = element_text(size=16),
        axis.text.x = element_text(angle=90, hjust=1))

# Print the histogram plots
print(histogram_plot)
```

    ## Warning: Removed 36 rows containing missing values (`geom_bar()`).

![](8---Resident-species_files/figure-gfm/resident%20species-1.png)<!-- -->

``` r
# Calculate t-test statistic and p-value
#t_test_result <- t.test(abundance.difference ~ State, data = final_df)

# Check assumptions
levene_test <- leveneTest(abundance.difference ~ State, data = final_df)
```

    ## Warning in leveneTest.default(y = y, group = group, ...): group coerced to
    ## factor.

``` r
shapiro_test_good <- shapiro.test(final_df$abundance.difference[final_df$State == "Good"])
shapiro_test_none <- shapiro.test(final_df$abundance.difference[final_df$State == "None"])

# Perform Wilcoxon rank-sum test
wilcox_test_result <- wilcox.test(abundance.difference ~ State, data = final_df)

# Create the histogram plot
box_plot <- ggplot(final_df, aes(y = abundance.difference.trans)) +
  geom_boxplot(aes(fill = State)) +
  scale_fill_manual(values = c(Good = "#00BFC4", None = "#EC7063")) +
  labs(y = "Difference in relative abundance (square root)", x = "State") +
  theme_bw() +
  #scale_y_log10() +
  
  # Add t-test results to the plot
  geom_text(data = data.frame(y = max(final_df$abundance.difference.trans) * 1.1),
            aes(x = 0, y = c(y-0.1),
                label = ifelse(wilcox_test_result$p.value < 0.001, "***",
                                ifelse(wilcox_test_result$p.value < 0.01, "**",
                                       ifelse(wilcox_test_result$p.value < 0.05, "*", ""))),
                vjust = -0.5), size = 6) +
  geom_text(data = data.frame(y = max(final_df$abundance.difference.trans) * 1.15),
            aes(x = 0, y = c(y-0.1),
                label = paste("p =", format.pval(wilcox_test_result$p.value, digits = 2)),
                vjust = -0.5), size = 4) +
    
  theme(text = element_text(size=16))

# Print the histogram plot
print(box_plot)
```

![](8---Resident-species_files/figure-gfm/resident%20species-2.png)<!-- -->

``` r
# Load necessary packages
leveneTest(abundance.difference ~ State, data = final_df)
```

    ## Warning in leveneTest.default(y = y, group = group, ...): group coerced to
    ## factor.

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##         Df F value Pr(>F)
    ## group    1  0.0029 0.9572
    ##       1542

``` r
# Check normality assumption (optional)
shapiro_test_good <- shapiro.test(final_df$abundance.difference[final_df$State == "Good"])
shapiro_test_none <- shapiro.test(final_df$abundance.difference[final_df$State == "None"])

# Assuming your data frame is called 'final_df'
# Make sure 'State' and 'Category' columns are defined as factors
final_df$State <- as.factor(final_df$State)

t.test(abundance.difference ~ State, data = final_df)
```

    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  abundance.difference by State
    ## t = 0.19682, df = 1221.4, p-value = 0.844
    ## alternative hypothesis: true difference in means between group Good and group None is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.001679102  0.002053570
    ## sample estimates:
    ## mean in group Good mean in group None 
    ##       0.0008310957       0.0006438617

``` r
# Perform Wilcoxon rank-sum test
wilcox_test_result <- wilcox.test(abundance.difference ~ State, data = final_df)
wilcox_test_result
```

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  abundance.difference by State
    ## W = 282582, p-value = 0.1161
    ## alternative hypothesis: true location shift is not equal to 0

``` r
combined_df <- data.frame()

# Iterate over each element in the list
for (i in 1:length(listofresults)) {
  # Get the current data frame from the list
  current_df <- listofresults[[i]]
  
  # Add an extra column with the list element name
  current_df$ListElementName <- names(listofresults)[i]
  
  # Bind the current data frame to the combined data frame
  combined_df <- rbind(combined_df, current_df)
}

# Print the combined data frame
# print(combined_df)

# Filter the combined data frame based on the conditions
filtered_df <- combined_df %>%
  filter(if_else(!is.na(Week14), Week14 %in% c("Resident", "Host Transient", "Species loss"),
                 if_else(!is.na(Week10), Week10 %in% c("Resident", "Host Transient", "Species loss"),
                         if_else(!is.na(Week8), Week8 %in% c("Resident", "Host Transient", "Species loss"), FALSE))))

# Print the filtered data frame
# print(filtered_df)

# Separate the "patient.characteristics" column into separate columns
separated_df <- filtered_df %>%
  separate(patient.characteristics, into = c("subject_id", "State", "Donor", "Pretreatment", "Age", "Sex"), sep = " - ")

# Print the separated data frame
# print(separated_df)

# Get the abundance values per mOTU and subject_id
abundance_data <- phyloseq::otu_table(physeq_mOTU.new) %>%
  t() %>%
    as.data.frame() %>% 
  rownames_to_column(var = "id") %>%
  separate(id, into = c("subject_id", "timepoint"), sep = "-") %>% 
  filter(timepoint == "FMT1" )

abundance_data$timepoint[abundance_data$timepoint == 'FMT1'] <- 'PreFMT'

#abundance_data %>% colnames()
column_names <- colnames(abundance_data)
# Get the index of the last column
last_column_index <- ncol(abundance_data)
# Get the name of the last column
last_column_name <- column_names[last_column_index]

abundance_data_long <- gather(abundance_data, mOTUs, abundance, ref_mOTU_v3_00077:ext_mOTU_v3_33798, factor_key=TRUE)
abundance_data_wide <- reshape(abundance_data_long, idvar = c("mOTUs","subject_id"), timevar = "timepoint", direction = "wide")

# Join the abundance data to the separated_df
final_df <- left_join(separated_df, abundance_data_wide, by = c("mOTUs","subject_id"))

# Create the abundances.final column based on the conditions
final_df$abundance.State <- ifelse(!is.na(final_df$Week14),
                                          final_df$Week14,
                                          ifelse(!is.na(final_df$Week10),
                                                 final_df$Week10,
                                                 ifelse(!is.na(final_df$Week8),
                                                        final_df$Week8,
                                                        NA)))

# Remove rows with NA values in the abundances.final column
final_df <- final_df[!is.na(final_df$abundance.State), ]
final_df$abundance.PreFMT.trans <- sign(final_df$abundance.PreFMT) * sqrt(abs(final_df$abundance.PreFMT))

# Print the updated dataframe
#head(final_df)

# remove patient 118, because no resident species
final_df <- final_df[!(final_df$subject_id=="118"),]

final_df$abundance.State <- factor(final_df$abundance.State, levels = c("Resident", "Host Transient", "Species loss"))

# Define colors for the levels
green_tints <- c("#00BFC4", "#00B8E5", "#00A5FF") %>% rev()
red_tints <- c("#F3B3A0", "#E1876E", "#EC7063") %>% rev()

# Map colors to each combination of abundance.State and State
final_df$color <- ifelse(final_df$State == "Good",
                         ifelse(final_df$abundance.State == "Resident", green_tints[1],
                                ifelse(final_df$abundance.State == "Host Transient", green_tints[2],
                                       ifelse(final_df$abundance.State == "Species loss", green_tints[3], NA))),
                         ifelse(final_df$abundance.State == "Resident", red_tints[1],
                                ifelse(final_df$abundance.State == "Host Transient", red_tints[2],
                                       ifelse(final_df$abundance.State == "Species loss", red_tints[3], NA))))

# Reorder the factor levels for abundance.State
final_df$abundance.State <- factor(final_df$abundance.State, levels = c("Resident", "Host Transient", "Species loss"))

# Create the histogram plot
histogram_plot2 <- ggplot(final_df, aes(x = abundance.PreFMT.trans, fill = color)) +
  geom_histogram(bins = 30, color = "black") +
  scale_fill_manual(values = c(green_tints, red_tints), 
                    labels = c("Resident - Good", "Host Transient - Good", "Species loss - Good",
                               "Resident - None", "Host Transient - None", "Species loss - None")) +
  facet_wrap(~State, scales = "free") +
  labs(x = "Relative abundance pre-FMT (square root)", y = "Count") +
  theme_bw() +
  #scale_x_log10() +
  facet_wrap(~factor(subject_id, levels = c('101', '104', '106', '107', '108', '113', '115', '120', '125',
                                           '102', '103', '105', '110', '111', '112', '119', '121', '122', '124')),
             ncol = 3)  +
    theme(text = element_text(size=16),
        axis.text.x = element_text(angle=90, hjust=1))

# Print the histogram plots
print(histogram_plot2)
```

![](8---Resident-species_files/figure-gfm/host%20species-1.png)<!-- -->

``` r
# Check assumptions
shapiro_test <- shapiro.test(final_df$abundance.PreFMT)
levene_test <- leveneTest(abundance.PreFMT ~ State, data = final_df)
```

    ## Warning in leveneTest.default(y = y, group = group, ...): group coerced to
    ## factor.

``` r
anno_df = compare_means(abundance.PreFMT ~ abundance.State, group.by = "State", data = final_df) %>%
  mutate(y_pos = c(0.6, 0.8, 0.5, 0.6, 0.8, 0.5))

#### Boxplot
# Create the histogram plot
box_plot2 <- ggplot(final_df, aes(y = abundance.PreFMT.trans, x = abundance.State)) +
  geom_boxplot(aes(fill = color)) +
    #scale_x_continuous(trans=weird) +
  facet_wrap(~ State) + 
  #scale_fill_manual(values = c(Good = "#82E0AA", None = "#EC7063")) +
  scale_fill_identity() +  # Use the color mapping as-is without any transformations
  labs(y = "Relative abundance pre-FMT (square root)", 
       x = "Category") +
  theme_bw() +
  #scale_y_log10() +
  geom_signif(data = anno_df, aes(xmin=group1,
                                  xmax=group2,
                                  annotations=p.signif,
                                  y_position = y_pos),
              manual = TRUE)  +
    theme(text = element_text(size=16))
```

    ## Warning in geom_signif(data = anno_df, aes(xmin = group1, xmax = group2, :
    ## Ignoring unknown aesthetics: xmin, xmax, annotations, and y_position

``` r
  #coord_flip() +

# Print the histogram plots
print(box_plot2)
```

![](8---Resident-species_files/figure-gfm/host%20species-2.png)<!-- -->

``` r
# Make sure 'State' and 'Category' columns are defined as factors
final_df$State <- as.factor(final_df$State)
final_df$abundance.State <- as.factor(final_df$abundance.State)

# Perform t-tests for each category compared between states
category_list <- levels(final_df$abundance.State)
results <- list()

for (category in category_list) {
  final_df_category <- final_df[final_df$abundance.State == category, ]
  #t_test_result <- t.test(abundance.PreFMT ~ State, data = final_df_category)
  wilcox_test_result <- wilcox.test(abundance.PreFMT ~ State, data = final_df_category)

  results[[category]] <- wilcox_test_result #t_test_result
}

# Extract and print the results
for (category in category_list) {
  cat(paste("Category:", category, "\n"))
  cat("-------------------------------\n")
  cat("Test statistic:", results[[category]]$statistic, "\n")
  cat("P-value:", results[[category]]$p.value, "\n")
  cat("\n")
}
```

    ## Category: Resident 
    ## -------------------------------
    ## Test statistic: 297125.5 
    ## P-value: 0.9256246 
    ## 
    ## Category: Host Transient 
    ## -------------------------------
    ## Test statistic: 71287 
    ## P-value: 0.5680374 
    ## 
    ## Category: Species loss 
    ## -------------------------------
    ## Test statistic: 140243 
    ## P-value: 7.630309e-07

``` r
# perform t-test pairwise between categories within a state
# Create a function
perform_pairwise_t_tests <- function(data, state) {
  category_list <- levels(data$abundance.State)
  results <- list()
  
  for (category1 in category_list) {
    for (category2 in category_list) {
      if (category1 != category2) {
        data_category1 <- data[data$abundance.State == category1 & data$State == state, ]
        data_category2 <- data[data$abundance.State == category2 & data$State == state, ]
        
        t_test_result <- t.test(data_category1$abundance.PreFMT, data_category2$abundance.PreFMT)
        
        result_label <- paste("Category:", category1, "vs", category2)
        results[[result_label]] <- t_test_result
      }
    }
  }
  
  return(results)
}

# Perform pairwise t-tests within each state
state_list <- levels(final_df$State)
all_results <- list()

for (state in state_list) {
  state_results <- perform_pairwise_t_tests(final_df, state)
  all_results[[state]] <- state_results
}

# Print the results
for (state in state_list) {
  cat("State:", state, "\n")
  cat("-------------------------------\n")
  state_results <- all_results[[state]]
  for (result_label in names(state_results)) {
    cat(result_label, "\n")
    cat("Test statistic:", state_results[[result_label]]$statistic, "\n")
    cat("P-value:", state_results[[result_label]]$p.value, "\n")
    cat("\n")
  }
}
```

    ## State: Good 
    ## -------------------------------
    ## Category: Resident vs Host Transient 
    ## Test statistic: 5.849521 
    ## P-value: 6.252394e-09 
    ## 
    ## Category: Resident vs Species loss 
    ## Test statistic: 8.509327 
    ## P-value: 6.784199e-17 
    ## 
    ## Category: Host Transient vs Resident 
    ## Test statistic: -5.849521 
    ## P-value: 6.252394e-09 
    ## 
    ## Category: Host Transient vs Species loss 
    ## Test statistic: 2.193716 
    ## P-value: 0.02864809 
    ## 
    ## Category: Species loss vs Resident 
    ## Test statistic: -8.509327 
    ## P-value: 6.784199e-17 
    ## 
    ## Category: Species loss vs Host Transient 
    ## Test statistic: -2.193716 
    ## P-value: 0.02864809 
    ## 
    ## State: None 
    ## -------------------------------
    ## Category: Resident vs Host Transient 
    ## Test statistic: 5.206749 
    ## P-value: 2.392041e-07 
    ## 
    ## Category: Resident vs Species loss 
    ## Test statistic: 4.545848 
    ## P-value: 6.192827e-06 
    ## 
    ## Category: Host Transient vs Resident 
    ## Test statistic: -5.206749 
    ## P-value: 2.392041e-07 
    ## 
    ## Category: Host Transient vs Species loss 
    ## Test statistic: -1.213345 
    ## P-value: 0.2253138 
    ## 
    ## Category: Species loss vs Resident 
    ## Test statistic: -4.545848 
    ## P-value: 6.192827e-06 
    ## 
    ## Category: Species loss vs Host Transient 
    ## Test statistic: 1.213345 
    ## P-value: 0.2253138
