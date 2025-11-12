##**************************************
##  
## Assignment01 -- BINF6210 
## 
## Rameen Jafri 
##
## 2025-10-12
##
## ===================================== ##
##  Step 0: Dependencies                 ##
## ===================================== ##

library(tidyverse)
library(janitor)
library (ggplot2)
library(scales)
library(viridis)
library(ggpubr) #Will need for statistical testing

# ===================================== #
#  Step 1: read + basic data prep       #
# ===================================== #

getwd()

#Read in TSV file and standardize column name.
lyco_full <- read_tsv("../data/lycosidae.tsv") %>% clean_names()

#check data to see what you're working with 
dim(lyco_full)
head(lyco_full)

# Keep only the variables needed downstream
lyco_sub <- lyco_full[, c("inst", "nuc", "nuc_basecount", "marker_code")]

#Turn all -999 values into NAs (missing values)
lyco_sub[lyco_sub == -999] <- NA 

# Drop rows missing nuc or inst (required for all later analyses)
lyco_cleaned <- drop_na(lyco_sub, nuc, inst) 

#check to see how many rows were dropped
nrow(lyco_full) - nrow(lyco_cleaned)

# Quick NA sanity checks (should be all FALSE now)
summary(is.na (lyco_cleaned$inst))
summary(is.na(lyco_cleaned$nuc)) 

#__# LV Added this code to determine the top 10 institutes in base R, makes code 
# rationale flow better
head(sort(table(lyco_cleaned$inst), decreasing = TRUE), n =10)

#--# LV Code above in tidyverse format, pick whichever you prefer!
lyco_cleaned %>%
  group_by(inst) %>%
  count(inst) %>%
  arrange(desc(n))

#rename top 10 depositories in lyco_cleaned (for neat labels when plotting later) 
depo_names <- c(
  "Centre for Biodiversity Genomics" = "CBG",
  "Mined from GenBank, NCBI" = "GenBank",
  "Natural History Museum, London" = "NHM London",
  "Bernardino Rivadavia Natural Sciences Museum" = "BRNSM",
  "University of Guelph" = "UofG",
  "University of Helsinki" = "Helsinki",
  "University of Oulu, Zoological Museum" = "Oulu Zoo") 

#--# LV Creating a new column with the appropriate shortened top 10 depositories
lyco_cleaned <- lyco_cleaned %>%
  mutate(depo_short = recode(inst, !!!depo_names))

# ===================================== #
#  Step 2: sequence QC/feature columns  #
# ===================================== #

# Create helper columns:
# - nuc_scrubbed: remove dashes and spaces before counting characters
# - nuc_scrubbed_count: character count of scrubbed sequence
# - nuc_same: whether our count matches provided basecount (curiosity check)
# - nuc_ranges: bin sequences by length (Short/Medium/Long)

#--# LV Instead of using mutate 4 times, you can just combine all your new columns under 1 mutate function.
lyco_cleaned <- lyco_cleaned %>%
  mutate(nuc_scrubbed = gsub("-","",nuc),
         nuc_scrubbed_count = nchar(nuc_scrubbed),
         nuc_same=nuc_scrubbed_count==nuc_basecount,
         nuc_ranges = cut(nuc_scrubbed_count, breaks = c(0,620,700,Inf), labels = c("Short","Medium", "Long"), include.lowest = TRUE))

# ===================================== #
#  Step 2a: quick exploratory plots     #
# ===================================== #

# Quick view of marker composition (e.g., COI vs ITS2, etc.)
lyco_cleaned %>% 
  count(marker_code, name = "count") %>% 
  arrange(desc(count))

# Base R hist to inspect exact length counts
#--# LV Use ggplot instead of base R to plot histogram
ggplot(lyco_cleaned, aes(x =nuc_scrubbed_count)) +
  geom_histogram(fill = "pink", colour = "black") +
  labs(x = "Sequence length (bp)", y = "Frequency", title = "Distribution of sequence length markers") +
  theme_minimal()
# Shows the majority of sequences are 600–700 bp (expected for COI).
# ITS2 averages ~300 bp for animals; dataset includes many ITS2 entries,
# but lacks a strong spike near 300 bp — proportions will clarify this pattern.

# Distribution of sequence length bins
ggplot(lyco_cleaned, aes(x=nuc_ranges)) + 
  geom_bar () +
  labs(title = "Nuc Seq by Depo", x = "Sequence ranges in BP", y = "Count") + 
  theme_minimal()
# COI typically ~660 bp, so chosen bins reflect ±40 bp:
# Short (0–620), Medium (621–700), and Long (700+).

#Distribution of Sequence Lengths by Marker Type - expected to see clustering at ~660 bp for COI markers and ~300bp for ITS2 markers. Saw clustering at ~660bp for COI as expected, but also saw incredible variation (from under 200bp all the way to 1500bp). Didn't see clustering around 300bp for ITS2.
ggplot(lyco_cleaned, aes(x = marker_code, y = nuc_scrubbed_count)) +
  geom_boxplot(fill = "pink", alpha = 0.7) +
  labs(
    title = "Distribution of Sequence Lengths by Marker Type",
    x = "Marker",
    y = "Sequence Length (bp)"
  ) +
  theme_minimal()
#this plot doesnt provide info about density (for example, CYTB only has 5 counts, but with a lot of variation so it takes up a lot of space) either the dataset needs to be bigger or markers with low counts need to be removed for this to be useful)

#check: number of records per marker
lyco_cleaned %>% count(marker_code, sort = TRUE) %>% head()

# =============================================================== #
#  Step 3: depositors by bin (Short/Medium/Long) with proportions #
# =============================================================== #

#how many in each bin before subsetting
table(lyco_cleaned$nuc_ranges)

# Subset shorter sequences into a new dataframe 
lyco_short <- lyco_cleaned[c(grep("Short", lyco_cleaned$nuc_ranges)), ] 

##--## LV Can use dplyr to make new dataframe same as above, up to you :)
lyco_short <- lyco_cleaned %>%
  filter(nuc_ranges == "Short")

# Count contributions by depositor for short sequences
##--## LV Created a function to count contributions by depositor for short sequences
summary_depo <- function(dataframe, length) { dataframe %>%
    count(inst, name = "n_records") %>% 
    arrange(desc(n_records)) %>% 
    mutate(inst_grouped = if_else(row_number(desc(n_records)) <= 4, inst, "Other")) %>%
    mutate(nuc_ranges = length) %>% 
    group_by(nuc_ranges, inst_grouped) %>%
    summarise(n_records = sum(n_records), .groups = "drop") %>% # %>%
    mutate(prop = n_records / sum(n_records)) %>% 
    arrange(desc(n_records))}

##--## LV Calling on summary_depo function to apply on short sequences
short_by_depo <- summary_depo(lyco_short, "Short")
head(short_by_depo)

# Creates a summary dataframe:
# counts depositor submissions, groups all but top 4 as "Other",
# and calculates each group's proportional contribution to the bin.

#Repeated calling function summary_depo for medium and long sequences below.
lyco_med <- lyco_cleaned[c(grep("Medium", lyco_cleaned$nuc_ranges)), ]
med_by_depo <- summary_depo(lyco_med, "Medium") 
head(med_by_depo)

lyco_long <- lyco_cleaned[c(grep("Long", lyco_cleaned$nuc_ranges)), ]
long_by_depo <- summary_depo(lyco_long, "Long")
head(long_by_depo)

# Combine all bins into one dataframe for plotting
plot_df <- bind_rows(short_by_depo, med_by_depo, long_by_depo)
# Combining bins this way ensures proportions remain relative to their bin totals,
# preserving the accuracy of “Other” vs. major depositors.

summary(plot_df)
head(plot_df)

# ===================================== #
#  Step 4: plotting results             #
# ===================================== #


## ===================================== ##
##  Graph 1: Bar graph                   

#total records contributed per depositor and bin
ggplot(plot_df, aes(x = nuc_ranges, y = n_records, fill = inst_grouped)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  scale_fill_viridis_d(option = "turbo") +
  labs(title = "Nucleotide Sequence Length Distribution by Depository",
       x = "Sequence Range (bp)",
       y = "Number of Records",
       fill = "Depository") +
  theme_minimal()

## ===================================== ##
##  Graph 2: Box Plot                    

top_5_markers <- lyco_cleaned %>%
  count(marker_code, sort = TRUE) %>%
  slice_max(n, n = 5) %>%
  pull(marker_code)

marker_plot_data <- lyco_cleaned %>%
  filter(marker_code %in% top_5_markers)

##--## LV Can add a statistical test here to test hypothesis
#The hypothesis: If dataset represents standardized barcoding practices, then majority of depositor submissions should cluster tightly around the expected COI barcode length (~660 bp)
#Plotting qqplot to determine distribution
qqnorm(marker_plot_data$nuc_scrubbed_count) 
qqline(marker_plot_data$nuc_scrubbed_count, col = "pink") #S shaped plot

##--## LV Since it is not normally distributed, we can use the Kruskal-Wallis test
ggplot(marker_plot_data, aes(x = marker_code, y = nuc_scrubbed_count, fill = marker_code)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
  scale_fill_viridis_d(option = "magma") +
  stat_compare_means(method = "kruskal", label.y = 1400) +
  labs(
    title = "Sequence Length Distribution by Marker Type",
    x = "Marker Type",
    y = "Sequence Length (bp)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1))

##--## LV p-value is p<2.2e-16, so statistically should reject null hypothesis

## ===================================== ##
##  Graph 3: Heatmap                     

#Subset dataset to only include COI-5P marker records (COI-5P is the most common barcode region for Lycosidae and will be used to evaluate sequence-length quality across different data depositories)

lyco_coi <- lyco_cleaned %>%
  filter(marker_code == "COI-5P")

# Before defining bins, inspect the sequence length distribution. table(cut()) shows how many sequences fall within each 10-bp window between 600–700 bp. This helps decide logical breakpoints for binning below.
summary(lyco_coi$nuc_scrubbed_count)
table(cut(lyco_coi$nuc_scrubbed_count, 
          breaks = c(seq(from = 600, to = 650, by = 10), 655, seq(from = 660, to = 700, by = 10))))

breaks <- c(0, 500, 600, 650, 700, 800, Inf)
labels <- c("0–500","500–600","600–650", "650–700","700–800","800+")

# Identify the top 4 depositories by total number of COI-5P records. These will be highlighted as individual rows; all others grouped as "Other".
top4 <- lyco_coi %>%
  count(depo_short, sort = TRUE) %>%
  slice_head(n = 4) %>%
  pull(depo_short)

# Create the summarized dataframe used for plotting the heatmap: Assign each sequence to a length_bin using the cut() function. Create a new column 'inst_grouped' that labels each depositor as either one of the top four or as "Other". Count how many sequences each depositor contributed to each bin (n). Convert those counts into proportions within each depositor row (so rows will sum to 1). This normalizes differences in total submissions.
heatmap_df <- lyco_coi %>%
  mutate(length_bin = cut(nuc_scrubbed_count, breaks = breaks, labels = labels, right = FALSE),
         inst_grouped = ifelse(depo_short %in% top4, depo_short, "Other")) %>%
  group_by(inst_grouped, length_bin) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(inst_grouped) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# Plot the normalized data as a heatmap.
ggplot(heatmap_df, aes(x = length_bin, y = inst_grouped, fill = prop)) +
  geom_tile(color = "grey80") +
  scale_fill_viridis(option = "cividis", name = "Proportion\nwithin group", limits = c(0, 1)) +
  labs(
    title = "COI-5P Sequence Length Distribution by Depository (Top 4 + Other)",
    subtitle = "Each row normalized by total submissions (proportion within depository)",
    x = "Sequence Length Bin (bp)",
    y = "Depository") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid = element_blank(), plot.title = element_text(face = "bold"))

#*********************** More plots (for interest mostly, I made these while deciding my project direction, and then ultimately chose other plots that I felt answered my research question better) 

# Graph 2: Proportional counts — depositor contribution within each bin
ggplot(plot_df, aes(x = nuc_ranges, y = prop, fill = inst_grouped)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  labs(title = "Relative Contribution by Depository",
       x = "Nucleotide Length Bin",
       y = "Proportion of Bin",
       fill = "Depository" ) +
  scale_fill_viridis_d(option = "plasma") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal()

# Graph 3: Distribution of sequence length bins per depositor
# Shows each top depositor’s internal submission proportions across bins

top_10_depos <- lyco_cleaned %>%
  count(inst, sort = TRUE) %>%
  slice_max(n, n = 10) %>%
  pull(inst)

depo_props <- lyco_cleaned %>%
  filter(inst %in% top_10_depos) %>% 
  count(inst, nuc_ranges, name = "n") %>%
  group_by(inst) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

ggplot(depo_props, aes(x = reorder(inst, -prop), y = prop, fill = nuc_ranges)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_viridis_d(option = "cividis") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Proportion of Sequence Length Submissions by Depositor",
       x = "Depositor",
       y = "Percentage of Submissions",
       fill = "Length Bin (nt)" ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 65, hjust = 1))

