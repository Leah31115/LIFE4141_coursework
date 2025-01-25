# Load in libraries
library(ggplot2)
library(tidyverse)


# This script makes plots of fst scans for each population comparison against 
# the LAB population to identify 1% empirical outliers

setwd("C:/Users/leahe/Documents/LIFE4141_Comparative_and_Evolutionary_Genomics/LIFE4141_coursework_resources/evo_data_files")

# FST scan for LAB vs NEN
fst_1000_win <- read.table("lab_nen_fst.1000.out.windowed.weir.fst", sep = "\t",
                           header = T)

# Sort by weighted fst
fst_1000_win %>%
  arrange(desc(WEIGHTED_FST))

# Make a mean bin column
fst_1000_win$MEAN_BIN = rowMeans(fst_1000_win[,c("BIN_START", "BIN_END")], na.rm=T)

# Add snp column.
# Assign a number to each SNP.
fst_1000_win_len <- dim(fst_1000_win)[1]
fst_1000_win$SNP <- paste('SNP', 1:fst_1000_win_len)

# Convert any negative weighted fst values to 0
fst_1000_win$WEIGHTED_FST[fst_1000_win$WEIGHTED_FST < 0] <- 0

# Fst outliers
# identify the 99% percentile
outlier_threshold <- quantile(fst_1000_win$WEIGHTED_FST, 0.99, na.rm = T)

# Make a new column of outliers
fst_1000_win$outlier <- ifelse(fst_1000_win$WEIGHTED_FST > outlier_threshold, "outlier", "background")

# Gene name labels with fst values greater than the 99th percentile
fst_high <- fst_1000_win %>%
  filter(WEIGHTED_FST >= outlier_threshold)

# Gene name labels for outliers
fst_outliers <- fst_1000_win %>%
  filter(outlier == "outlier")

# FST plot for LAB VS NEN
# Plot with  > 1% FST labels without outliers highlighted by coloured text
ggplot(fst_1000_win, aes(x=MEAN_BIN, y=WEIGHTED_FST)) +
  geom_point() +
  geom_point(data = fst_high, aes(x = MEAN_BIN, y = WEIGHTED_FST), color = "blue") +
  xlab("Chromosome 5 location (Mb)") +
  ylab("Weighted FST") +
  ggtitle("FST scan of chromosome 5 for LAB vs NEN Cochlearia pyrenaica populations") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = outlier_threshold, linetype="dashed", color="red") +
  annotate("rect", xmin = 4350000, xmax = 4500000, ymin = 0, ymax = 0.7,
           alpha = 0.15, fill = "grey") +
  annotate("rect", xmin = 4800000, xmax = 5100000, ymin = 0, ymax = 0.7,
           alpha = 0.15, fill = "grey") +
  annotate("rect", xmin = 5550000, xmax = 5750000, ymin = 0, ymax = 0.7,
           alpha = 0.15, fill = "grey") +
  annotate("rect", xmin = 5800000, xmax = 6100000, ymin = 0, ymax = 0.7,
           alpha = 0.15, fill = "grey") +
  annotate("rect", xmin = 6350000, xmax = 6450000, ymin = 0, ymax = 0.7,
           alpha = 0.15, fill = "grey") +
  geom_text(data = fst_high, aes(label=SNP), size=2.5, color = "blue") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) 

# Save outlier fst LAB NEN
lab_nen_out_fst <- fst_1000_win %>%
  filter(outlier == "outlier") %>%
  select("CHROM", "BIN_START", "BIN_END")
write.table(lab_nen_out_fst, "lab_nen_out_fst.bed", sep = "\t")


# For LAB vs ODN
fst_ODN_LAB <- read.table("ODN_LAB_out.windowed.weir.fst", sep="\t", header=T)

# Sort by weighted fst
fst_ODN_LAB %>%
  arrange(desc(WEIGHTED_FST))

# Add snp column.
# Assign a number to each SNP.
fst_ODN_LAB_len <- dim(fst_ODN_LAB)[1]
fst_ODN_LAB$SNP <- paste('SNP', 1:fst_ODN_LAB_len)

# Make a mean bin column
fst_ODN_LAB$MEAN_BIN = rowMeans(fst_ODN_LAB[,c("BIN_START", "BIN_END")], na.rm=T)

# Convert any negative fst values to 0
fst_ODN_LAB$WEIGHTED_FST[fst_ODN_LAB$WEIGHTED_FST < 0] <- 0

# Fst outliers
# Identify the 99% percentile, 1 tail
OL_outlier_threshold <- quantile(fst_ODN_LAB$WEIGHTED_FST, 0.99, na.rm = T)

# Make a new column of outliers
fst_ODN_LAB$outlier <- ifelse(fst_ODN_LAB$WEIGHTED_FST > OL_outlier_threshold, "outlier", "background")

# Gene name labels with fst values greater than 99th percentile
OL_fst_high <- fst_ODN_LAB %>%
  filter(WEIGHTED_FST >= OL_outlier_threshold)

# Gene name labels for outliers
OL_fst_outliers <- fst_ODN_LAB %>%
  filter(outlier == "outlier")

# Plot with  99% outlier labels without outliers highlighted
ggplot(fst_ODN_LAB, aes(x=MEAN_BIN, y=WEIGHTED_FST)) +
  geom_point() +
  geom_point(data = OL_fst_outliers, aes(x = MEAN_BIN, y = WEIGHTED_FST), color = "blue") +
  xlab("Chromosome 5 location (Mb)") +
  ylab("Weighted FST") +
  ggtitle("FST scan of chromosome 5 for LAB vs ODN Cochlearia pyrenaica populations") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = OL_outlier_threshold, linetype="dashed", color="red") +
  annotate("rect", xmin = 4200000, xmax = 4250000, ymin = 0, ymax = 1,
           alpha = 0.15, fill = "grey") +
  annotate("rect", xmin = 4350000, xmax = 4400000, ymin = 0, ymax = 1,
           alpha = 0.15, fill = "grey") +
  annotate("rect", xmin = 4450000, xmax = 4750000, ymin = 0, ymax = 1,
           alpha = 0.15, fill = "grey") +
  annotate("rect", xmin = 4900000, xmax = 5100000, ymin = 0, ymax = 1,
           alpha = 0.15, fill = "grey") +
  annotate("rect", xmin = 5250000, xmax = 5500000, ymin = 0, ymax = 1,
           alpha = 0.15, fill = "grey") +
  annotate("rect", xmin = 5600000, xmax = 5700000, ymin = 0, ymax = 1,
           alpha = 0.15, fill = "grey") +
  annotate("rect", xmin = 5900000, xmax = 6250000, ymin = 0, ymax = 1,
           alpha = 0.15, fill = "grey") +
  annotate("rect", xmin = 6650000, xmax = 6750000, ymin = 0, ymax = 1,
           alpha = 0.15, fill = "grey") +
  geom_text(data = OL_fst_high, aes(label = SNP), size = 2.5, color = "blue") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))


# Save outlier fst LAB ODN
lab_odn_out_fst <- fst_ODN_LAB %>%
  filter(outlier == "outlier") %>%
  select("CHROM", "BIN_START", "BIN_END")
write.table(lab_odn_out_fst, "lab_odn_out_fst.bed", sep = "\t")
