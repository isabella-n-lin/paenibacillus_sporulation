# load packages
library(ggplot2)
library(ggtree)
library(ggnewscale)
library(phangorn)
library(readxl)
library(bio3d)
library(ggbreak)
library(reshape)
library(tidyverse)
library(ggheatmap)



### --- Figure 2: Phylogenetic Tree --- ###

# Read in Spo0B-TM data
spo0T_TMRs = data.frame(read_excel("~/Library/CloudStorage/Box-Box/Feaga Lab/Cassidy Prince/Bella/Bioinformatics for Paper/predicted_tmrs.xlsx"))
rownames(spo0T_TMRs) = spo0T_TMRs$assembly
spo0T_TMRs$TMRs = ifelse(spo0T_TMRs$Number_of_predicted_TMRs == 0, "Spo0B", "Spo0B-TM")

# Read in Spo0A and Spo0F data
spo0A_spo0F = data.frame(read_excel("~/Library/CloudStorage/Box-Box/Feaga Lab/Cassidy Prince/Bella/Bioinformatics for Paper/spo0F_spo0A.xlsx"))
rownames(spo0A_spo0F) = spo0A_spo0F$assembly
spo0A_spo0F[, 2] <- ifelse(spo0A_spo0F[, 2] == 1, "Spo0A", "absent")
spo0A_spo0F[, 3] <- ifelse(spo0A_spo0F[, 3] == 1, "Spo0F", "absent")

# Import unrooted GToTree of 1460 Paenibacillus genomes
unrooted_tree = read.tree("~/Library/CloudStorage/Box-Box/Feaga Lab/Cassidy Prince/Bella/Bioinformatics for Paper/gto_paeni_tree.tre")

# Midpoint root the tree
tree = midpoint(unrooted_tree)

# Modify the tip labels so they are NCBI RefSeq Assembly database accessions (beginning wih "GCF_").
tips = as.character(lapply(tree$tip.label, function(x) {
  tip_split = strsplit(x, "_")[[1]]
  paste(tip_split[1], tip_split[2], sep = "_")}))
tree$tip.label = tips

# Drop the genome GCF_047444735.1 becasue the genome was suppressed.
tree <- drop.tip(tree, "GCF_047444735.1")

#Plot tree with gene presence/absence
p = ggtree(tree, layout="circular", size = 0.1)  
p2 = p + new_scale_fill() 
p3 = gheatmap(p2, spo0A_spo0F["spo0F"], colnames = FALSE, width = 0.07, color=NA) + 
  scale_fill_manual(values=c("Spo0F" = "#006400", "absent" = "gray"), name = "", 
                    guide = guide_legend(order = 1))
p4 = p3 + new_scale_fill()
p5 = gheatmap(p4, spo0T_TMRs["TMRs"], colnames = FALSE, width = 0.07, color=NA, offset=0.04) +
  scale_fill_manual(values=c("Spo0B-TM" = "#27408B", "Spo0B" = "#6F8ED1"), name = "", 
                    breaks = c("Spo0B-TM", "Spo0B"), guide = guide_legend(order = 2)) 
p6 = p5 + new_scale_fill()
p7 = gheatmap(p6, spo0A_spo0F["spo0A"], colnames = FALSE, width = 0.07, color=NA, offset = 0.08) +
  scale_fill_manual(values=c("Spo0A" = "#542788", "absent" = "gray"), name = "", 
                    guide = guide_legend(order = 3)) 
p8 <- p7 + theme(
  legend.text = element_text(size = 10),
  legend.title = element_text(size = 12))

ggsave("~/Library/CloudStorage/Box-Box/Feaga Lab/Cassidy Prince/Bella/Bioinformatics for Paper/BLASTp/paeni_tree.png", p8, dpi=300, width=7, units="in")



### --- Figure 3A: Hydropathy Plot --- ###

# Import Kyte-Doolittle hyrdropathy plot values (from ProtScale on Expasy server)
hydrophobicity <- read_xlsx("~/Library/CloudStorage/Box-Box/Feaga Lab/Cassidy Prince/Bella/Bioinformatics for Paper/hydrophobicity_values.xlsx")

ggplot(data = hydrophobicity, aes(x=Position, y=Score)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  theme_classic()  +
  scale_x_continuous(expand = c(0, 0)) +
  labs(y="Hydropathy Score") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) + 
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))

ggsave("~/Library/CloudStorage/Box-Box/Feaga Lab/Cassidy Prince/Bella/Bioinformatics for Paper/hydrophobicity_plot.png", dpi=300, width = 4.5, units = "in")



### --- Figure 4B: Paenibacillus B-gal plot --- ###

# Import B-gal activity
df_activity <- read_xlsx("~/Library/CloudStorage/Box-Box/Feaga Lab/Cassidy Prince/Bella/BACTH/B-gal/all_paeni_bgal.xlsx")

# Calculate mean and SEM of B-gal activity for each transformation
activity <- df_activity %>% 
  group_by(transformation) %>%
  summarize(mean_B_gal_activity = mean(B_gal_activity),
    sem_B_gal_activity  = sd(B_gal_activity) / sqrt(n()),
    .groups = "drop") %>%
  mutate(protein = ifelse(transformation %in% 1:8, "Spo0A", "Spo0F"),
    transformation = factor(transformation, levels = c(4, 7, 5, 6, 8, 12, 15, 13, 14, 16)),
    protein = factor(protein, levels = c("Spo0F", "Spo0A")))

# Spo0F statistics
anova_spo0F <- aov(B_gal_activity ~ as.factor(transformation), 
                   data = df_activity %>% filter(transformation %in% c(12, 13, 14, 15, 16))) 
tukey_spo0F <- TukeyHSD(anova_spo0F)
print(tukey_spo0F)

# Spo0A statistics 
anova_spo0A <- aov(B_gal_activity ~ as.factor(transformation), 
                   data = df_activity %>% filter(transformation %in% c(4, 5, 6, 7, 8)))
tukey_spo0A <- TukeyHSD(anova_spo0A)
print(tukey_spo0A)

# Plotting function
plot_activity <- function(df, protein_name, colors, labels, annot_label) {
  max_y <- max(df$mean_B_gal_activity + df$sem_B_gal_activity)
  
  ggplot(df, aes(x = transformation, y = mean_B_gal_activity, fill = transformation)) +
    geom_col() +
    geom_errorbar(aes(ymin = mean_B_gal_activity - sem_B_gal_activity,
                      ymax = mean_B_gal_activity + sem_B_gal_activity),
                  width = 0.2) +
    scale_fill_manual(values = colors) +
    scale_y_continuous(expand = expansion(mult = 0)) +
    scale_x_discrete(labels = labels) +
    coord_cartesian(ylim = c(0, max_y * 1.1)) +
    ylab("Miller units") +
    theme_classic() +
    theme(
      text = element_text(size = 12),
      legend.position = "none",
      axis.text = element_text(color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
      axis.title.x = element_blank()) +
    # Use annotate instead of geom_segment
    annotate("segment", x = 0.6, xend = 3.4, y = max_y * 1.03, yend = max_y * 1.03, size = 1) +
    annotate("text", x = 2,   y = max_y * 1.07, label = paste0(protein_name, "-T25"), size = 4.25) +
    annotate("segment", x = 3.6, xend = 5.4, y = max_y * 1.03, yend = max_y * 1.03, size = 1) +
    annotate("text", x = 4.5, y = max_y * 1.07, label = "T25", size = 4.25)}

# Plot
spo0F_colors <- c("#006400", "#66a366", "gray", "gray", "gray")
spo0A_colors <- c("#542788", "#9a73c2", "gray", "gray", "gray")
x_labels <- c("Spo0B-TM-T18", "Spo0B-TM-∆TM-T18", "T18", "Spo0B-TM-T18", "Spo0B-TM-∆TM-T18")

p_spo0F <- plot_activity(filter(activity, protein == "Spo0F"), "Spo0F", spo0F_colors, x_labels, "Spo0F-T25")
p_spo0A <- plot_activity(filter(activity, protein == "Spo0A"), "Spo0A", spo0A_colors, x_labels, "Spo0A-T25") +
  theme(axis.title.y = element_blank())

combined_plot <- p_spo0F + p_spo0A

ggsave("~/Library/CloudStorage/Box-Box/Feaga Lab/Cassidy Prince/Bella/Bioinformatics for Paper/paeni_bgal_plot.png", combined_plot, dpi=300)



### --- Figure 6: Percent Identity Matrix --- ###

# Spo0F matrix
F_matrix <- as.data.frame(read_excel("~/Library/CloudStorage/Box-Box/Feaga Lab/Cassidy Prince/Bella/Bioinformatics for Paper/percent identity/spo0F_blast_matrix.xlsx"))
F_df <- melt(F_matrix)

names(F_df)[1] <- "species_1"
names(F_df)[2] <- "species_2"
names(F_df)[3] <- "identity_val"

F_df$species_1 = with(F_df, factor(species_1, levels = unique(F_df$species_1)))
F_df$species_2 = with(F_df, factor(species_2, levels = rev(unique(F_df$species_2))))

p <- ggplot(F_df, aes(x = species_1, y = species_2, fill = identity_val)) + 
  geom_tile() +
  scale_fill_gradient(low="white", high="#006400") +
  geom_text(aes(label = round(identity_val, digits = 1),
                color = ifelse(identity_val > 70, "white", "black")), size = 4) +
  scale_color_identity() +
  theme_minimal() +
  labs(fill = "Spo0F Percent Identity") +
  theme(text = element_text(size = 13),
        axis.text = element_text(face="italic"),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())

# Spo0B matrix
B_matrix <- as.data.frame(read_excel("~/Library/CloudStorage/Box-Box/Feaga Lab/Cassidy Prince/Bella/Bioinformatics for Paper/percent identity/spo0B_blast_matrix.xlsx"))
B_df <- melt(B_matrix)

names(B_df)[1] <- "species_1"
names(B_df)[2] <- "species_2"
names(B_df)[3] <- "identity_val"

B_df$species_1 = with(B_df, factor(species_1, levels = unique(B_df$species_1)))
B_df$species_2 = with(B_df, factor(species_2, levels = rev(unique(B_df$species_2))))

p2 <- ggplot(B_df, aes(x = species_1, y = species_2, fill = identity_val)) + 
  geom_tile() +
  scale_fill_gradient(low="white", high="#27408B") +
  geom_text(aes(label = round(identity_val, digits = 1),
                color = ifelse(identity_val > 70, "white", "black")), size = 4) +
  scale_color_identity() +
  theme_minimal() +
  labs(fill = "Spo0B Percent Identity") +
  theme(text = element_text(size = 13),
        axis.text = element_text(face="italic"),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())

# Spo0A matrix
A_matrix <- as.data.frame(read_excel("~/Library/CloudStorage/Box-Box/Feaga Lab/Cassidy Prince/Bella/Bioinformatics for Paper/percent identity/spo0A_blast_matrix.xlsx"))
A_df <- melt(A_matrix)

names(A_df)[1] <- "species_1"
names(A_df)[2] <- "species_2"
names(A_df)[3] <- "identity_val"

A_df$species_1 = with(A_df, factor(species_1, levels = unique(A_df$species_1)))
A_df$species_2 = with(A_df, factor(species_2, levels = rev(unique(A_df$species_2))))

p3 <- ggplot(A_df, aes(x = species_1, y = species_2, fill = identity_val)) + 
  geom_tile() +
  scale_fill_gradient(low="white", high="#542788") +
  geom_text(aes(label = round(identity_val, digits = 1),
                color = ifelse(identity_val > 70, "white", "black")), size = 4) +
  scale_color_identity() +
  theme_minimal() +
  labs(fill = "Spo0A Percent Identity") +
  theme(text = element_text(size = 13),
        axis.text = element_text(face="italic"),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())

# Combine plots
combined <- p / p2 / p3 

ggsave("~/Library/CloudStorage/Box-Box/Feaga Lab/Cassidy Prince/Bella/Bioinformatics for Paper/percent identity/combined_matrix.svg", combined, dpi = 300)



### --- Figure 7: Conservation Plot --- ###

# Impot MAFFT alignment
aln <- read.fasta("~/Library/CloudStorage/Box-Box/Feaga Lab/Cassidy Prince/Bella/Bioinformatics for Paper/Spo0B_align_op3_ep_3.fasta")

# Score conservation
scores <- conserv(x=aln$ali, method="identity", sub.matrix="bio3d")
df <- data.frame(position = seq_along(scores), conservation = scores)

cplot <- ggplot(df, aes(x = position, y = conservation)) +
  geom_area(mapping = aes(y= ifelse(conservation>0.75, conservation, 0)), fill = "#27408B") +
  geom_line() +
  labs(x = "Alignment Position", y = "Conservation") +
  scale_y_continuous(expand = expansion(mult = c(0,0), add = c(0, 0.04)), limits = c(0, 1.06)) +
  scale_x_continuous(expand = c(0,0)) +
  theme_classic() +
  annotate("segment", x = 120, xend = 127, y = 1.03, yend = 1.03, colour = "red", linewidth = 1.2) + # phosphorylatable region
  annotate("segment", x = 32, xend = 88, y = 0.7, yend = 0.7, colour = "#66A542", linewidth = 1.2) + # transmembrane region
  theme(text = element_text(size = 14),
        axis.text.x = element_text(size = 12, color="black"),
        axis.text.y = element_text(size = 12, color="black"))

ggsave("~/Library/CloudStorage/Box-Box/Feaga Lab/Cassidy Prince/Bella/Bioinformatics for Paper/conservation_plot.png", cplot, dpi = 300, width = 4.5, units = "in", height = 3.5)



### --- Figure S1 --- ###

# Import gene presence matrix
metadata_df <- read_excel("~/Library/CloudStorage/Box-Box/Feaga Lab/Cassidy Prince/Bella/Bioinformatics for Paper/metadata.xlsx")
metadata_df <- metadata_df[, -c(2:24)]
metadata_df <- t(metadata_df)
colnames(metadata_df) <- metadata_df[1,]
metadata_df <- metadata_df[-1,]
metadata_df <- tibble::rownames_to_column(as.data.frame(metadata_df), "gene")

# Shuffle the columns to simulate the order of sequencing
set.seed(32)
shuffled_columns <- sample(colnames(metadata_df)[-1])
shuffled_df <- metadata_df %>% select(c("gene", all_of(shuffled_columns)))

# Initialize variables to track new genes
cumulative_new_genes <- numeric()
unique_genes <- character()

for (genome in shuffled_columns) {
  new_genes <- setdiff(shuffled_df$gene[shuffled_df[[genome]] == 1], unique_genes)
  unique_genes <- union(unique_genes, new_genes)
  cumulative_new_genes <- c(cumulative_new_genes, length(new_genes))
}

# Plot
plot_data <- data.frame(
  Number_of_Genomes = 1:length(cumulative_new_genes),
  Number_of_New_Genes = cumulative_new_genes)

plot <- ggplot(plot_data, aes(x = Number_of_Genomes, y = Number_of_New_Genes)) +
  geom_point(size= 1, color="black") +
  geom_line(color="black") +
  scale_x_cut(108, which = c(1,2), scales=c(1, 0.9), space = 0.2) +
  scale_y_continuous(expand = c(0, 0)) + 
  labs(x = "Number of Genomes", y = "Number of New Genes") +
  theme_classic() +
  theme(axis.text = element_text(color="black", size = 12),
        axis.ticks = element_line(color = "black"),
        axis.text.x = element_text(hjust = 0.5),
        axis.labels = element_text(size = 12))

ggsave("~/Library/CloudStorage/Box-Box/Feaga Lab/Cassidy Prince/Bella/Bioinformatics for Paper/genes_per_seq_genome.png", plot, width = 7, dpi = 300, units = "in")



### --- Figure S3 --- ###

# Import B-gal activity
df_activity <- read_xlsx("~/Library/CloudStorage/Box-Box/Feaga Lab/Cassidy Prince/Bella/BACTH/B-gal/all_bs_bgal.xlsx")

# Calculate mean and SEM of B-gal activity for each transformation
activity <- df_activity %>% 
  group_by(transformation) %>%
  summarize(mean_B_gal_activity = mean(B_gal_activity),
            sem_B_gal_activity  = sd(B_gal_activity) / sqrt(n()),
            .groups = "drop") %>%
  mutate(species = ifelse(transformation %in% c(1,2,3,9,10,11), "B. subtilis", "P. polymyxa"),
         protein = ifelse(transformation %in% 1:8, "Spo0A", "Spo0F"),
         transformation = factor(sort(as.numeric(transformation))),
         protein = factor(protein, levels = c("Spo0F", "Spo0A")))

# Spo0F statistics
anova_spo0F <- aov(B_gal_activity ~ as.factor(transformation), 
                   data = df_activity %>% filter(transformation %in% c(9, 10, 11))) 
tukey_spo0F <- TukeyHSD(anova_spo0F)
print(tukey_spo0F)

# Spo0A statistics 
anova_spo0A <- aov(B_gal_activity ~ as.factor(transformation), 
                   data = df_activity %>% filter(transformation %in% c(1, 2, 3)))
tukey_spo0A <- TukeyHSD(anova_spo0A)
print(tukey_spo0A)

# Plotting function
plot_activity <- function(df, colors, x_labels, protein_name) {
  max_y <- max(df$mean_B_gal_activity + df$sem_B_gal_activity)
  
  ggplot(df, aes(x = transformation, y = mean_B_gal_activity, fill = transformation)) +
    geom_col() +
    geom_errorbar(
      aes(ymin = mean_B_gal_activity - sem_B_gal_activity,
          ymax = mean_B_gal_activity + sem_B_gal_activity),
      width = 0.2) +
    scale_fill_manual(values = colors) +
    scale_y_continuous(expand = expansion(mult = 0)) +
    scale_x_discrete(labels = x_labels) +
    coord_cartesian(ylim = c(0, max_y * 1.1)) +
    ylab("Miller units") +
    theme_classic() +
    theme(
      text = element_text(size = 12),
      legend.position = "none",
      axis.text = element_text(color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
      axis.title.x = element_blank()) +
    # Annotate
    annotate("segment", x = 0.6, xend = 2.4, y = max_y * 1.03, yend = max_y * 1.03, size = 1) +
    annotate("text", x = 1.5, y = max_y * 1.06, label = paste0(protein_name, "-T25"), size = 3.5) +
    annotate("segment", x = 2.6, xend = 3.4, y = max_y * 1.03, yend = max_y * 1.03, size = 1) +
    annotate("text", x = 3, y = max_y * 1.06, label = "T25", size = 3.5)}


#Plot
spo0F_colors <- c("#006400", "gray", "gray")
spo0A_colors <- c("#542788", "gray", "gray")

spo0F_labels <- c("T18-Spo0B", "T18", "T18-Spo0B")
spo0A_labels <- c("Spo0B-T18", "T18", "Spo0B-T18")

p_spo0F <- plot_activity(filter(activity, protein == "Spo0F"), spo0F_colors, spo0F_labels, "Spo0F")
p_spo0A <- plot_activity(filter(activity, protein == "Spo0A"), spo0A_colors, spo0A_labels, "Spo0A") +
  theme(axis.title.y = element_blank())

combined_plot <- p_spo0F + p_spo0A
combined_plot

ggsave("~/Library/CloudStorage/Box-Box/Feaga Lab/Cassidy Prince/Bella/Bioinformatics for Paper/bsub_bgal_plot.png",combined_plot, dpi = 300, width = 4.5, units = "in")

