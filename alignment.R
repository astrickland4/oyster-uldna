# Load packages
library(ggplot2)
library(dplyr)

theme_ecogen <- function(){
  theme_classic()+
    theme(axis.title = element_text(size = 14), 
          axis.text.x = element_text(size=10, color="black"), 
          axis.ticks = element_line (colour = "black", size = 1),  
          axis.line = element_line(color = "black", size = 1),
          axis.text.y = element_text(size = 14, color = 'black'))
}

# Read the data
coords <- read.table("alignment2.coords.txt", header = FALSE, fill = TRUE, sep = "\t",
                     col.names = c("S1", "E1", "S2", "E2", "LEN1", "LEN2", "IDY", "LEN_R", "LEN_Q", 
                                   "COV_R", "COV_Q", "ref", "query"))

coords <- coords %>%
  mutate(across(c(S1, E1, S2, E2, LEN1, LEN2, IDY, LEN_R, LEN_Q, COV_R, COV_Q), as.numeric))


# Optional: filter or clean
coords <- coords %>% filter(LEN1 > 1000)  # keep long alignments only

coords <- coords %>%
  mutate(ref_mid = (S1 + E1) / 2,
         query_mid = (S2 + E2) / 2,
         direction = ifelse((E1 - S1)*(E2 - S2) > 0, "forward", "reverse"))

ggplot(coords, aes(x = ref_mid, y = query_mid, color = direction)) +
  geom_point(size = 0.5, alpha = 0.6) +
  facet_grid(query ~ ref, scales = "free", space = "free") +
  theme_bw() +
  theme(strip.text.y = element_text(angle = 0)) +
  labs(x = "Reference Genome Position", y = "Query Contig Position")

ggplot(coords, aes(x = S1, xend = E1, y = S2, yend = E2, color = direction)) +
  geom_segment(size = 0.3, alpha = 0.6) +
  facet_grid(query ~ ref, scales = "free", space = "free") +
  theme_bw() +
  labs(x = "Reference Genome Position", y = "Query Contig Position")


# Rename REF and QUERY if needed
colnames(coords)[12:13] <- c("ref", "query")
coords_subset <- coords %>% filter(ref == "CM008241.1")


# Plot: dot plot with REF on x-axis and QUERY on y-axis
ggplot(coords_subset, aes(x = ref_mid, y = query_mid, color = direction)) +
  geom_point(alpha = 0.4, size = 0.5) +
  facet_grid(query ~ ref, scales = "free", space = "free") +
  scale_color_manual(values = c("forward" = "salmon", "reverse" = "skyblue")) +
  labs(
    title = "Genome Alignment: Hifiasm vs Reference",
    x = "Reference Genome Position",
    y = "Query Contig Position"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size = 6),
    panel.spacing = unit(0.1, "lines"),
    legend.position = "right"
  )

library(tidyverse)

coords %>%
  mutate(
    direction = ifelse((E1 - S1)*(E2 - S2) > 0, "forward", "reverse"),
    ref_label = factor(ref, levels = unique(ref)),
    query_label = factor(query, levels = unique(query))
  ) %>%
  ggplot(aes(x = ref_label, y = query_label, color = direction)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_color_manual(values = c("forward" = "#FF6666", "reverse" = "#66B2FF")) +
  labs(
    title = "Genome Alignment: Hifiasm Contigs vs Reference Chromosomes",
    x = "Reference Chromosome",
    y = "Hifiasm Contig"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

coords %>%
  mutate(
    direction = ifelse((E1 - S1)*(E2 - S2) > 0, "forward", "reverse"),
    query = factor(query, levels = unique(query))
  ) %>%
  ggplot(aes(x = S1, xend = E1, y = S2, yend = E2, color = direction)) +
  geom_segment(alpha = 0.4) +
  facet_wrap(~ query, scales = "free_y") +
  scale_color_manual(values = c("forward" = "firebrick", "reverse" = "steelblue")) +
  labs(
    title = "Alignment Blocks per Contig",
    x = "Reference Position",
    y = "Query Position"
  ) +
  theme_bw(base_size = 10)

# Step 1: Get chromosome lengths
chrom_lengths <- coords %>%
  group_by(ref) %>%
  summarise(length = max(c(S1, E1), na.rm = TRUE)) %>%
  arrange(ref)  # optional: enforce order

# Step 2: Compute cumulative start positions
chrom_starts <- chrom_lengths %>%
  mutate(start = lag(cumsum(length), default = 0))

# Step 3: Create named vector of starts for mapping
chrom_starts_vec <- setNames(chrom_starts$start, chrom_starts$ref)


coords %>%
  mutate(
    direction = ifelse((E1 - S1)*(E2 - S2) > 0, "forward", "reverse"),
    ref_mid = (S1 + E1) / 2,
    contig_row = as.integer(factor(query))  # numeric contig row
  ) %>%
  ggplot(aes(x = ref_mid, y = contig_row, color = direction)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("forward" = "#FF6666", "reverse" = "#66B2FF")) +
  labs(
    title = "Genome Alignment: Hifiasm Contigs vs Reference (bp scale)",
    x = "Reference Genome Position (bp)",
    y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank()
  )+ geom_vline(xintercept = chrom_starts$start, color = "gray80", linetype = "dotted")


# Example structure
chrom_sizes <- data.frame(
  ref = c("CM008241.1", "CM008242.1", "CM008243.1", "CM008244.1", "CM008245.1", 
          "CM008246.1", "CM008247.1", "CM008248.1", "CM008249.1", "CM008250.1"),
  length = c(65668440, 61752955, 77061148, 59691872, 98698416, 
             51258098, 57830854, 75944018, 104168038, 32650045)  # actual chromosome sizes
)

chrom_sizes <- chrom_sizes %>%
  mutate(start = c(0, cumsum(length)[-length(length)]))

coords_offset <- coords %>%
  left_join(chrom_sizes, by = c("ref" = "ref")) %>%
  mutate(
    direction = ifelse((E1 - S1)*(E2 - S2) > 0, "forward", "reverse"),
    ref_mid = (S1 + E1) / 2 + start,
    contig_row = as.integer(factor(query))  # or `query`
  )

chrom_labels <- chrom_sizes %>%
  mutate(mid = start + length / 2)

x_labs <- c("Chr 1", "Chr 2", "Chr 3", "Chr 4", "Chr 5", "Chr 6", "Chr 7", "Chr 8", "Chr 9", "Chr 10")


ggplot(coords_offset, aes(x = ref_mid, y = 0.05, color = direction)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_vline(xintercept = chrom_sizes$start, color = "gray80", linetype = "dotted") +
  scale_color_manual(values = c("forward" = "#FF6666", "reverse" = "#66B2FF")) +
  labs(
    x = "Reference Genome Position (bp)",
    y = NULL,
    color = "Alignment Direction"
  ) +
  ylim(0,1)+
  theme_ecogen() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank()
  )+ scale_x_continuous(
    breaks = chrom_labels$mid,
    labels = x_labs
  )

ggsave("alignment.png")


