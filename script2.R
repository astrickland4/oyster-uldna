library(readr)
library(dplyr)
library(ggplot2)

# Define column names for PAF format
cols <- c("query", "qlen", "qstart", "qend", "strand",
          "target", "tlen", "tstart", "tend", 
          "nmatch", "alnlen", "mapq")

# Load the PAF file
paf <- read_tsv("alignment.paf", col_names = cols)

# Optional: filter good alignments
filtered <- paf %>%
  filter(alnlen > 10000, mapq >= 20) %>%
  mutate(query = factor(query, levels = unique(query)))  # keep contigs in order

ggplot(filtered) +
  geom_segment(aes(x = tstart, xend = tend,
                   y = query, yend = query,
                   color = target),
               size = 1.2, alpha = 0.8) +
  scale_x_continuous(labels = scales::comma) +
  facet_wrap(~ target, scales = "free_x", nrow = 1) +  # one panel per chromosome
  labs(title = "Contig Coverage Across Reference Genome",
       x = "Reference Genome Position",
       y = "Hifiasm Contigs",
       color = "Chromosome") +
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(size = 10),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.spacing = unit(1, "lines")
  )

# OPTIONAL: simplify to a continuous x-axis with labeled chromosomes
chrom_order <- unique(filtered$target)
chrom_lengths <- filtered %>%
  group_by(target) %>%
  summarise(chr_len = max(tend))

chrom_starts <- chrom_lengths %>%
  mutate(start = cumsum(lag(chr_len, default = 0)))

# Add offset to coordinates for genome-wide plotting
filtered_offset <- filtered %>%
  left_join(chrom_starts, by = c("target")) %>%
  mutate(tstart_adj = tstart + start,
         tend_adj = tend + start)

ggplot(filtered_offset) +
  geom_segment(aes(x = tstart_adj, xend = tend_adj,
                   y = query, yend = query),
               size = 1.2, color = "steelblue") +
  geom_vline(data = chrom_starts, aes(xintercept = start), color = "gray80") +
  scale_x_continuous(
    breaks = chrom_starts$start + chrom_lengths$chr_len / 2,
    labels = chrom_starts$target
  ) +
  theme_bw() +
  labs(x = "Reference Genome (chromosomes)", y = "Hifiasm Contigs") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

ggplot(filtered_offset) +
  geom_tile(aes(x = (tstart_adj + tend_adj)/2, y = query, width = tend_adj - tstart_adj, height = 0.5),
            fill = "steelblue", alpha = 0.8) +
  theme_minimal() +
  labs(x = "Reference Genome", y = "Contigs")


head(filtered_offset[, c("target", "tstart", "tend", "start", "tstart_adj", "tend_adj")])
unique(filtered$target)
unique(chrom_starts$target)
filtered$target <- trimws(filtered$target)
chrom_starts$target <- trimws(chrom_starts$target)
filtered_offset %>% filter(query %in% top_contigs$query)


length(top_contigs$query)
top_contigs <- filtered %>%
  group_by(query) %>%
  summarise(total_aligned = sum(alnlen)) %>%
  arrange(desc(total_aligned)) %>%
  slice(1:30)  # adjust as needed

filtered_offset <- filtered_offset %>%
  filter(query %in% top_contigs$query)

summary(filtered_offset$tend_adj - filtered_offset$tstart_adj)

ggplot(filtered) +
  geom_segment(aes(x = tstart, xend = tend,
                   y = query, yend = query,
                   color = target),
               size = 1.2, alpha = 0.8) +
  labs(title = "Contig Alignments on Reference",
       x = "Reference Position", y = "Contigs") +
  theme_minimal() +
  theme(axis.text.y = element_blank())

