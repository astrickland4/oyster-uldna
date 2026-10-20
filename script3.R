library(readr)
library(dplyr)
library(ggplot2)

# Load PAF with appropriate columns
cols <- c("query", "qlen", "qstart", "qend", "strand",
          "target", "tlen", "tstart", "tend",
          "nmatch", "alnlen", "mapq")

paf <- read_tsv("alignment.paf", col_names = cols)

# Filter alignments to remove very short ones
filtered <- paf %>%
  filter(alnlen > 10000, mapq >= 20)

# Get chromosome (reference target) lengths
chrom_lengths <- filtered %>%
  group_by(target) %>%
  summarise(chr_len = max(tend, na.rm = TRUE)) %>%
  mutate(start = cumsum(lag(chr_len, default = 0)))

# Join to get offset for plotting
filtered_offset <- filtered %>%
  left_join(chrom_lengths, by = "target") %>%
  mutate(
    tstart_adj = tstart + start,
    tend_adj = tend + start
  )

# Plot
ggplot(filtered_offset) +
  geom_segment(aes(x = tstart_adj, xend = tend_adj, y = query, yend = query),
               color = "steelblue", alpha = 0.8, size = 1) +
  geom_vline(data = chrom_lengths, aes(xintercept = start), color = "gray80") +
  scale_x_continuous(
    breaks = chrom_lengths$start + chrom_lengths$chr_len / 2,
    labels = chrom_lengths$target
  ) +
  labs(
    title = "Contig Coverage Across Reference Genome",
    x = "Reference Genome (Chromosomes)",
    y = "Hifiasm Contigs"
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

summary(filtered_offset$tstart_adj)
summary(filtered_offset$tend_adj)
summary(filtered_offset$tend_adj - filtered_offset$tstart_adj)

filtered_offset$query <- factor(filtered_offset$query, levels = unique(filtered_offset$query))

top_contigs <- filtered_offset %>%
  count(query, sort = TRUE) %>%
  slice(1:30) %>%  # pick top 30 contigs
  pull(query)

ggplot(filtered_offset %>% filter(query %in% top_contigs)) +
  geom_segment(aes(x = tstart_adj, xend = tend_adj, y = query, yend = query),
               color = "steelblue", size = 1.2) +
  geom_vline(data = chrom_lengths, aes(xintercept = start), color = "gray80") +
  scale_x_continuous(
    breaks = chrom_lengths$start + chrom_lengths$chr_len / 2,
    labels = chrom_lengths$target
  ) +
  theme_bw() +
  labs(
    title = "Top 30 Contig Alignments Across Reference Genome",
    x = "Reference Genome (Chromosomes)", y = "Hifiasm Contigs"
  ) +
  theme(
    axis.text.y = element_text(size = 7),
    panel.grid.major.y = element_blank()
  )


ggplot(filtered_offset %>% filter(query == "ptg000001l")) +
  geom_segment(aes(x = tstart_adj, xend = tend_adj, y = 1, yend = 1), size = 3)


ggplot(filtered_offset) +
  geom_segment(aes(x = tstart_adj, xend = tend_adj, y = 1, yend = 1), alpha = 0.2) +
  geom_vline(data = chrom_lengths, aes(xintercept = start), color = "gray80") +
  scale_x_continuous(
    breaks = chrom_lengths$start + chrom_lengths$chr_len / 2,
    labels = chrom_lengths$target
  ) +
  labs(
    title = "Assembly Coverage Across Reference Genome",
    x = "Reference Genome", y = NULL
  ) +
  theme_bw() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

ggplot(filtered_offset) +
  geom_segment(aes(x = tstart_adj, xend = tend_adj, y = query, yend = query)) +
  facet_wrap(~ target, scales = "free_x") +
  theme_minimal()


# Make sure query is a factor for y-axis ordering
filtered_offset$query <- factor(filtered_offset$query, levels = unique(filtered_offset$query))

# Plot top 30 contigs
top_contigs <- filtered_offset %>%
  group_by(query) %>%
  summarise(total_aligned = sum(tend_adj - tstart_adj)) %>%
  arrange(desc(total_aligned)) %>%
  slice_head(n = 30) %>%
  pull(query)

ggplot(filtered_offset %>% filter(query %in% top_contigs)) +
  geom_segment(aes(x = tstart_adj, xend = tend_adj, y = query, yend = query, color = target),
               size = 1.5, alpha = 0.9) +
  geom_vline(data = chrom_lengths, aes(xintercept = start), color = "gray80") +
  scale_x_continuous(
    breaks = chrom_lengths$start + chrom_lengths$chr_len / 2,
    labels = chrom_lengths$target
  ) +
  labs(
    title = "Top 30 Hifiasm Contigs Aligned to Reference Genome",
    x = "Reference Chromosomes", y = "Contigs"
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

ggplot(filtered_offset %>% filter(query %in% top_contigs)) +
  geom_segment(aes(x = tstart_adj, xend = tend_adj, y = query, yend = query),
               size = 2.5, color = "steelblue", lineend = "round") +
  geom_vline(data = chrom_lengths, aes(xintercept = start), color = "gray80") +
  scale_x_continuous(
    breaks = chrom_lengths$start + chrom_lengths$chr_len / 2,
    labels = chrom_lengths$target
  ) +
  labs(
    title = "Top 30 Hifiasm Contigs Aligned to Reference Genome",
    x = "Reference Chromosomes", y = "Contigs"
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggplot(filtered_offset) +
  geom_segment(
    aes(x = tstart_adj, xend = tend_adj, y = 1, yend = 1),
    color = "steelblue",  # or "black"
    size = 2,             # thickness of lines
    alpha = 1             # fully opaque
  ) +
  geom_vline(
    data = chrom_lengths,
    aes(xintercept = start),
    color = "gray80",
    linewidth = 0.5
  ) +
  scale_x_continuous(
    breaks = chrom_lengths$start + chrom_lengths$chr_len / 2,
    labels = chrom_lengths$target
  ) +
  coord_cartesian(ylim = c(0.95, 1.05)) +
  scale_y_continuous(expand = c(0, 0))+
  labs(
    title = "Hifiasm Contigs Aligned to Reference Genome",
    x = "Reference Genome (Chromosomes)", y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )




    
