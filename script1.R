library(dplyr)
library(readr)

# PAF has no headers; specify manually
cols <- c("query", "qlen", "qstart", "qend", "strand", "target", "tlen", "tstart", "tend", "matches", "aln_len", "mapq")

paf <- read_tsv("alignment.paf", col_names = cols)

# Optional: Filter by alignment length or map quality
filtered <- paf %>%
  filter(aln_len > 10000, mapq > 20) %>%
  mutate(y_query = 1, y_ref = 2)

summary(filtered$qlen)

sum(is.na(filtered$qlen))
sum(filtered$qlen == 0)

max_qlen <- max(filtered$qlen, na.rm = TRUE)
max_tlen <- max(filtered$tlen, na.rm = TRUE)

ggplot(filtered) +
  # Reference bar
  geom_segment(aes(x = 0, xend = max_tlen, y = 2, yend = 2), size = 4, color = "skyblue") +
  # Assembly bar (fix here)
  geom_segment(aes(x = 0, xend = max_qlen, y = 1, yend = 1), size = 4, color = "orange") +
  # Curves
  geom_curve(aes(x = qstart, xend = tstart, y = 1, yend = 2), curvature = -0.2, alpha = 0.4) +
  geom_curve(aes(x = qend, xend = tend, y = 1, yend = 2), curvature = 0.2, alpha = 0.4) +
  theme_minimal() +
  xlab("Position") + ylab("") +
  scale_y_continuous(breaks = NULL) +
  ggtitle("Assembly-to-Reference Alignment")



library(ggplot2)

ggplot(filtered) +
  # Reference genome bar
  geom_segment(aes(x = 0, xend = max(tlen), y = 2, yend = 2), size = 4, color = "skyblue") +
  
  # Assembly genome bar
  geom_segment(aes(x = 0, xend = max(qlen), y = 1, yend = 1), size = 4, color = "orange") +
  
  # Aligned segments as curves
  geom_curve(aes(x = qstart, xend = tstart, y = 1, yend = 2), curvature = -0.2, alpha = 0.4) +
  geom_curve(aes(x = qend, xend = tend, y = 1, yend = 2), curvature = 0.2, alpha = 0.4) +
  
  theme_minimal() +
  xlab("Position") + ylab("") +
  scale_y_continuous(breaks = NULL) +
  ggtitle("Assembly-to-Reference Alignment")
