library(readr)
library(dplyr)
library("ggsci")
library("ggplot2")
library("gridExtra")

library(ggrepel)
library(forcats)
library(stringr)

KEGG_up <- read.csv("KEGG-UP.csv", sep = '\t')
KEGG_down <- read.csv("KEGG-DOWN.csv", sep = '\t')
KEGG_full <- rbind(KEGG_up, KEGG_down)

KEGG_full$Direction <- c(rep("UP", nrow(KEGG_up)), rep("DOWN", nrow(KEGG_down)) )

KEGG_full <- KEGG_full[, c("Term", "Count", "PValue", "Direction")]

#KEGG_full$Term <- str_replace(KEGG_full$Term, 'dre\\d{5}:', '')
KEGG_full$negLogPValue <- -log(KEGG_full$PValue, base = 10)
head(KEGG_full)

KEGG_full$Direction <- factor(KEGG_full$Direction, levels = c("DOWN", "UP"))
KEGG_full <-KEGG_full %>% arrange(Direction, negLogPValue)
values <- as.vector(KEGG_full$Term)
KEGG_full$Term <- factor(KEGG_full$Term, levels = values)


ggplot(data = KEGG_full, aes(x = Term, y = negLogPValue, fill = Direction)) +
  geom_col() +
  scale_fill_manual(values=c(  "#56B4E9", "red")) +
  ylab("Negative Log10 of P-value") +
  coord_flip()+
  theme_classic() +
  theme(axis.text.x = element_text( size = 16),
        axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 17),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank()
  )

ggsave("KEGG_pathways.png", dpi = 300, height = 10, width = 10)

