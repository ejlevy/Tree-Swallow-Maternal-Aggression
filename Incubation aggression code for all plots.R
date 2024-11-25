# Code for all figures


library(gridExtra)
library(readxl)
library(ggplot2)
library(dplyr)
library(UpSetR)
library(forcats)
library(ggExtra)


# Aggression score correlation (Figure 1) --------------------------------------------

df <- read_xlsx("Data frames for figures.xlsx", sheet = 1)

ggplot(df, aes(x=agg5min, y = agg30min)) + 
  geom_point(size = 3) + 
  xlim(0,60) + ylim(0,60) + 
  theme(panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor = element_line(colour = "grey90"),
        text = element_text(size = 13, colour = "black"), 
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 12, colour = "black")) +
  labs(x="5 Minute Aggression Score", y = "30 Mininute Aggression Score")
ggsave("Aggscores.jpeg", width = 4, height = 4, units= "in") 
ggsave("Aggscores.eps", width = 4, height = 4, units= "in") 



# Number of DEGs (Figure 2) -----------------------------------------------

## Colors:
# orange, upregulated: #F57340
# light orange, upregulated uncorrected DEGs: #F5B69A
# pink, downregulated: #B538A2
# light pink, downregulated uncorrected DEGs: #B77DAC

fig <- read_xlsx("Data frames for figures.xlsx", sheet = 2)

## Combined into one - DEGs only
p1 <- ggplot(subset(fig, Type == "DEG"), aes(x=reorder(Trait, Order), y=vSBN, 
                      fill = ifelse(vSBN>0 & Order > 1, "#B538A2","#F57340"))) +
  scale_fill_manual(values = c("#F57340", "#B538A2")) +
  geom_col(position = "stack") +
  ylim(-225, 225) + 
  geom_hline(aes(yintercept = 0)) +
  geom_text(size = 5, aes(label = ifelse(vSBN >= 0, vSBN, -vSBN), 
                          vjust = ifelse(vSBN >= 0, -.25, 1.15))) + 
  labs(x="", y = "Count \n\n", subtitle = "ventral SBN") +
  theme(panel.background = element_blank(), 
        text = element_text(size = 16, colour = "black"), 
        axis.text = element_text(size = 14, colour = "black"),
        axis.text.x = element_text(angle = 40, vjust = 1, hjust=1),
        #axis.text.y = element_blank(),
        #axis.ticks.y = element_blank(),
        legend.position = "none")

p2 <- ggplot(subset(fig, Type == "DEG"), aes(x=reorder(Trait, Order), y=dSBN, 
                      fill = ifelse(dSBN>0 & Order > 1, "#B538A2", "#F57340"))) +
  geom_col(position = "stack") +
  scale_fill_manual(values = c("#F57340", "#B538A2")) +
  ylim(-225, 225) + 
  geom_hline(aes(yintercept = 0)) +
  geom_text(size = 5, aes(label = ifelse(dSBN >= 0, dSBN, -dSBN), 
                          vjust = ifelse(dSBN >= 0, -.25, 1.15))) + 
  labs(x="", y = "\n\n", subtitle = "dorsal SBN") +
  guides(fill="none") +
  theme(panel.background = element_blank(), 
        text = element_text(size = 16, colour = "black"), 
        axis.text = element_text(size = 14, colour = "black"),
        axis.text.x = element_text(angle = 40, vjust = 1, hjust=1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")

grid.arrange(p1, p2, nrow = 1)
# I added the "upregulated" and "downregulated" labels in illustrator.

## Here are the transcripts associated with trait before correcting for multiple comparisons ( raw p < 0.05)
p3 <- ggplot(subset(fig, Type == "raw"), aes(x=reorder(Trait, Order), y=vSBN, 
                                             fill = ifelse(vSBN>0, "#B77DAC","#F5B69A"))) +
  scale_fill_manual(values = c("#F5B69A", "#B77DAC")) +
  geom_col(position = "stack") +
  ylim(-1000, 1400) + 
  geom_hline(aes(yintercept = 0)) +
  geom_text(size = 5, aes(label = ifelse(vSBN >= 0, vSBN, -vSBN), 
                          vjust = ifelse(vSBN >= 0, -.25, 1.15))) + 
  labs(x="", y = "Count \n\n", subtitle = "ventral SBN") +
  theme(panel.background = element_blank(), 
        text = element_text(size = 16, colour = "black"), 
        axis.text = element_text(size = 14, colour = "black"),
        axis.text.x = element_text(angle = 40, vjust = 1, hjust=1),
        #axis.text.y = element_blank(),
        #axis.ticks.y = element_blank(),
        legend.position = "none")
             
p4 <- ggplot(subset(fig, Type == "raw"), aes(x=reorder(Trait, Order), y=dSBN, 
                                             fill = ifelse(dSBN>0, "#B77DAC","#F5B69A"))) +
  scale_fill_manual(values = c("#F5B69A", "#B77DAC")) +
  geom_col(position = "stack") +
  ylim(-1000, 1400) + 
  geom_hline(aes(yintercept = 0)) +
  geom_text(size = 5, aes(label = ifelse(vSBN >= 0, dSBN, -dSBN), 
                          vjust = ifelse(vSBN >= 0, -.25, 1.15))) + 
  labs(x="", y = "\n\n", subtitle = "dorsal SBN") +
  theme(panel.background = element_blank(), 
        text = element_text(size = 16, colour = "black"), 
        axis.text = element_text(size = 14, colour = "black"),
        axis.text.x = element_text(angle = 40, vjust = 1, hjust=1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")
grid.arrange(p3, p4, nrow = 1)

  

# Upset plot (Figure 3) --------------------------------------------------------------

# Colors:
# blue, treatment uncorrected DEGs #297ED9
# lime green, aggression DEGs: #A4D444
# teal, overlapping DEGs: #3CBC9B

up <- read_xlsx("Data frames for figures.xlsx", sheet = 3)
up <- as.data.frame(up)

upset(up, sets = c("Aggression (5 min) dSBN", "Aggression (30 min) vSBN", "Aggression (5 min) vSBN", "Treatment dSBN*", "Treatment vSBN*"), 
      order.by = "freq", 
      main.bar.color = c("#297ED9", "#297ED9", "#297ED9", "#A4D444", "#A4D444", "#297ED9", "#3CBC9B", "#A4D444", 
                         "#3CBC9B", "#A4D444", "#3CBC9B", "#3CBC9B", "#3CBC9B", "#3CBC9B", "#3CBC9B", 
                         "#3CBC9B", "#A4D444", "#3CBC9B"),
      sets.bar.color = c("#A4D444", "#A4D444", "#A4D444", "#297ED9", "#297ED9"),
      keep.order = T, text.scale = 2) 



#####  Fisher's exact test to test whether overlap is greater than expected #### 
# Lower p-value = less likely to be by chance
# 14608 transcripts total

## Overlap between Treatment vSBN x Treatment dSBN
upset(up, sets = c("Treatment dSBN", "Treatment vSBN"), 
      order.by = "freq", keep.order = T, text.scale = 2) # Default version kicks out testosterone, cap latency, incubation day, and body mass which is fine.
# 14608 transcripts total
# Unique to Treatment vSBN = 313
# Unique to Treatment dSBN = 244
# Shared = 43

fisher.test(matrix(c((14608-313-244-43), 313, 244, 43), nrow=2), alternative="greater")
# Fisher's Exact Test for Count Data
# 
# data:  matrix(c((14608 - 313 - 244 - 43), 313, 244, 43), nrow = 2)
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
#  5.788384      Inf
# sample estimates:
# odds ratio 
#   7.884468 

# Way more overlapping that you would expect by chance


## Overlap between Agg5min vSBN x Agg30min dSBN
upset(up, sets = c("Agg5min vSBN", "Agg30min vSBN"), 
      order.by = "freq", keep.order = T, text.scale = 2) # Default version kicks out testosterone, cap latency, incubation day, and body mass which is fine.
# 14608 transcripts total
# Unique to Agg5 vSBN = 292
# Unique to Agg30 vSBN = 122
# Shared = 124

fisher.test(matrix(c((14608-292-122-124), 292, 122, 124), nrow=2), alternative="greater")
# Fisher's Exact Test for Count Data
# 
# data:  matrix(c((14608 - 292 - 122 - 124), 292, 122, 124), nrow = 2)
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
#  38.41916      Inf
# sample estimates:
# odds ratio 
#   48.90752 

# Waaaaay more overlapping that you would expect by chance

## Overlap between Treatment vSBN x Agg5min vSBN
upset(up, sets = c("Treatment vSBN", "Agg5min vSBN"), 
      order.by = "freq", keep.order = T, text.scale = 2) # Default version kicks out testosterone, cap latency, incubation day, and body mass which is fine.
# 14608 transcripts total
# Unique to Treatment vSBN = 339
# Unique to Agg5min vSBN = 399
# Shared = 17

fisher.test(matrix(c((14608-339-399-17), 339, 399, 17), nrow=2), alternative="greater")
# Fisher's Exact Test for Count Data
# 
# data:  matrix(c((14608 - 339 - 399 - 17), 339, 399, 17), nrow = 2)
# p-value = 0.02658
# alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
#  1.088021      Inf
# sample estimates:
# odds ratio 
#   1.740998 

# Relatively weak, but still more overlapping that we may expect by chance (using a 0.05 cut-off)


## Overlap between Treatment vSBN x Agg30min vSBN
upset(up, sets = c("Treatment vSBN", "Agg30min vSBN"), 
      order.by = "freq", keep.order = T, text.scale = 2) # Default version kicks out testosterone, cap latency, incubation day, and body mass which is fine.
# 14608 transcripts total
# Unique to Treatment vSBN = 351
# Unique to Agg30min vSBN = 241
# Shared = 5

fisher.test(matrix(c((14608-351-241-5), 351, 241, 5), nrow=2), alternative="greater")
# Fisher's Exact Test for Count Data
# 
# data:  matrix(c((14608 - 351 - 241 - 5), 351, 241, 5), nrow = 2)
# p-value = 0.7198
# alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
#  0.3217754       Inf
# sample estimates:
# odds ratio 
#  0.8281701 

# These two gene sets are not more overlapping that you would expect by chance.



# GO & Gene expression (Figure 4) -----------------------------------------

go <- read_xlsx("Data frames for figures.xlsx", sheet = 3)

## Colors:
# orange, upregulated: #F57340
# pink, downregulated: #B538A2

#### Main GO barplots ####


## 5min agg in vSBN DEGs -- POSITIVE LFC (upregulated in more aggressive females)
# Specific GO terms (at tips of nested terms)
# 5-minute aggression in vSBN \nDEGs with Positive LFC \nTop 10 Specific GO Terms
go %>% 
  filter(brain_region == "vSBN" & trait == "5minutePosLFC" & most_specific == "y") %>%
  arrange(desc(FDR_p)) %>%
  slice_tail(n=10) %>%
  mutate(GO_term=factor(GO_term, levels=GO_term)) %>% # This is a dumb line that makes the order  correct in the graph because my reorder functions within ggplot weren't working all of a sudden.
  ggplot(., aes(x=FDR_p, y = GO_term)) + 
  geom_point(color = "#F57340", size = 3) +
  xlim(0,0.05) +
  labs(x="P-Value (FDR)", y = "GO Term") +
  theme(panel.background = element_blank(), 
        text = element_text(size = 13, colour = "black"), 
        axis.text.y = element_text(size = 10, colour = "black", margin = margin(r = -8)),
        axis.text.x = element_text(size = 10, colour = "black"),
        axis.ticks.y = element_blank())
ggsave("vSBN 5min Positive LFC specific GO terms.jpeg", height = 4, width = 8.5, units = "in")
ggsave("vSBN 5min Positive LFC specific GO terms.eps", height = 4, width = 8.5, units = "in")
ggsave("vSBN 5min Positive LFC specific GO terms.pdf", height = 4, width = 8.5, units = "in")


## 5min agg in vSBN DEGs -- NEGATIVE LFC (upregulated in less aggressive females)
# Specific GO terms (at tips of nested terms)
# 5-minute aggression in vSBN \nDEGs with Negative LFC \nAll Specific GO Terms
go %>% 
  filter(brain_region == "vSBN" & trait == "5minuteNegLFC" & most_specific == "y") %>%
  arrange(desc(FDR_p)) %>%
  slice_tail(n=10) %>%
  mutate(GO_term=factor(GO_term, levels=GO_term)) %>% # This is a dumb line that makes the order  correct in the graph because my reorder functions within ggplot weren't working all of a sudden.
  ggplot(., aes(x=FDR_p, y = GO_term)) + 
  geom_point(color = "#B538A2", size = 3) +
  #scale_size_continuous(breaks = c(3, 6, 9, 12, 20, 30)) +
  xlim(0,0.05) +
  labs(x="P-Value (FDR)", y = "GO Term") +
  theme(panel.background = element_blank(), 
        text = element_text(size = 13, colour = "black"), 
        axis.text.y = element_text(size = 10, colour = "black", margin = margin(r = -8)),
        axis.text.x = element_text(size = 10, colour = "black"),
        axis.ticks.y = element_blank())
ggsave("vSBN 5min Negative LFC all specific GO terms2.jpeg", height = 3, width = 6.15, units = "in")
ggsave("vSBN 5min Negative LFC all specific GO terms2.eps", height = 3, width = 6.15, units = "in")
ggsave("vSBN 5min Negative LFC all specific GO terms2.pdf", height = 3, width = 6.15, units = "in")



#### Gene expression plots for immune and synapse genes ####

# I pulled GE and agg score data for some genes that exemplify the enrichment observed in the vSBN, 5min agg deseq results
# These are DEGs from vSBN 5 min aggression

ge_agg_long <- read_xlsx("Data frames for figures.xlsx", sheet = 4)

# IL18 (tree_swallow018772t1_sg0.0), C3 (tree_swallow000637t4_sg0.3)

ge_agg_long %>% filter(GeneName == "C3" | GeneName == "IL18") %>%
  ggplot(., aes(x=Agg5min, y=NormGeneExpression)) +
  geom_smooth(method = "lm", color = "darkgrey")+
  geom_point(size = 2, color = "#F57340") +
  labs(x="5 Minute Aggression Score", y = "Counts") +
  xlim(0,60) +
  theme(panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor = element_blank(),
        text = element_text(size = 13, colour = "black"), 
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 12, colour = "black"),
        strip.text.x = element_text(size = 13),
        strip.background = element_blank()) +
  facet_wrap(~GeneName, ncol = 2, scales = "free")
ggsave("Immune Gene Expression plot_12Jul24.jpeg", height = 3, width = 4, units = "in")
ggsave("Immune Gene Expression plot_12Jul24.pdf", height = 3, width = 4, units = "in")

# Grin2A (tree_swallow000664t1_sg0.0), KCNH8 (tree_swallow001286t1_sg0.0)
ge_agg_long %>% filter(GeneName == "GRIN2A" | GeneName == "KCNH8") %>%
  ggplot(., aes(x=Agg5min, y=NormGeneExpression)) +
  geom_smooth(method = "lm", color = "darkgrey")+
  geom_point(size = 2, color = "#B538A2") +
  labs(x="5 Minute Aggression Score", y = "Counts") +
  xlim(0,60) +
  theme(panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor = element_blank(),
        text = element_text(size = 13, colour = "black"), 
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 12, colour = "black"),
        strip.text.x = element_text(size = 13),
        strip.background = element_blank()) +
  facet_wrap(~GeneName, ncol = 2, scales = "free")
ggsave("Synaptic plascticy Gene Expression plot_12Jul24.jpeg", height = 3, width = 4, units = "in")
ggsave("Synaptic plascticy Gene Expression plot_12Jul24.pdf", height = 3, width = 4, units = "in")


# Overlapping DEGs for aggression DEGs AND genes associated with treatment before p-val correction (Figure 5) --------
# A version that highlights the overlapping genes between aggression and treatement 

shared <- read_xlsx("Data frames for figures.xlsx", sheet = 5)

ggplot(shared, aes(x=LFC_agg.vSBN, y = LFC_treat.vSBN, color = Trait.x)) + 
  geom_point(size = 2) + 
  geom_abline(slope = 0, intercept = 0) + 
  geom_vline(xintercept = 0) + 
  ylim(-1,3) + xlim(-0.02,0.225) +
  scale_color_manual(breaks = c("5 min", "30 min"), values = c("midnightblue", "coral")) +
  labs(x="Log2 Fold Change for Aggression Score", 
       y = "Log2 Fold Change for Treatment",
       color = "Aggression Score") +
  theme(panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor = element_line(colour = "grey90"),
        text = element_text(size = 13, colour = "black"), 
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 12, colour = "black"))
ggsave("overlapping genes in vSBN_24Jun24.jpeg", height = 3.5, width = 5, units = "in")


# Log2 fold change for all aggression DEGs (supplemental) --------------------------------


all_agg_DEG <- read_xlsx("Data frames for figures.xlsx", sheet = 6)

# Transcripts that were DEGs for aggression generally, both brain regions. 
p2a <- ggplot(subset(all_agg_DEG, Trait == "Agg5min" & BrainRegion == "vSBN"), 
              aes(x= LFC_agg5.vSBN, y = LFC_treat.vSBN)) + 
  geom_point(alpha = 0.3) + 
  geom_abline(slope = 0, intercept = 0) + 
  ylim(-1.2,3) +
  xlim(-.1, .24) +
  geom_vline(xintercept = 0) + 
  geom_label(x = 0.08, y = 2, label = "149", label.size = 0, size = 4) +
  geom_label(x = -0.05, y = 2, label = "18", label.size = 0, size = 4) +
  geom_label(x = -0.05, y = -1.1, label = "177", label.size = 0, size = 4) +
  geom_label(x = 0.08, y = -1.1, label = "72", label.size = 0, size = 4) +
  labs(x="5 Minute Aggression Score Log2 Fold Change", y="Treatment Log2 Fold Change",
       title = "ventral SBN \nDEGs for 5 Minute Aggression Score") +
  theme_minimal()
p2b <- ggplot(subset(all_agg_DEG, Trait == "Agg30min" & BrainRegion == "vSBN"), 
              aes(x= LFC_agg30.vSBN, y = LFC_treat.vSBN)) + 
  geom_point(alpha = 0.3) + 
  geom_abline(slope = 0, intercept = 0) + 
  ylim(-1.2,3) +
  xlim(-.1, .24) +
  geom_vline(xintercept = 0) + 
  geom_label(x = 0.08, y = 2, label = "143", label.size = 0, size = 4) +
  geom_label(x = -0.05, y = 2, label = "6", label.size = 0, size = 4) +
  geom_label(x = -0.05, y = -1.1, label = "26", label.size = 0, size = 4) +
  geom_label(x = 0.08, y = -1.1, label = "71", label.size = 0, size = 4) +
  labs(x="30 Minute Aggression Score Log2 Fold Change", y="Treatment Log2 Fold Change",
       title = "ventral SBN \nDEGs for 30 Minute Aggression Score") +
  theme_minimal()

p2c <- ggplot(subset(all_agg_DEG, Trait == "Agg5min" & BrainRegion == "dSBN"), 
              aes(x= LFC_agg5.dSBN, y = LFC_treat.dSBN)) + 
  geom_point(alpha = 0.8) + 
  geom_abline(slope = 0, intercept = 0) + 
  ylim(-1.2,3) +
  xlim(-.1, .24) +
  geom_vline(xintercept = 0) + 
  geom_label(x = 0.08, y = 2, label = "13", label.size = 0, size = 4) +
  geom_label(x = -0.05, y = 2, label = "0", label.size = 0, size = 4) +
  geom_label(x = -0.05, y = -1.1, label = "0", label.size = 0, size = 4) +
  geom_label(x = 0.08, y = -1.1, label = "5", label.size = 0, size = 4) +
  labs(x="5 Minute Aggression Score Log2 Fold Change", y="Treatment Log2 Fold Change",
       title = "dorsal SBN \nDEGs for 5 Minute Aggression Score") +
  theme_minimal()
p2d <- ggplot(subset(all_agg_DEG, Trait == "Agg30min" & BrainRegion == "dSBN"), 
              aes(x= LFC_agg30.dSBN, y = LFC_treat.dSBN)) + 
  geom_point(alpha = 0.3) + 
  geom_abline(slope = 0, intercept = 0) + 
  ylim(-1.2,3) +
  xlim(-.1, .24) +
  geom_vline(xintercept = 0) + 
  labs(x="30 Minute Aggression Score Log2 Fold Change", y="Treatment Log2 Fold Change",
       title = "dorsal SBN \nNo DEGs for 30 Minute Aggression Score") +
  theme_minimal()
grid.arrange(p2a, p2b, p2c, p2d, nrow = 2)




