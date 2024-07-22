# Testosterone analysis for incubation TRES transcriptome analysis
# When I ran this, the concentrations were imputed for anything that fell below Standard 8 (which affected 6 samples) - see printout of sample concentrations at the end.

library(dplyr)
library(ggplot2)
library(readxl)

traits <- read_xlsx("Supplement0_Traits.xlsx", sheet = 1)

### Testosterone t-test
# Sti v Con

# Plot first to see if we need to transform or use a non-parametric test
ggplot(traits, aes(x=testosterone, fill = factor(treatment_0con_1sti))) + geom_density(alpha = 0.5)
ggplot(traits, aes(x=log(testosterone), fill = factor(treatment_0con_1sti))) + geom_density(alpha = 0.5) # Logged looks way better.

# Separate the two conditions
sti <- traits %>% filter(treatment_0con_1sti == 1)
con <- traits %>% filter(treatment_0con_1sti == 0)

# Shapiro-Wilk test for normality - p-val < 0.05 indicates it's NOT normally distributed.
shapiro.test(traits$testosterone) # p = 1.888e-05
shapiro.test(log(traits$testosterone)) # p = 0.03918, much better, still not perfect.

shapiro.test(sti$testosterone) # p = 0.06 - normalish, but not great
shapiro.test(log(sti$testosterone)) # p = 0.26, that's better

shapiro.test(con$testosterone) # p = 0.0002...
shapiro.test(log(con$testosterone)) # 0.175, that's better.

# Running t-test on logged values
  # The var.equal says we don't assume variance is the same between the two populations. 
    # Results are identical either way (which prob means var are equal). Because this is false, it's a "welch two sample t-test".
t.test(log(con$testosterone), log(sti$testosterone), alternative = "two.sided", var.equal = FALSE) 
# t = 0.19554, df = 17.34, p-value = 0.8472

### Testosterone linear model
# Among STI birds, look at relationship between STI score and T concentration and Latency to capture and T concentration

# First with 5-minute aggression score
mod5min <- lm(log(testosterone) ~ agg5min, data = sti)
summary(mod5min)
#               Estimate    Std. Error    t value   Pr(>|t|)    
# (Intercept)   4.296043    0.585341      7.339     8.08e-05 ***
#   agg5min     -0.007221   0.020699      -0.349    0.736    

# Now with 30-minute aggression score
mod30min <- lm(log(testosterone) ~ agg30min, data = sti)
summary(mod30min)
#               Estimate    Std. Error    t value   Pr(>|t|)    
# (Intercept)   4.51819     0.75648       5.973     0.000333 ***
#   agg30min    -0.01049    0.01842       -0.569    0.584719    

latency <- lm(log(testosterone) ~ latency_sti_start_to_cap_min, data = sti)
summary(latency)
#                         Estimate    Std. Error    t value     Pr(>|t|)   
# (Intercept)             6.69897     1.78361       3.756       0.00558 **
#   sti_start_to_cap_min  -0.01684    0.01152       -1.462      0.18180 

### Sample concentrations - these are rounded a bit for some reason, but good enough to recreate if needed.

traits %>% select(id, testosterone) # Returns...
# id                testosterone
# 1 PATSF21         64.8
# 2 PATSF22         73  
# 3 PATSF23         69.7
# 4 PATSF24         63.3
# 5 PATSF25         82.4
# 6 PATSF26         60.4
# 7 PATSF27         73.6
# 8 PATSF28         22.5
# 9 PATSF29         22.5
# 10 PATSF30        22.5
# 11 PATSF31        22.5
# 12 PATSF32        441. 
# 13 PATSF33        81.6
# 14 PATSF34        22.5
# 15 PATSF35        22.5
# 16 PATSF36        221. 
# 17 PATSF37        188. 
# 18 PATSF38        60  
# 19 PATSF39        130. 
# 20 PATSF40        99.6