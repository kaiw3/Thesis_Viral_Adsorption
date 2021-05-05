# Data Reading and Tidying
library(tidyverse)
library(agricolae)

plates <- read.csv("raw_data.csv")

plates_t <- reshape(plates, 
                    direction = "long",
                    varying = list(names(plates)[3:9]),
                    v.names = "Lysed_Wells",
                    idvar = c("Day", "Plate"),
                    timevar = "Column",
                    times = 1:7)


plates_t <- separate(plates_t, Plate, c("Ad", "Treatment", "Replicate"), sep = ";")

plates_t$AdT <- paste(plates_t$Ad, plates_t$Treatment)

plates_t$AdTCol <- paste(plates_t$AdT, plates_t$Column)

plates_t <- plates_t %>% mutate(Lysed_Wells_perc = (Lysed_Wells/ 8) * 100)

plates_t$Day <- as.character(plates_t$Day)
plates_t$Column <- as.character(plates_t$Column)
plates_t$Lysed_Wells_perc <- as.numeric(plates_t$Lysed_Wells_perc)

str(plates_t)

plates_t <- filter(plates_t, Column != 7)

# Assumptions

hist(plates_t$Lysed_Wells)

# Create Theme

theme_lysebar <- function(){
  theme_bw() +
    theme(text = element_text(family = "Helvetica Light"),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 16),
          axis.line.x = element_line(color="black"),
          axis.line.y = element_line(color="black"),
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),
          plot.margin = unit(c(1, 3, 1, 1), units = , "cm"),
          plot.title = element_text(size = 18, vjust = 1, hjust = 0),
          legend.text = element_text(size = 12),
          legend.position = c(1.05, 0.2),
          legend.key = element_blank(),
          legend.background = element_rect(color = "black",
                                           fill = "transparent",
                                           size = 2, linetype = "blank"))
  
}


# Adsorped vs Desorped Sediment Effect Day 14 Pure Lysate

all.aov_data_14 <- filter(plates_t, Day == "14", Column == "2")
all.aov_14 <- aov(Lysed_Wells_perc ~ Treatment * Ad, all.aov_data_14)
summary(all.aov_14)
turkey_all_14 <- HSD.test(all.aov_14, trt = c("Treatment", "Ad"))


aov_data_summ <- all.aov_data_14 %>% 
  group_by(AdT) %>% 
  summarise(mean_trt = mean(Lysed_Wells_perc), sd_trt = sd(Lysed_Wells_perc), n_trt = n(), se_trt = sd(Lysed_Wells_perc)/sqrt(n())) %>% 
  ungroup() %>% 
  mutate(Ad = c("A", "A", "A", "D", "D", "D"), Sed = c("C", "L", "S", "C", "L", "S"))

ggplot(aov_data_summ, aes(fill = Ad, x = reorder(Sed, -mean_trt), y = mean_trt)) +
  geom_bar(stat = "summary", position = "dodge", alpha = 0.8) +
  geom_errorbar(aes(ymin = mean_trt - se_trt, ymax = mean_trt + se_trt), position = position_dodge(0.9), width = 0.3) +
  labs( x = "\nSediment Treatment", y = "Average Wells Lysed (%)", fill = "Desorption Treatment") +
  scale_x_discrete(labels=c("C" = "Control", "S" = "Small", "L" = "Large")) +
  scale_fill_manual(labels = c("A" = "Untreated", "D" = "Treated"), values = c("#2bcfcf", "#87bc4d")) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  theme_lysebar()

plot(all.aov_14, 1)
leveneTest(Lysed_Wells ~ Treatment*Ad, all.aov_data_14)

plot(all.aov_14, 2)
shapiro.test(all.aov_14$coefficients)


# Adsorped vs Desorped Sediment Effect Day 17 Pure Lysate

all.aov_data_17 <- filter(plates_t, Day == "17", Column == "2")
all.aov_17 <- aov(Lysed_Wells ~ Treatment * Ad, all.aov_data_17)
summary(all.aov_17)
turkey_all_17 <- HSD.test(all.aov_17, trt = c("Treatment", "Ad"))

(Ad_17_bar <- ggplot(all.aov_data_17, aes(fill = Ad, x = Treatment, y = Lysed_Wells)) +
    geom_bar(stat = "summary", position = "dodge"))


# MPN Analysis

mpn_plates <- read.csv("MPN_Data.csv")

mpn_plates <- separate(mpn_plates, Plate, c("Ad", "Treatment", "Replicate"), sep = ";")
mpn_plates$AdT <- paste(mpn_plates$Ad, mpn_plates$Treatment)

mpn_sed.aov <- lm(MPN ~ Treatment, mpn_plates)
summary(mpn_sed.aov)                             # Not significant

(mpn_bar <- ggplot(mpn_plates, aes(x = reorder(Treatment, -MPN), y = MPN)) +
    geom_bar(stat = "summary", position = "dodge"))


# Stepwise Linear Regression of Dilutions ----

#Get Data

plates_dils <- read.csv("raw_data_2.csv")

plates_dils_t <- reshape(plates_dils, 
                         direction = "long",
                         varying = list(names(plates_dils)[3:9]),
                         v.names = "Lysed_Wells",
                         idvar = c("Day", "Plate"),
                         timevar = "Dilution",
                         times = c("Blank", "0.1", "0.01", "0.001", "0.0001", "0.00001", "0.000001"))


plates_dils_t <- separate(plates_dils_t, Plate, c("Ad", "Treatment", "Replicate"), sep = ";")

plates_dils_t$AdT <- paste(plates_dils_t$Ad, plates_dils_t$Treatment)

plates_dils_t$AdTCol <- paste(plates_dils_t$AdT, plates_dils_t$Dilution)

plates_dils_t <- plates_dils_t %>% mutate(Lysed_Wells_perc = (Lysed_Wells/ 8) * 100)

plates_dils_t$Day <- as.character(plates_dils_t$Day)
plates_dils_t$Dilution <- as.character(plates_dils_t$Dilution)
plates_dils_t$Lysed_Wells_perc <- as.numeric(plates_dils_t$Lysed_Wells_perc)
str(plates_dils_t)

plates_dils_t <- plates_dils_t %>% filter(Dilution != "0.000001")

plates_dils_t_14 <- plates_dils_t %>% filter(Day == 14)
plates_dils_t_17 <- plates_dils_t %>% filter(Day == 17)

# Theme

theme_lyseline <- function(){
  theme_bw() +
    theme(text = element_text(family = "Helvetica Light"),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 16),
          axis.line.x = element_line(color="black"),
          axis.line.y = element_line(color="black"),
          panel.border = element_blank(),
          plot.margin = unit(c(1, 3.5, 1, 1), units = , "cm"),
          plot.title = element_text(size = 18, vjust = 1, hjust = 0),
          legend.text = element_text(size = 12),
          legend.position = c(1.05, 0.2),
          legend.key = element_blank(),
          legend.background = element_rect(color = "black",
                                           fill = "transparent",
                                           size = 2, linetype = "blank"))
  
}

# Graph Dilutions

ggplot(plates_dils_t, aes(x = reorder(Dilution, desc(Dilution)), y = Lysed_Wells_perc, group = interaction(Treatment,Ad))) +
  geom_line(aes(color = Treatment, linetype = Ad), stat = "summary", size = 1, alpha = 0.9) +
  labs( x = "\nDilution Factor", y = "Average Wells Lysed (%)", color = "Sediment Treatment", linetype = "Desorption Treatment") +
  scale_linetype_discrete(labels = c("A" = "Untreated", "D" = "Treated")) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  scale_color_manual(labels = c("C" = "Control", "L" = "Large", "S" = "Small"), values = c("#2bcfcf", "#87bc4d", "#ff781f")) +
  theme_lyseline() +
  facet_wrap(vars(Day), labeller = labeller(Day = c("10" = "Day 10", "14" = "Day 14", "17" = "Day 17")))


# Regression Lines for Dilutions

plates_dils_t_14_reg <- plates_dils_t_14 %>% filter(Dilution != "Blank")

ggplot(plates_dils_t_14_reg, aes(x = reorder(Dilution, desc(Dilution)), y = Lysed_Wells_perc, group = interaction(Treatment,Ad))) +
  geom_smooth(method = "lm", aes(color = Treatment, linetype = Ad), alpha = 0.2) +
  labs( x = "Dilution Factor", y = "Average Wells Lysed (%)", color = "Sediment Treatment", linetype = "Desorption Treatment") +
  scale_linetype_discrete(labels = c("A" = "Untreated", "D" = "Treated")) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  ylim(c(0, 100)) +
  scale_color_manual(labels = c("C" = "Control", "L" = "Large", "S" = "Small"), values = c("#2bcfcf", "#87bc4d", "#ff781f")) +
  theme_lyseline() 

# Linear Models for Gradients - Day 14

Tr_AC <- subset(plates_dils_t_14, AdT=="A C")
Tr_AC.lm <- aov(Lysed_Wells ~ Dilution, Tr_AC)
summary(Tr_AC.lm)

Tr_DC <- subset(plates_dils_t_14, AdT=="D C")
Tr_DC.lm <- aov(Lysed_Wells ~ Dilution, Tr_DC)
summary(Tr_DC.lm)

Tr_AS <- subset(plates_dils_t_14, AdT=="A S")     # NOT SIGNIFICANT
Tr_AS.lm <- aov(Lysed_Wells ~ Dilution, Tr_AS)
summary(Tr_AS.lm)

Tr_DS <- subset(plates_dils_t_14, AdT=="D S")
Tr_DS.lm <- aov(Lysed_Wells ~ Dilution, Tr_DS)
summary(Tr_DS.lm)

Tr_AL <- subset(plates_dils_t_14, AdT=="A L")     # NOT SIGNIFICANT
Tr_AL.lm <- aov(Lysed_Wells ~ Dilution, Tr_AL)
summary(Tr_AL.lm)

Tr_DL <- subset(plates_dils_t_14, AdT=="D L")     # NOT SIGNIFICANT
Tr_DL.lm <- aov(Lysed_Wells ~ Dilution, Tr_DL)
summary(Tr_DL.lm)


# Linear Models for Gradients - Day 17

Tr_AC_17 <- subset(plates_dils_t_17, AdT=="A C")
Tr_AC_17.lm <- aov(Lysed_Wells ~ Dilution, Tr_AC_17)
summary(Tr_AC_17.lm)

Tr_DC_17 <- subset(plates_dils_t_17, AdT=="D C")
Tr_DC_17.lm <- aov(Lysed_Wells ~ Dilution, Tr_DC_17)
summary(Tr_DC_17.lm)

Tr_AS_17 <- subset(plates_dils_t_17, AdT=="A S")    
Tr_AS_17.lm <- aov(Lysed_Wells ~ Dilution, Tr_AS_17)
summary(Tr_AS_17.lm)

Tr_DS_17 <- subset(plates_dils_t_17, AdT=="D S")
Tr_DS_17.lm <- aov(Lysed_Wells ~ Dilution, Tr_DS_17)
summary(Tr_DS_17.lm)

Tr_AL_17 <- subset(plates_dils_t_17, AdT=="A L")     # NOT SIGNIFICANT
Tr_AL_17.lm <- aov(Lysed_Wells ~ Dilution, Tr_AL_17)
summary(Tr_AL_17.lm)

Tr_DL_17 <- subset(plates_dils_t_17, AdT=="D L")
Tr_DL_17.lm <- aov(Lysed_Wells ~ Dilution, Tr_DL_17)
summary(Tr_DL_17.lm)



