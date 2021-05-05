library(dplyr)
library(MPN)

#Plate 1: AC1
AC1 <- mpn(positive = c(7, 3, 2, 1), tubes = c(8, 8, 8, 8), amount = c(.01, .001, .0001, .00001))
AC1 <- as.data.frame(AC1)

#Plate 2: AC2
AC2 <- mpn(positive = c(8, 5, 2, 2), tubes = c(8, 8, 8, 8), amount = c(.01, .001, .0001, .00001))
AC2 <- as.data.frame(AC2)

#Plate 3: AC3
AC3 <- mpn(positive = c(6, 3, 3, 2), tubes = c(8, 8, 8, 8), amount = c(.01, .001, .0001, .00001))
AC3 <- as.data.frame(AC3)

#Plate 4: DC1
DC1 <- mpn(positive = c(5, 4, 4, 1), tubes = c(8, 8, 8, 8), amount = c(.01, .001, .0001, .00001))
DC1 <- as.data.frame(DC1)

#Plate 5: DC2
DC2 <- mpn(positive = c(7, 5, 5, 4), tubes = c(8, 8, 8, 8), amount = c(.01, .001, .0001, .00001))
DC2 <- as.data.frame(DC2)

#Plate 6: DC3
DC3 <- mpn(positive = c(7, 4, 1, 3), tubes = c(8, 8, 8, 8), amount = c(.01, .001, .0001, .00001))
DC3 <- as.data.frame(DC3)

#Plate 7: AS1
AS1 <- mpn(positive = c(6, 5, 2, 2), tubes = c(8, 8, 8, 8), amount = c(.01, .001, .0001, .00001))
AS1 <- as.data.frame(AS1)

#Plate 8: AS2
AS2 <- mpn(positive = c(3, 1, 2, 2), tubes = c(8, 8, 8, 8), amount = c(.01, .001, .0001, .00001))
AS2 <- as.data.frame(AS2)

#Plate 9: AS3
AS3 <- mpn(positive = c(4, 3, 2, 2), tubes = c(8, 8, 8, 8), amount = c(.01, .001, .0001, .00001))
AS3 <- as.data.frame(AS3)

#Plate 10: DS1
DS1 <- mpn(positive = c(5, 3, 3, 4), tubes = c(8, 8, 8, 8), amount = c(.01, .001, .0001, .00001))
DS1 <- as.data.frame(DS1)

#Plate 11: DS2
DS2 <- mpn(positive = c(5, 4, 3, 4), tubes = c(8, 8, 8, 8), amount = c(.01, .001, .0001, .00001))
DS2 <- as.data.frame(DS2)

#Plate 12: DS3
DS3 <- mpn(positive = c(6, 5, 4, 4), tubes = c(8, 8, 8, 8), amount = c(.01, .001, .0001, .00001))
DS3 <- as.data.frame(DS3)

#Plate 13: AL1
AL1 <- mpn(positive = c(0, 2, 0, 1), tubes = c(8, 8, 8, 8), amount = c(.01, .001, .0001, .00001))
AL1 <- as.data.frame(AL1)

#Plate 14: AL2
AL2 <- mpn(positive = c(4, 1, 2, 2), tubes = c(8, 8, 8, 8), amount = c(.01, .001, .0001, .00001))
AL2 <- as.data.frame(AL2)

#Plate 15: AL3
AL3 <- mpn(positive = c(3, 3, 3, 3), tubes = c(8, 8, 8, 8), amount = c(.01, .001, .0001, .00001))
AL3 <- as.data.frame(AL3)

#Plate 16: DL1
DL1 <- mpn(positive = c(2, 3, 3, 2), tubes = c(8, 8, 8, 8), amount = c(.01, .001, .0001, .00001))
DL1 <- as.data.frame(DL1)

#Plate 17: DL2
DL2 <- mpn(positive = c(4, 3, 2, 1), tubes = c(8, 8, 8, 8), amount = c(.01, .001, .0001, .00001))
DL2 <- as.data.frame(DL2)

#Plate 18: DL3
DL3 <- mpn(positive = c(0, 3, 2, 3), tubes = c(8, 8, 8, 8), amount = c(.01, .001, .0001, .00001))
DL3 <- as.data.frame(DL3)


plates_all <- bind_rows(AC1, AC2, AC3, DC1, DC2, DC3, AS1, AS2, AS3, DS1, DS2, DS3, AL1, AL2, AL3, DL1, DL2, DL3) %>% 
  mutate(Treatment = c("AC1", "AC2", "AC3", "DC1", "DC2", "DC3", "AS1", "AS2", "AS3", "DS1", "DS2", "DS3", "AL1", "AL2", "AL3", "DL1", "DL2", "DL3"), .before = 1) %>% 
  mutate(Sed_ad = c("AC", "AC", "AC", "DC", "DC", "DC", "AS", "AS", "AS", "DS", "DS", "DS", "AL", "AL", "AL", "DL", "DL", "DL"), .before = 1) %>% 
  mutate(Sed = c("C", "C", "C", "C", "C", "C", "S", "S", "S", "S", "S", "S", "L", "L", "L", "L", "L", "L"), .before = 1) %>% 
  mutate(Ad = c("A", "A", "A", "D", "D", "D", "A", "A", "A", "D", "D", "D", "A", "A", "A", "D", "D", "D"), .before = 1) %>% 
  mutate(Replicate = c("1", "2", "3", "1", "2", "3", "1", "2", "3", "1", "2", "3", "1", "2", "3", "1", "2", "3"), .before = 1)



plates_all_t <- plates_all %>%
  t() %>%
  as.data.frame(stringsAsFactors = F) %>%
  rownames_to_column("value") %>%
  `colnames<-`(.[1,]) %>%
  .[-1,] %>%
  `rownames<-`(NULL) %>% 
  rename(Values = Treatment)

a <- plates_all[-c(1),]
boxplot(MPN ~ Sed_ad, data = a)

sed.aov <- aov(MPN ~ Sed_ad, plates_all)
summary(sed.aov)
turkey <- HSD.test(sed.aov, trt = "Sed_ad")

ad.aov <- aov(MPN ~ Ad, a)
summary(ad.aov)

  
