#flow extremes shift the balance of stochastic and deterministic processes across aquatic biomes

#code written by Jill Tupitza

#open libs
library(vegan)
library(dplyr)
library(ggplot2)
library(BiodiversityR)
library(rpart)
library(patchwork)

####################################
#Barataria Bay, LA 

#read in data

all_comm <-read.csv("~/Documents/Dissertation/Data/community_nozero.csv") #species comm + env, zeros removed
comm_zero <-read.csv("~/Documents/Dissertation/Data/community.csv") #species comm + env
all_fam_nozero <-read.csv("~/Documents/Dissertation/Data/all_families_nozero.csv") #family +feeding guild + env, zeros removed
all_fam <-read.csv("~/Documents/Dissertation/Data/all_seasons_families.csv") #family + feeding guild + env

#make useful subsets
#all with site 4 and zero abundance sites intact
sp.community <-comm_zero[1:249, 5:113] #subset community
sp.community[is.na(sp.community)] <- 0 #replace NA with 0 in df
sp.env <-data.frame(comm_zero$Season, comm_zero$Salinity, comm_zero$Salinity_class, comm_zero$Temperature,
                    comm_zero$pH, comm_zero$percent_clay, comm_zero$percent_silt, comm_zero$percent_sand,
                    comm_zero$percent_LOI, comm_zero$SAL_anom)
colnames(sp.env) <-c("season", "salinity", "salinity_class", "temperature", "pH", "percent_clay",
                     "percent_silt", "percent_sand", "percent_LOI", "sal_anom")
sp.env <-sp.env[1:249,] #environmental dataset without pesky blank last row

#family subsets
fam.community <-all_fam[1:249,5:69] #with hippogriffs - zeros abundance rows intact
fam.community[is.na(fam.community)] <-0 #replace NAs with 0 in df

fam.community.nozero <-all_fam_nozero[1:195,5:68] #without hippogriffs -zero abundance rows deleted
fam.community.nozero[is.na(fam.community.nozero)] <-0 #replace NAs with 0 in df


#feeding guild subsets
fg.community <-all_fam[1:249,82:94] #zero abundance rows intact
fg.community[is.na(fg.community)] <- 0 #replace NAs with 0 in df

fg.community.nozero <-all_fam_nozero[1:195,82:93] # zero abundance rows deleted
fg.community.nozero[is.na(fg.community.nozero)] <-0 #replace NAs with 0s in df

######################################################
#PERMANOVAs- Table S1 - electronic supplementary

#########################################
#PERMANOVAs
#all with site 4 intact
sp.community <-comm_zero[1:249, 5:113] #subset community
sp.community[is.na(sp.community)] <- 0 #replace NA with 0 in df
bray = vegdist(sp.community, method = "bray")
sp.env <-data.frame(comm_zero$Season, comm_zero$Salinity, comm_zero$Salinity_class, comm_zero$Temperature,
                    comm_zero$pH, comm_zero$percent_clay, comm_zero$percent_silt, comm_zero$percent_sand,
                    comm_zero$percent_LOI, comm_zero$SAL_anom)
colnames(sp.env) <-c("season", "salinity", "salinity_class", "temperature", "pH", "percent_clay",
                     "percent_silt", "percent_sand", "percent_LOI", "sal_anom")
sp.env <-sp.env[1:249,] #environmental dataset

permanova_total=adonis2(sp.community ~ percent_LOI*percent_clay*percent_sand*sal_anom*salinity_class*season+temperature, data = sp.env, permutations =999, method ="bray")
permanova_total #every interaction in the book is happening. Parse out by season and salinity class
summary(permanova_total)


#with salinity anomaly 
permanova_sal =adonis2(sp.community.nozero ~ percent_LOI*percent_clay*percent_sand*salinity_class*season*sal_anomaly+temperature, data=sp.env.nozero, permutations=999, method="bray")

#PERMANOVA Spring- overall, interactive effects between percent clay: percent OM in all seasons and in both meso/polyhaline regimes
#ALL : interactive effects between salinity class:percent clay, salinity class:percent LOI, percent LOI:percent clay, 
#and a three-way interaction between salinity class:LOI:clay
#sal anom is significant interaction
spring_comm=comm_zero[comm_zero$Season=="Spring", 5:114]
spring_env=sp.env[sp.env$season == "Spring",]
permanova_spring=adonis2(spring_comm ~salinity_class*percent_clay*percent_LOI*sal_anom+temperature, data=spring_env, permutations=999, method="bray")
permanova_spring
summary(permanova_spring)

#Df SumOfSqs       R2       F Pr(>F)    
#salinity_class                           2  2.94676 0.118264 6.77174  0.001 ***
#  percent_clay                             1  0.72668 0.029164 3.33987  0.001 ***
# percent_LOI                              1  0.43956 0.017641 2.02024  0.017 *  
#  temperature                              1  0.57600 0.023117 2.64732  0.002 ** 
#  salinity_class:percent_clay              2  0.64725 0.025977 1.48741  0.082 .  
#salinity_class:percent_LOI               2  1.13591 0.045588 2.61035  0.001 ***
#  percent_clay:percent_LOI                 1  0.92011 0.036927 4.22889  0.001 ***
#  salinity_class:percent_clay:percent_LOI  1  0.55352 0.022215 2.54404  0.008 ** 
#  Residual                                78 16.97104 0.681108                   
#Total                                   89 24.91683 1.000000                  

#PERMANOVA - Spring + polyhaline regime
#interactive effects between percent clay: percent LOI : sal anom and all with each other
#sig main effect of temp
spring_poly_comm=comm_zero[comm_zero$Season=="Spring" & comm_zero$Salinity_class=="polyhaline",5:114]
spring_poly_env=spring_env[spring_env$salinity_class=="polyhaline",]
permanova_spring_poly=adonis2(spring_poly_comm ~percent_clay*percent_LOI*sal_anom+temperature, data=spring_poly_env, permutations=999, method="bray")
permanova_spring_poly
summary(permanova_spring_poly)

#Df SumOfSqs      R2      F Pr(>F)    
#  percent_clay                       1   0.6751 0.05163 3.0738  0.001 ***
#  percent_LOI                        1   0.9641 0.07374 4.3898  0.001 ***
#  sal_anom                           1   0.8292 0.06342 3.7757  0.001 ***
#  temperature                        1   0.6259 0.04787 2.8500  0.003 ** 
#  percent_clay:percent_LOI           1   0.8939 0.06837 4.0700  0.001 ***
#  percent_clay:sal_anom              1   0.6641 0.05079 3.0239  0.001 ***
#  percent_LOI:sal_anom               1   0.5309 0.04061 2.4174  0.003 ** 
#  percent_clay:percent_LOI:sal_anom  1   0.6440 0.04926 2.9323  0.002 ** 
#  Residual                          33   7.2478 0.55432                  
#Total                             41  13.0752 1.00000    

#PERMANOVA - Spring + mesohaline regime- interactive effect between percent clay:percent LOI (p=0.003)
spring_meso_comm=comm_zero[comm_zero$Season=="Spring" & comm_zero$Salinity_class=="mesohaline",5:114]
spring_meso_env=spring_env[spring_env$salinity_class=="mesohaline",]
permanova_spring_meso=adonis2(spring_meso_comm ~percent_clay*percent_LOI*sal_anom+temperature, data=spring_meso_env, permutations=999, method="bray")
permanova_spring_meso
summary(permanova_spring_meso)

#Df SumOfSqs      R2      F Pr(>F)    
#percent_clay                       1   0.2245 0.03486 1.6746  0.134    
#percent_LOI                        1   0.4023 0.06247 3.0006  0.006 ** 
#  sal_anom                           1   0.3747 0.05818 2.7948  0.010 ** 
#  temperature                        1   0.6200 0.09627 4.6242  0.001 ***
#  percent_clay:percent_LOI           1   0.8279 0.12856 6.1755  0.001 ***
#  percent_clay:sal_anom              1   0.2662 0.04133 1.9855  0.070 .  
#percent_LOI:sal_anom               1   0.2171 0.03372 1.6197  0.132    
#percent_clay:percent_LOI:sal_anom  1   0.2895 0.04495 2.1592  0.052 .  
#Residual                          24   3.2176 0.49964                  
#Total                             32   6.4397 1.00000 

#PERMANOVA - spring + oligohaline regime -significant main effects of % clay (p=0.005) and sal anom (p=0.005)
spring_oligo_comm=comm_zero[comm_zero$Season=="Spring" & comm_zero$Salinity_class=="oligohaline",5:114]
spring_oligo_env=spring_env[spring_env$salinity_class=="oligohaline",]
permanova_spring_oligo=adonis2(spring_oligo_comm ~percent_clay*percent_LOI*sal_anom+temperature, data=spring_oligo_env, permutations=999, method="bray")
permanova_spring_oligo
summary(permanova_spring_oligo)
#Df SumOfSqs      R2       F Pr(>F)   
#percent_clay  1  0.59601 0.24276  8.9492  0.005 **
#  percent_LOI   1  0.10509 0.04280  1.5779  0.138   
#sal_anom      1  0.99709 0.40612 14.9714  0.005 **
#  temperature   1  0.09098 0.03706  1.3661  0.230   
#Residual     10  0.66600 0.27126                  
#Total        14  2.45517 1.00000               


#PERMANOVA Summer -multiple interactions, parsing out by salinity class should straighten this out a bit
summer_comm=comm_zero[comm_zero$Season=="Summer", 5:101]
summer_env=sp.env[sp.env$season == "Summer",]
permanova_summer=adonis2(summer_comm ~salinity_class*percent_clay*percent_LOI*sal_anom+temperature, data=summer_env, permutations=999, method="bray")
permanova_summer
summary(permanova_summer)
summer_env_nozero =sp.env.nozero[sp.env.nozero$season=="Summer",]

#Df SumOfSqs      R2      F Pr(>F)    
#salinity_class                           2   2.6422 0.10550 5.6122  0.001 ***
#  percent_clay                             1   0.6475 0.02585 2.7505  0.003 ** 
#  percent_LOI                              1   0.8087 0.03229 3.4354  0.001 ***
#  temperature                              1   0.2878 0.01149 1.2228  0.231    
#salinity_class:percent_clay              1   0.5450 0.02176 2.3152  0.002 ** 
#  salinity_class:percent_LOI               1   0.5258 0.02099 2.2335  0.006 ** 
#  percent_clay:percent_LOI                 1   0.2450 0.00978 1.0409  0.374    
#salinity_class:percent_clay:percent_LOI  1   0.5113 0.02042 2.1721  0.018 *  
#  Residual                                80  18.8317 0.75192                  
#Total                                   89  25.0449 1.00000  

#PERMANOVA - Summer + polyhaline regime: only one site. omit from analysis
summer_poly_comm=comm_zero[comm_zero$Season=="Summer" & comm_zero$Salinity_class=="polyhaline",5:101]
summer_poly_env=summer_env[summer_env$salinity_class=="polyhaline",]
permanova_summer_poly=adonis2(summer_poly_comm ~percent_clay*percent_LOI*sal_anom+temperature, data=summer_poly_env, permutations=999, method="bray")
permanova_summer_poly

#PERMANOVA - Summer + mesohaline regime: interactive effect between percent clay: percent LOI and sal anomaly
summer_meso_comm=all_comm[all_comm$Season=="Summer" & all_comm$Salinity_class=="mesohaline",5:101]
summer_meso_env=summer_env_nozero[summer_env_nozero$salinity_class=="mesohaline",]
permanova_summer_meso=adonis2(summer_meso_comm ~percent_clay*percent_LOI*sal_anomaly+temperature, data=summer_meso_env, permutations=999, method="bray")
permanova_summer_meso
#Df SumOfSqs      R2      F Pr(>F)    
#percent_clay                          1   0.9351 0.04340 2.4531  0.001 ***
#  percent_LOI                           1   0.7973 0.03700 2.0916  0.006 ** 
#  sal_anomaly                           1   0.7474 0.03469 1.9607  0.008 ** 
#  temperature                           1   0.6145 0.02852 1.6121  0.047 *  
#  percent_clay:percent_LOI              1   0.7319 0.03397 1.9201  0.007 ** 
#  percent_clay:sal_anomaly              1   0.6939 0.03220 1.8203  0.012 *  
#  percent_LOI:sal_anomaly               1   0.9696 0.04500 2.5436  0.001 ***
#  percent_clay:percent_LOI:sal_anomaly  1   0.8099 0.03759 2.1246  0.003 ** 
# Residual                             40  15.2479 0.70764                  
#Total                                48  21.5476 1.00000  

#PERMANOVA - Summer + oligohaline regime - interactive effects between everything, but no 3-way interaction
summer_oligo_comm=all_comm[all_comm$Season=="Summer" & all_comm$Salinity_class=="oligohaline",5:113]
summer_oligo_env=summer_env_nozero[summer_env_nozero$salinity_class=="oligohaline",]
permanova_summer_oligo=adonis2(summer_oligo_comm ~percent_clay*percent_LOI*sal_anomaly+temperature, data=summer_oligo_env, permutations=999, method="bray")
permanova_summer_oligo
#Df SumOfSqs      R2      F Pr(>F)    
#percent_clay                          1   0.6212 0.07281 2.4709  0.008 ** 
#  percent_LOI                           1   1.2743 0.14936 5.0689  0.001 ***
#  sal_anomaly                           1   0.5812 0.06812 2.3117  0.009 ** 
#  temperature                           1   0.8886 0.10414 3.5344  0.002 ** 
#  percent_clay:percent_LOI              1   0.5530 0.06482 2.1997  0.003 ** 
#  percent_clay:sal_anomaly              1   0.7380 0.08650 2.9355  0.002 ** 
#  percent_LOI:sal_anomaly               1   0.5069 0.05941 2.0163  0.021 *  
#  percent_clay:percent_LOI:sal_anomaly  1   0.3521 0.04127 1.4006  0.156    
#Residual                             12   3.0169 0.35359                  
#Total                                20   8.5322 1.00000  

#PERMANOVA Fall- no oligohaline sites, interactions between salinity class, % clay and %LOI
fall_comm=all_comm[all_comm$Season=="Fall", 5:113]
fall_comm[is.na(fall_comm)] <- 0 #replace NAs with zero in df
fall_env=sp.env.nozero[sp.env.nozero$season=="Fall",]
permanova_fall=adonis2(fall_comm ~salinity_class*percent_clay*percent_LOI*sal_anomaly+temperature, data=fall_env, permutations=999, method="bray")
permanova_fall

#Df SumOfSqs      R2      F Pr(>F)   
#salinity_class                           2   0.4270 0.04903 1.8745  0.009 **
#percent_clay                             1   0.1527 0.01754 1.3409  0.192   
#percent_LOI                              1   0.2620 0.03008 2.3002  0.019 * 
#temperature                              1   0.1816 0.02086 1.5949  0.092 . 
#salinity_class:percent_clay              2   0.5121 0.05881 2.2481  0.002 **
#salinity_class:percent_LOI               2   0.3096 0.03555 1.3591  0.136   
#percent_clay:percent_LOI                 1   0.1150 0.01321 1.0101  0.434   
#salinity_class:percent_clay:percent_LOI  2   0.3701 0.04250 1.6247  0.045 * 
#Residual                                56   6.3779 0.73242                 
#Total                                   68   8.7079 1.00000  

#PERMANOVA Fall + polyhaline regime: percent LOI is significant main effect (p=0.048)
fall_poly_comm=all_comm[all_comm$Season=="Fall" & all_comm$Salinity_class=="polyhaline",5:113]
fall_poly_env=fall_env[fall_env$salinity_class=="polyhaline",]
permanova_fall_poly=adonis2(fall_poly_comm ~percent_clay*percent_LOI*sal_anomaly+temperature, data=fall_poly_env, permutations=999, method="bray")
permanova_fall_poly

#Df SumOfSqs      R2      F Pr(>F)  
#percent_clay  1   0.1854 0.03988 0.4240  0.951  
#percent_LOI   1   0.8222 0.17689 1.8809  0.048 *
#  sal_anomaly   1   0.4146 0.08919 0.9484  0.496  
#temperature   1   0.1659 0.03569 0.3795  0.959  
#Residual      7   3.0600 0.65834                
#Total        11   4.6481 1.00000                
---
  
  #PERMANOVA Fall + mesohaline regime - significant interactive effect between percent LOI and sal anomaly
  fall_meso_comm=all_comm[all_comm$Season=="Fall" & all_comm$Salinity_class=="mesohaline",5:113]
fall_meso_comm[is.na(fall_meso_comm)] <- 0 
fall_meso_env=fall_env[fall_env$salinity_class=="mesohaline",]
permanova_fall_meso=adonis2(fall_meso_comm ~percent_clay*percent_LOI*sal_anomaly+temperature, data=fall_meso_env, permutations=999, method="bray")
permanova_fall_meso
#Df SumOfSqs      R2      F Pr(>F)   
#percent_clay                          1   0.7450 0.10177 1.9155  0.033 * 
#  percent_LOI                           1   0.6320 0.08633 1.6250  0.051 . 
#sal_anomaly                           1   0.3516 0.04803 0.9040  0.557   
#temperature                           1   0.5971 0.08156 1.5352  0.080 . 
#percent_clay:percent_LOI              1   0.3942 0.05384 1.0134  0.450   
#percent_clay:sal_anomaly              1   0.5454 0.07450 1.4022  0.097 . 
#percent_LOI:sal_anomaly               1   0.7366 0.10062 1.8939  0.009 **
#  percent_clay:percent_LOI:sal_anomaly  1   0.2075 0.02835 0.5336  0.915   
#Residual                              8   3.1116 0.42501                 
#Total                                16   7.3212 1.00000  

######################################################
#Canonical Correspondence Analysis (CCA) 
########################################
#Canonical Correspondence Analysis
#subset of driving env variables
sp.env.sub <-data.frame(sp.env.nozero$season, sp.env.nozero$salinity_class, sp.env.nozero$salinity, sp.env.nozero$sal_anomaly, sp.env.nozero$temperature, sp.env.nozero$percent_clay, sp.env.nozero$percent_LOI)
colnames(sp.env.sub) <-c("season", "salinity_class", "salinity", "salinity_anomaly", "temp", "percent_clay", "percent_LOI")


##cca with all seasons
################################################
sp.community.cca <-cca(sp.community.nozero~., sp.env.sub) #formula interface with all seasons
summary(sp.community.cca) #percentage of inertia by constrained (environmental) variables -9.92%
anova.cca(sp.community.cca) # test significant of model. Model is significant (p<0.001)
anova.cca(sp.community.cca, by ="terms") #test significance of env. variables: %LOI, %clay, salinity, salinity class, season. big residual. NOT temp
anova.cca(sp.community.cca, by="axis") #test significance of axes: CCA1, CCA2 both p<0.001.

#envfit- for env vectors
fit_all <-envfit(sp.community.nozero, sp.env.sub, permutations = 999, display="sites") #fit those env variables that went into the model

# correlation between axis 1 and organic matter (OM)
cca.score <-scores(sp.community.cca)
axes <-cca.score$sites
cor(axes[,1], sp.env.sub$percent_LOI)
#cor is -91.8
cor(axes[,1], sp.env.sub$salinity)
#cor is -46.8
cor(axes[,1], sp.env.sub$salinity_anomaly)

cor(axes[,2], sp.env.sub$percent_clay)
#cor is -0.41
cor(axes[,2], sp.env.sub$salinity)
#cor is -0.13
cor(axes[,2], sp.env.sub$salinity_anomaly)
#cor is 0.274
cor(axes[,2], sp.env.sub$temp)
#cor is 0.07
cor(axes[,2], sp.env.sub$percent_LOI)

#CCA w/o stochastic sites
determ <-all_comm[all_comm$Site==2 | all_comm$Site==5 | all_comm$Site ==6 | all_comm$Site==7 | all_comm$Site==11 | all_comm$Site ==12 | all_comm$Site==13 | all_comm$Site==14,]
determ_comm <-determ[1:98,5:113]
determ_comm[is.na(determ_comm)] <- 0  #deterministic community

#env of only deterministic sites
determ_env <-data.frame(determ$Season, determ$Salinity, determ$percent_clay, determ$percent_LOI, determ$Salinity_class, determ$Temperature)
colnames(determ_env) <- c("season", "salinity", "percent_clay", "percent_LOI", "salinity_class", "temp")

#CCA with deterministic sites
determ.cca <-cca(determ_comm~., determ_env)
summary(determ.cca) #only 14.25% explained by constrained vars

######### ####################
#Figure 4
#CCA species plot -  by salinity class ################################ Figure 4A
par(mfrow=c(2,2))
sp.env.sub$salinity_class <-as.factor(sp.env.sub$salinity_class)
with(sp.env.sub, levels(salinity_class))
colsalinity <- c("mediumvioletred","gray", "gray13") #color vector
scl <- 3 #set up scaling
plot(sp.community.cca, type = "n", scaling = scl)
with(sp.env.sub, points(sp.community.cca, display = "sites", col = colsalinity[salinity_class],
                        scaling = scl, pch = 21, bg = colsalinity[salinity_class],))
with(sp.env.sub, legend("topleft", legend = levels(salinity_class), bty = "n",
                        col = colsalinity, pch = 21, pt.bg = colsalinity, cex=1.2))
title(cex.lab=3)
scores_all <-scores(sp.community.cca, display = "sites")
fit_all <-envfit(scores_all, sp.env.sub[,4:6], permutations =999, display = "sites")
plot(fit_all, cex =1.0, col = "navy")
text(x=1.5, y=2.5, "A", cex=2)


##CCA - family plot ################# by salinity class - Figure 4B
par(mar = c(2, 2, 2, 2))

sp.env.sub$salinity_class <-as.factor(sp.env.sub$salinity_class)
with(sp.env.sub, levels(salinity_class))
colsalinity <- c("mediumvioletred", "gray", "gray13") #color vector
scl <- 3 #set up scaling
plot(fam.community.cca, type = "n", scaling = scl)
with(sp.env.sub, points(fam.community.cca, display = "sites", col = colsalinity[salinity_class],
                        scaling = scl, pch = 21, bg = colsalinity[salinity_class],))
with(sp.env.sub, legend("topleft", legend = levels(salinity_class), bty = "n",
                        col = colsalinity, pch = 21, pt.bg = colsalinity))
fam_scores_all <-scores(fam.community.cca, display = "sites")
fam_fit_all <-envfit(fam_scores_all, sp.env.sub[,3:6], permutations =999, display = "sites")
plot(fam_fit_all, cex =1, col = "navy")
text(x=1.5,y=1.75, "B", cex=2)

#feeding guilds by salinity class ####### - Figure 4C

sp.env.sub$salinity_class <-as.factor(sp.env.sub$salinity_class)
with(sp.env.sub, levels(salinity_class))
colsalinity <- c("mediumvioletred", "gray", "gray13") #color vector
scl <- 3 #set up scaling
plot(fg.community.cca, type = "n", scaling = scl)
with(sp.env.sub, points(fg.community.cca, display = "sites", col = colsalinity[salinity_class],
                        scaling = scl, pch = 21, bg = colsalinity[salinity_class],))
with(sp.env.sub, legend("topleft", legend = levels(salinity_class), bty = "n",
                        col = colsalinity, pch = 21, pt.bg = colsalinity))
fg_scores_all <-scores(fg.community.cca, display = "sites")
fg_fit_all <-envfit(fg_scores_all, sp.env.sub[,4:6], permutations =999, display = "sites")
plot(fg_fit_all, cex =1, col = "navy")
text(x=1.5,y=2.5, "C", cex=2)


#M####################################################
#Multivariate Regression Trees (Figure 5-6)
#De'ath 2002
#Multivariate regression trees with all sites and all seasons
#Step 1: subset community & env data - all zero abundance rows removed
devtools::install_github("cran/mvpart") #mvpart no longer active on CRAN- download from archive here
library(mvpart)
library(rpart)
all_comm <-read.csv("~/Documents/Dissertation/Data/community_nozero.csv")
everybody <-all_comm[1:195, 5:113] #community dataset
everybody[is.na(everybody)] <- 0 ##replace NAs with 0 in data frame

fam.community.nozero #family community dataset

fg.community.nozero #feeding guild community dataset

#env dataset
env_everybody <-data.frame(all_comm$percent_LOI, all_comm$percent_clay, all_comm$Season,
                           all_comm$percent_sand, all_comm$Salinity_class, all_comm$Salinity,
                           all_comm$pH, all_comm$percent_sand, all_comm$percent_silt, all_comm$SAL_anomaly)
#rename columns
env_everybody<-env_everybody[1:195,]
colnames(env_everybody) <-c("percent_LOI", "percent_clay", "season", "percent_sand", 
                            "salinity_class", "salinity", "pH", "percent_sand", "percent_silt", "salinity_anomaly")

#Step 2: perform Hellinger transformation on community data
everybody_hel <-decostand(everybody, "hellinger")

family_hel <-decostand(fam.community.nozero, "hellinger")

fg_hel <-decostand(fg.community.nozero, "hellinger")

#Step 3: use this  function to renumber clusters sequentially
renumber.cl <- function(gr) {
  aa <- 1
  gr2 <- rep(1,length(gr))
  for (i in 2:length(gr)) {
    if (gr[i]!=gr[i-1]) aa <- aa+1
    gr2[i] <- aa
  }
  gr2
}

#Step 4: Make tree, select number of branches or use xv="min" to pick the best-fitting tree
#output should pop up in XQuartz  is xv="pick" is selected

#species tree 
everybody.hel.seq <- mvpart(as.matrix(everybody_hel) ~ percent_LOI+percent_clay+season+salinity_class+salinity_anomaly+salinity,
                            data=env_everybody, xv="pick",
                            margin=0.08, xvmult=100, legend=F) 
#family tree
family.hel.seq <- mvpart(as.matrix(family_hel) ~ percent_LOI+percent_clay+season++salinity_class+salinity,
                         data=env_everybody, xv="pick",
                         margin=0.08, xvmult=100, legend=F) 

#feeding guild tree
fg.hel.seq <- mvpart(as.matrix(fg_hel) ~ percent_LOI+percent_clay+season++salinity_class+salinity,
                     data=env_everybody, xv="pick",
                     margin=0.08, xvmult=100, legend=TRUE) 


#cowplot these together 
plot_grid(everybody.hel.seq,                                    # Create grid of plots
          fg.hel.seq,
          ncol = 2,
          hjust = - 1.5,
          labels = c("species", "feeding guilds"))

#extract relative abundances from tree
#species tree, no zeros 

library(ggplot2)
sp3<-fg.hel.seq$frame #tree as data frame
leafy<-filter(sp3, var=="<leaf>") # filter rows by those that represent one leaf of the tree
leafy_comm <-as.data.frame(leafy[,9]) #column 9 is a matrix of spp abundances at each leaf
colnames(leafy_comm) <-names(fg.community.nozero) #rename according to original spp data frame
leafy_commt <-as.data.frame(t(leafy_comm)) #transpose df for easier interpretation
leafy_commt

sp3<-everybody.hel.seq$frame #tree as data frame
leafy<-filter(sp3, var=="<leaf>") # filter rows by those that represent one leaf of the tree
leafy_comm <-as.data.frame(leafy[,9]) #column 9 is a matrix of spp abundances at each leaf
colnames(leafy_comm) <-names(everybody) #rename according to original spp data frame
leafy_commt <-as.data.frame(t(leafy_comm)) #transpose df for easier interpretation
leafy_commt



######################################################
#figure S1 - electronic supplementary

#############################

#Species accumulation curve 
sp.community[is.na(sp.community)] = 0
sp1 <-specaccum(sp.community)
sp2 <- specaccum(sp.community, method = "collector")
spec_accum <-plot(sp2, xlab = "samples", ylab = "species")

###################################################
#Whittakers
#Whittaker plots by salinity regime
#omit hippogriffs
salinity_class <-as.factor(sp.env$salinity_class)
sp.env$salinity_class <-as.factor(sp.env$salinity_class)
par(mfrow=c(3,1))

#spring whittaker
spring <-rep(c("Spring"),each=46)
spring_whittaker <-rankabuncomp(spring_comm[,1:109], y=spring_env, factor='salinity_class',
                                scale='abundance', specnames=c(1:2), legend =TRUE)
spring_rankabun <-data.frame(spring_whittaker$Grouping, spring_whittaker$species, spring_whittaker$rank, spring_whittaker$abundance, spring)
colnames(spring_rankabun) <-c("salinity_class", "species", "rank", "abundance", "season")


#summer whittaker
summer <-rep(c("Summer"), each=78)
summer_whittaker <-rankabuncomp(summer_comm[,1:96], y=summer_env, factor='salinity_class',
                                scale='abundance', specnames=c(1:2), legend =TRUE)
summer_rankabun <-data.frame(summer_whittaker$Grouping, summer_whittaker$species, summer_whittaker$rank, summer_whittaker$abundance, summer)
colnames(summer_rankabun) <-c("salinity_class", "species", "rank", "abundance", "season")

#fall whittaker
fall <-rep(c("Fall"), each=34)
fall_env$salinity_class <-as.factor(fall_env$salinity_class)
fall_whittaker <-rankabuncomp(fall_comm[,1:96], y=fall_env, factor='salinity_class',
                              scale='abundance', legend =TRUE, specnames=c(1:2))
fall_rankabun <-data.frame(as.factor(fall_whittaker$Grouping), fall_whittaker$species, fall_whittaker$rank, fall_whittaker$abundance, fall)
colnames(fall_rankabun) <-c("salinity_class", "species", "rank", "abundance", "season")

#all seasons whittaker
all_whittaker <-rankabuncomp(sp.community, y=sp.env, factor='salinity_class',
                             scale='abundance', legend=TRUE)
all_rankabun <-rbind(spring_rankabun, summer_rankabun, fall_rankabun)

#ggplot Whittakers in one- rankabun per season
colsalinity <- c("mediumvioletred" "lightblue", "gray13")

whittaker <- ggplot(data=all_rankabun, aes(x = rank, y = abundance)) + 
  scale_x_continuous(expand=c(0, 1), sec.axis = dup_axis(labels=NULL, name=NULL, breaks=NULL)) +
  scale_y_continuous(expand=c(0, 1), sec.axis = dup_axis(labels=NULL, name=NULL, breaks=NULL)) +
  geom_line(aes(colour=salinity_class), size=0.5) +
  geom_point(aes(colour=salinity_class, shape=salinity_class), size=2, alpha=0.7) +
  theme_classic()+
  theme(axis.text=element_text(size=15), axis.title=element_text(size=15, face="bold"))+
  theme(legend.text=element_text(size=15), legend.title=element_text(size=15), legend.key.size=unit(0.75, "cm"), legend.position="bottom")+
  theme(strip.text=element_text(size=13, face="bold"), strip.background=element_rect(fill="grey80"))+
  scale_colour_manual(values=colsalinity)+
  facet_wrap(~season, scale='free')+
  labs(x = "rank", y = "abundance", colour = "salinity_class", shape = "salinity_class")
whittaker

#cowplot with spec accumulation curve
#convert specaccum object
my.res <- with(sp2, data.frame(sites, richness))
head(my.res)
accumul <-ggplot(data=my.res, aes(sites, richness))+
  theme_classic()+
  geom_step()+
  xlab("samples") + ylab("species")
accumul

library(cowplot)
plot_grid(accumul, whittaker, nrow=2, ncol=1)


#################################################################
#Data synthesis - null models & multi-model inference

library(MuMIn)

###1. Barataria Bay data
#null model works better without 0 abundance rows 
all_comm <-read.csv("~/Documents/Dissertation/Data/community_nozero.csv") #species + env, zeros removed
species_matrix <-all_comm[1:195,5:113] #subset community
species_matrix[is.na(species_matrix)] <- 0 #replace NA with 0 in df
sp.env.nozero <-data.frame(all_comm$Season, all_comm$Salinity, all_comm$Salinity_class, all_comm$Temperature,
                           all_comm$pH, all_comm$percent_clay, all_comm$percent_silt, all_comm$percent_sand,
                           all_comm$percent_LOI, all_comm$SAL_anomaly)
colnames(sp.env.nozero) <-c("season", "salinity", "salinity_class", "temperature", "pH", "percent_clay",
                            "percent_silt", "percent_sand", "percent_LOI", "sal_anomaly")
sp.env.nozero <-sp.env.nozero[1:195,] #environmental dataset without pesky blank last row

#salinity anomaly vs bc dissimilarity - data for Figure 3
xyplot(site.meanBB ~ abs(all_comm.SAL_anomaly), data=species_env_matrix) #BB data
flow.BB.bray <-meandist(braycurtis, species_env_matrix$all_comm.SAL_anomaly)
plot(flow.BB.bray, kind="histogram", ylim)
flow.mean.BB <-colMeans(flow.BB.bray)
bb.flow.df <- data.frame(flow.mean.BB)
sal.anom.bb <-as.numeric(rownames(bb.flow.df))
bb.flow.df1 <- data.frame(flow.mean.BB, sal.anom.bb)
bb.flow.df1 <-na.omit(bb.flow.df1)
model.bb <-lm(flow.mean.BB ~ abs(sal.anom.bb), data=bb.flow.df1)
plot(flow.mean.BB ~ abs(sal.anom.bb), data=bb.flow.df1)
lines(abs(bb.flow.df1$sal.anom.bb), predict(model.bb), col="cornflowerblue") #line for this LM
summary(model.bb) 

# Bray-Curtis dissimilarity for Barataria Bay - means for each site + stdev

#mean dissimilarity for site pairs - for Figure 1
mean.bray <-meandist(braycurtis, species_env_matrix$all_comm.Site)
plot(mean.bray, kind = "histogram", ylim)
site.meanBB <-colMeans(mean.bray) #mean BC between each site
site1.stdev <-sd(mean.bray[1:15,1]) #stdev of BC between each site
site2.stdev <-sd(mean.bray[1:15,2])
site3.stdev <-sd(mean.bray[1:15,3])
site4.stdev <-sd(mean.bray[1:15,4])
site5.stdev <-sd(mean.bray[1:15,5])
site6.stdev <-sd(mean.bray[1:15,6])
site7.stdev <-sd(mean.bray[1:15,7])
site8.stdev <-sd(mean.bray[1:15,8])
site9.stdev <-sd(mean.bray[1:15,9])
site10.stdev <-sd(mean.bray[1:15,10])
site11.stdev <-sd(mean.bray[1:15,11])
site12.stdev <-sd(mean.bray[1:15,12])
site13.stdev <-sd(mean.bray[1:15,13])
site14.stdev <-sd(mean.bray[1:15,14])
site15.stdev <-sd(mean.bray[1:15,15])

#vector of all stdevs
all.stdev <-c(site1.stdev,site2.stdev,site3.stdev,site4.stdev, site5.stdev,site6.stdev,site7.stdev,site8.stdev,
              site9.stdev, site10.stdev, site11.stdev, site12.stdev, site13.stdev, site14.stdev, site15.stdev)

#create null model for each site separately. Mean=Bexp (expected dissimilarity according to null model)
#Bdep =(Bobs - Bexp)/sd(Bexp)) -departure from null model
#Site 1
site1_null <-nullmodel(species_env_matrix[all_comm$Site==1,1:109], method="swsh_both_r") #using algorithm from Huttunen et al., 2017 (quantitative swap and shuffle model)
sm1 <-simulate(site1_null, nsim=1000) 
summary(sm1) #mean of simulated null model is 0.05515, min=0, max=8. Mean Bexp =0.05515.
sd(sm1) #sd(Bexp)
site.meanBB[1] #Bobs
Bdep1=(site.meanBB[1]-mean(sm1))/sd(sm1)
Bdep1 #1.93

#matrix of Bobs for site 1
Bobs1 <-as.matrix(braycurtis)[species_env_matrix$all_comm.Site==1,species_env_matrix$all_comm.Site==1] 

#matrix of Bdep (departure from null) to be used as response variable
Bdep1_matrix <-(Bobs1-Bexp1)/sd(sm1)
Bdep1_matrix
colMeans(Bdep1_matrix)
mod1 <-lm(Bdep1_matrix ~all_comm$Salinity[Site==1,])
mod1


#Site 2
site2_null <-nullmodel(species_env_matrix[all_comm$Site==2,1:109], method="swsh_samp_r")
sm2 <-simulate(site2_null, nsim=1000)
summary(sm2) #mean=0.02105, max=4. Bexp=0.02105
bray.2 <-vegdist(sm2, method="bray")
sd(sm2) #sd(Bexp)
site.meanBB[2] #Bobs
Bdep2=(site.meanBB[2]-mean(sm2))/sd(sm2)
Bdep2 #

#Site 3
site3_null <-nullmodel(species_env_matrix[all_comm$Site==3,1:109], method="swsh_both_r")
sm3 <-simulate(site3_null, nsim=1000)
summary(sm3) #mean=0.05742, max=10
sd(sm3) #sd(Bexp)
bray.3 <-vegdist(sm3, method="bray")
site.meanBB[3] #Bobs
Bdep3=(site.meanBB[3]-mean(sm3))/sd(sm3)
Bdep3 #1.16

#Site 4
site4_null <-nullmodel(species_env_matrix[all_comm$Site==4,1:109], method="swsh_both_r")
sm4 <-simulate(site4_null, nsim=1000)
summary(sm4) #mean=0.1206, max=49
sd(sm4) #sd(Bexp)
bray.4 <-vegdist(sm4, method="bray")
site.meanBB[4] #Bobs
Bdep4=(site.meanBB[4]-mean(sm4))/sd(sm4)
Bdep4 #0.3487299

#Site 5
site5_null <-nullmodel(species_env_matrix[all_comm$Site==5,1:109], method="swsh_both_r")
sm5 <-simulate(site5_null, nsim=1000)
summary(sm5) 
sd(sm5) #sd(Bexp)
site.meanBB[5] #Bobs
Bdep5=(site.meanBB[5]-mean(sm5))/sd(sm5)
Bdep5 #3.925037

#Site 6
site6_null <-nullmodel(species_env_matrix[all_comm$Site==6,1:109], method="swsh_both_r")
sm6 <-simulate(site6_null, nsim=1000)
summary(sm6) 
sd(sm6) #sd(Bexp)
site.meanBB[6] #Bobs
Bdep6=(site.meanBB[6]-mean(sm6))/sd(sm6)
Bdep6 

#Site 7
site7_null <-nullmodel(species_env_matrix[all_comm$Site==7,1:109], method="swsh_both_r")
sm7 <-simulate(site7_null, nsim=1000)
summary(sm7) 
sd(sm7) #sd(Bexp)
site.meanBB[7] #Bobs
Bdep7=(site.meanBB[7]-mean(sm7))/sd(sm7)
Bdep7 #

#Site 8
site8_null <-nullmodel(species_env_matrix[all_comm$Site==8,1:109], method="swsh_both_r")
sm8 <-simulate(site8_null, nsim=1000)
summary(sm8) 
sd(sm8) #sd(Bexp)
site.meanBB[8] #Bobs
Bdep8=(site.meanBB[8]-mean(sm8))/sd(sm8)
Bdep8 

#Site 9
site9_null <-nullmodel(species_env_matrix[all_comm$Site==9,1:109], method="swsh_both_r")
sm9 <-simulate(site9_null, nsim=1000)
summary(sm9) 
sd(sm9) #sd(Bexp)
site.meanBB[9] #Bobs
Bdep9=(site.meanBB[9]-mean(sm9))/sd(sm9)
Bdep9 

#Site 10
site10_null <-nullmodel(species_env_matrix[all_comm$Site==10,1:109], method="swsh_both_r")
sm10 <-simulate(site10_null, nsim=1000)
summary(sm10) 
sd(sm10) #sd(Bexp)
site.meanBB[10] #Bobs
Bdep10=(site.meanBB[10]-mean(sm10))/sd(sm10)
Bdep10 

#Site 11
site11_null <-nullmodel(species_env_matrix[all_comm$Site==11,1:109], method="swsh_both_r")
sm11 <-simulate(site11_null, nsim=1000)
summary(sm11) 
sd(sm11) #sd(Bexp)
site.meanBB[11] #Bobs
Bdep11=(site.meanBB[11]-mean(sm11))/sd(sm11)
Bdep11 

#Site 12
site12_null <-nullmodel(species_env_matrix[all_comm$Site==12,1:109], method="swsh_both_r")
sm12 <-simulate(site12_null, nsim=1000)
summary(sm4) 
sd(sm12) #sd(Bexp)
site.meanBB[12] #Bobs
Bdep12=(site.meanBB[12]-mean(sm12))/sd(sm12)
Bdep12 

#Site 13
site13_null <-nullmodel(species_env_matrix[all_comm$Site==13,1:109], method="swsh_both_r")
sm13 <-simulate(site13_null, nsim=1000)
summary(sm13) 
sd(sm13) #sd(Bexp)
site.meanBB[13] #Bobs
Bdep13=(site.meanBB[13]-mean(sm13))/sd(sm13)
Bdep13 

#Site 14
site14_null <-nullmodel(species_env_matrix[all_comm$Site==14,1:109], method="swsh_both_r")
sm14 <-simulate(site14_null, nsim=1000)
summary(sm14) 
sd(sm14) #sd(Bexp)
site.meanBB[14] #Bobs
Bdep14=(site.meanBB[14]-mean(sm14))/sd(sm14)
Bdep14 

#Site 15
site15_null <-nullmodel(species_env_matrix[all_comm$Site==15,1:109], method="swsh_both_r")
sm15 <-simulate(site15_null, nsim=1000)
summary(sm15) 
sd(sm15) #sd(Bexp)
site.meanBB[15] #Bobs
Bdep15=(site.meanBB[15]-mean(sm15))/sd(sm15)
Bdep15 

#SES (departure from null at each site)
Bdeps <-c(Bdep1, Bdep2, Bdep3, Bdep4, Bdep5, Bdep6, Bdep7, Bdep8, Bdep9, Bdep10, Bdep11, Bdep12, Bdep13, Bdep14, Bdep15)
mean.Bdep <-mean(Bdeps)
sd.BDEP <-sd(Bdeps)

#plot of SES vs site - for Figure 2
plot(Bdeps ~c(1:15), ylab = "departure from null expectation", xlab = "site", ylim=c(0,5),
     pch=19, col="cornflowerblue")
text(x=15,y=4.5, "B")
abline(1.96*sd.BDEP,0)

#compare BB
AIC(model.bb, model.bb.nls) #lm wins
##########################################################
#Texas data - null models

#read in comm and env files
montagna_env <-read.csv("~/Documents/Dissertation/Data/Montagna/Montagna_env1.csv")
montagna_comm <-read.csv("~/Documents/Dissertation/manuscripts/Montagna_comm.csv")

#replace NAs in df with 0
montagna_comm[is.na(montagna_comm)] <- 0
sp <-as.data.frame(montagna_comm[,5:150])

#view datasets
View(sp)
View(montagna_env)


#calculate salinity anomaly based on mean annual average salinity
montagna_env$SAL
avg_SAL1 <-mean(montagna_env[montagna_env$Year==2020 & montagna_env$Site ==1,montagna_env$SAL])
avg_SAL1

#env data frame with salinity anomaly to represent flow extreme
env <-data.frame(montagna_env$Season, montagna_env$SAL, montagna_env$TEMP, montagna_env$CLAY, montagna_env$TOC, montagna_env$Est, montagna_env$STA, montagna_env$SAL_anomaly)
colnames(env) <-c("season", "salinity", "temp", "clay", "TOC", "estuary", "station", "salinity_anomaly")

#Bray-Curtis and Bobs for each site pair
matrix1<-data.frame(sp, env)
#remove all no zero sum abundance rows and corresponding env vars
matrix <-matrix1[rowSums(matrix1[,1:108])>0,]
BC <- vegdist(matrix[,1:108], method="bray") #species matrix
matrix$Est_STA <- paste(matrix$estuary, matrix$station, sep="") #add column with merged est & sta 

#mean dissimilarity for salinity anomaly pairs
matrix$salinity_anomaly <-as.factor(matrix$salinity_anomaly)
flow.BC <-meandist(BC, matrix$salinity_anomaly) #full BC matrix
plot(flow.BC, kind="histogram", ylim)
flow.mean <-colMeans(flow.BC) #mean dissimilarity for each pair

flow.df <-data.frame(flow.mean) 
sal_anomaly <-as.numeric(rownames(flow.df)) #row names are salinity values
flow.df1 <-data.frame(flow.mean, sal_anomaly)
flow.df1 <-na.omit(flow.df1)#remove rows where flow mean is NA

#plot mean BC for each salinity anomaly pair
plot(flow.mean ~ abs(sal_anomaly), data =flow.df1)
points(flow.mean, abs(sal_anomaly), col="green", type="l", lwd=1)
model.tx <-lm(flow.mean ~ abs(sal_anomaly), data=flow.df1) #lm is basically 0. slope is -0.004. LM is not appropriate here
summary(model.tx)
#use nls model to fit negative logarithm (looks like exponential decay)
model.flow <- nls(flow.mean ~ a*exp(k*abs(sal_anomaly)), data = flow.df1, start=list(a=0.001, k=0.3), trace=T)
#the big number 60.13is the sum of squares for our initial start values (0.001,0.3). R goes on to adjust
#the last number 15.13 is the sum of squares for our start values a and k. a is 0.8009 and k is -0.005558

plot(resid(model.flow)~fitted(model.flow),xlab='Fitted Value',ylab='Residual') #residual v. fitted 

#plot model - y=ae^(kx) - don't use this nls. turns out you can just use an lm()
plot(flow.mean ~ abs(sal_anomaly), data =flow.df1) #plot original points
yFitted.tx <-predict(model.flow) #
lines(abs(flow.df1$sal_anomaly), predict(model.flow), col= "coral3") #plot trendline
summary(model.flow)

#which model - use AIC to compare lm() results and nls() results
#compare TX
AIC(model.tx, model.flow) #nls wins

#mean dissimilarity for site pairs - TX - for figure 1
mean.BC <-meandist(BC, matrix$Est_STA) #mean BC grouped by est + station
plot(mean.BC, kind = "histogram", ylim)
site.mean <-colMeans(mean.BC) #mean BC between each station
BRA.stdev <-sd(mean.BC[1:15,1]) #stdev of BC between each station
BRB.stdev <-sd(mean.BC[1:15,2])
BRC.stdev <-sd(mean.BC[1:15,3])
CBA.stdev <-sd(mean.BC[1:15,4])
CBB.stdev <-sd(mean.BC[1:15,5])
CBC.stdev <-sd(mean.BC[1:15,6])
LCA.stdev <-sd(mean.BC[1:15,7])
LCB.stdev <-sd(mean.BC[1:15,8])
LCC.stdev <-sd(mean.BC[1:15,9])
LCD.stdev <-sd(mean.BC[1:15,10])
RGA.stdev <-sd(mean.BC[1:15,11])
RGB.stdev <-sd(mean.BC[1:15,12])
RGC.stdev <-sd(mean.BC[1:15,13])
RGD.stdev <-sd(mean.BC[1:15,14])
RGE.stdev <-sd(mean.BC[1:15,15])
#vector of all stdevs
stdev <-c(BRA.stdev, BRB.stdev, BRC.stdev, CBA.stdev, CBB.stdev, CBC.stdev, LCA.stdev, LCB.stdev, LCC.stdev, 
          LCD.stdev, RGA.stdev, RGB.stdev, RGC.stdev, RGD.stdev, RGE.stdev)
               

################################################
#null models- TEXAS lagoons

#create null model for each est + station separately. Mean=Bexp (expected dissimilarity according to null model)
#Bdep =(Bobs - Bexp)/sd(Bexp)) -departure from null model

#BR - A
BR_A_null <-nullmodel(matrix[matrix$Est_STA=="BRA",1:142], method="swsh_both_r") #using algorithm from Huttunen et al., 2017 (quantitative swap and shuffle model)
smBR_A <-simulate(BR_A_null, nsim=1000) 
summary(smBR_A) #mean of simulated null model. Mean Bexp =0.12654.
sd(smBR_A) #sd(Bexp)
site.mean[1] #Bobs
Bdep_BR_A=(site.mean[1]-mean(smBR_A))/sd(smBR_A)
Bdep_BR_A #1.93

#matrix of Bexp for BR-A
simBR_A <-smBR_A[,,1]
rowSums(simBR_A) #check to see in there are any zero count rows
Bexp_BR_A <-as.matrix(vegdist(simBR_A, method="bray", na.rm=TRUE))#calculate bray-curtis for simulated matrix  

#matrix of Bobs for BR-A - don't need this
BobsBRA <-as.matrix(BC)[matrix$Est_STA=="BRA",matrix$Est_STA=="BRA"] 

#matrix of Bdep (departure from null) to be used as response variable
Bdep_BRA_matrix <-(BobsBRA-Bexp_BR_A)/sd(smBR_A)
Bdep_BRA_matrix
colMeans(Bdep_BRA_matrix)
mod_BRA <-lm(Bdep_BRA_matrix ~matrix$salinity[matrix$Est_STA=="BRA",])
mod_BRA

#BR-B
BR_B_null <-nullmodel(matrix[matrix$Est_STA=="BRB",1:142], method="swsh_both_r") #using algorithm from Huttunen et al., 2017 (quantitative swap and shuffle model)
smBR_B <-simulate(BR_B_null, nsim=1000) 
summary(smBR_B) #mean of simulated null model. Mean Bexp =0.18634.
sd(smBR_B) #sd(Bexp)
site.mean[2] #Bobs
Bdep_BR_B=(site.mean[2]-mean(smBR_B))/sd(smBR_B)
Bdep_BR_B #0.3900036

#BR-C
BR_C_null <-nullmodel(matrix[matrix$Est_STA=="BRC",1:142], method="swsh_both_r") #using algorithm from Huttunen et al., 2017 (quantitative swap and shuffle model)
smBR_C <-simulate(BR_C_null, nsim=1000) 
summary(smBR_C) #mean of simulated null model. Mean Bexp =0.22068.
sd(smBR_C) #sd(Bexp)
site.mean[3] #Bobs
Bdep_BR_C=(site.mean[3]-mean(smBR_C))/sd(smBR_C)
Bdep_BR_C #0.2808904

#CB-A
CB_A_null <-nullmodel(matrix[matrix$Est_STA=="CBA",1:142], method="swsh_both_r") #using algorithm from Huttunen et al., 2017 (quantitative swap and shuffle model)
smCB_A <-simulate(CB_A_null, nsim=1000) 
summary(smCB_A) #mean of simulated null model. Mean Bexp =0.67945.
sd(smCB_A) #sd(Bexp)
site.mean[4] #Bobs
Bdep_CB_A=(site.mean[4]-mean(smCB_A))/sd(smCB_A)
Bdep_CB_A #0.001501

#CB-B
CB_B_null <-nullmodel(matrix[matrix$Est_STA=="CBB",1:142], method="swsh_both_r") #using algorithm from Huttunen et al., 2017 (quantitative swap and shuffle model)
smCB_B <-simulate(CB_B_null, nsim=1000) 
summary(smCB_B) #mean of simulated null model. Mean Bexp =0.22134.
sd(smCB_B) #sd(Bexp)
site.mean[5] #Bobs
Bdep_CB_B=(site.mean[5]-mean(smCB_B))/sd(smCB_B)
Bdep_CB_B #0.57144

#CB-C
CB_C_null <-nullmodel(matrix[matrix$Est_STA=="CBC",1:142], method="swsh_both_r") #using algorithm from Huttunen et al., 2017 (quantitative swap and shuffle model)
smCB_C <-simulate(CB_C_null, nsim=1000) 
summary(smCB_C) #mean of simulated null model. Mean Bexp =0.88713.
sd(smCB_C) #sd(Bexp)
site.mean[6] #Bobs
Bdep_CB_C=(site.mean[6]-mean(smCB_C))/sd(smCB_C)
Bdep_CB_C #-0.08272

#LC-A #this one is being wonky
LC_A_null <-nullmodel(matrix[matrix$Est_STA=="LCA",1:142], method="swsh_both_r") #using algorithm from Huttunen et al., 2017 (quantitative swap and shuffle model)
smLC_A <-simulate(LC_A_null, nsim=1000) 
summary(smLC_A) #mean of simulated null model. Mean Bexp =0.88713.
sd(smLC_A) #sd(Bexp)
site.mean[7] #Bobs
Bdep_LC_A=(site.mean[7]-mean(smLC_A))/sd(smLC_A)
Bdep_LC_A #-0.08272

#LC-B
LC_B_null <-nullmodel(matrix[matrix$Est_STA=="LCB",1:142], method="swsh_both_r") #using algorithm from Huttunen et al., 2017 (quantitative swap and shuffle model)
smLC_B <-simulate(LC_B_null, nsim=1000) 
summary(smLC_B) #mean of simulated null model. Mean Bexp =0.88713.
sd(smLC_B) #sd(Bexp)
site.mean[8] #Bobs
Bdep_LC_B=(site.mean[8]-mean(smLC_B))/sd(smLC_B)
Bdep_LC_B #-0.08272

#LC-C
LC_C_null <-nullmodel(matrix[matrix$Est_STA=="LCC",1:142], method="swsh_both_r") #using algorithm from Huttunen et al., 2017 (quantitative swap and shuffle model)
smLC_C <-simulate(LC_C_null, nsim=1000) 
summary(smLC_C) #mean of simulated null model. Mean Bexp =0.88713.
sd(smLC_C) #sd(Bexp)
site.mean[9] #Bobs
Bdep_LC_C=(site.mean[9]-mean(smLC_C))/sd(smLC_C)
Bdep_LC_C #-0.08272

#LC-D
LC_D_null <-nullmodel(matrix[matrix$Est_STA=="LCD",1:142], method="swsh_both_r") #using algorithm from Huttunen et al., 2017 (quantitative swap and shuffle model)
smLC_D <-simulate(LC_D_null, nsim=1000) 
summary(smLC_D) #mean of simulated null model. Mean Bexp =0.88713.
sd(smLC_D) #sd(Bexp)
site.mean[10] #Bobs
Bdep_LC_D=(site.mean[10]-mean(smLC_D))/sd(smLC_D)
Bdep_LC_D #-0.08272

#RG-A
RG_A_null <-nullmodel(matrix[matrix$Est_STA=="RGA",1:142], method="swsh_both_r") #using algorithm from Huttunen et al., 2017 (quantitative swap and shuffle model)
smRG_A <-simulate(RG_A_null, nsim=1000) 
summary(smRG_A) #mean of simulated null model. Mean Bexp =0.88713.
sd(smRG_A) #sd(Bexp)
site.mean[11] #Bobs
Bdep_RG_A=(site.mean[11]-mean(smRG_A))/sd(smRG_A)
Bdep_RG_A #-0.08272

#RG-B
RG_B_null <-nullmodel(matrix[matrix$Est_STA=="RGB",1:142], method="swsh_both_r") #using algorithm from Huttunen et al., 2017 (quantitative swap and shuffle model)
smRG_B <-simulate(RG_B_null, nsim=1000) 
summary(smRG_B) #mean of simulated null model. Mean Bexp =0.88713.
sd(smRG_B) #sd(Bexp)
site.mean[12] #Bobs
Bdep_RG_B=(site.mean[12]-mean(smRG_B))/sd(smRG_B)
Bdep_RG_B #-0.08272

#RG-C
RG_C_null <-nullmodel(matrix[matrix$Est_STA=="RGC",1:142], method="swsh_both_r") #using algorithm from Huttunen et al., 2017 (quantitative swap and shuffle model)
smRG_C <-simulate(RG_C_null, nsim=1000) 
summary(smRG_C) #mean of simulated null model. Mean Bexp =0.88713.
sd(smRG_C) #sd(Bexp)
site.mean[13] #Bobs
Bdep_RG_C=(site.mean[13]-mean(smRG_C))/sd(smRG_C)
Bdep_RG_C #-0.08272

#RG-D
RG_D_null <-nullmodel(matrix[matrix$Est_STA=="RGD",1:142], method="swsh_both_r") #using algorithm from Huttunen et al., 2017 (quantitative swap and shuffle model)
smRG_D <-simulate(RG_D_null, nsim=1000) 
summary(smRG_D) #mean of simulated null model. Mean Bexp =0.88713.
sd(smRG_D) #sd(Bexp)
site.mean[14] #Bobs
Bdep_RG_D=(site.mean[14]-mean(smRG_D))/sd(smRG_D)
Bdep_RG_D #-0.08272

#RG-E
RG_E_null <-nullmodel(matrix[matrix$Est_STA=="RGE",1:142], method="swsh_both_r") #using algorithm from Huttunen et al., 2017 (quantitative swap and shuffle model)
smRG_E <-simulate(RG_E_null, nsim=1000) 
summary(smRG_E) #mean of simulated null model. Mean Bexp =0.88713.
sd(smRG_E) #sd(Bexp)
site.mean[15] #Bobs
Bdep_RG_E=(site.mean[15]-mean(smRG_E))/sd(smRG_E)
Bdep_RG_E #-0.08272

####################################################################################################
#New Zealand data 
#######################################################################################################
#read in data
stubbington_all <- read.csv("~/Documents/Dissertation/Data/Stubbington_2017/stubbington__selwyn.csv", header=TRUE)
selwyn_comm <-stubbington_all[, 5:93]
selwyn_env <-stubbington_all[,1:4]

#format data
env_selwyn <-data.frame(selwyn_env$Month, selwyn_env$Year, selwyn_env$Site, selwyn_env$Flow.intermittence....)
colnames(env_selwyn) <-c("month", "year", "site", "intermittance")

Bray-Curtis and Bobs for each unique flow scheme
#Stubbington 2017

matrix_selwyn<-data.frame(selwyn_comm, env_selwyn)
BC_selwyn <- vegdist(matrix_selwyn[,1:89], method="bray") #species matrix
matrix_selwyn$site 


#mean dissimilarity for site pairs
mean.BC.selwyn <-meandist(BC_selwyn, matrix_selwyn$site) #mean BC grouped by est + station
plot(mean.BC.selwyn, kind = "histogram", ylim)
site.mean.selwyn <-colMeans(mean.BC.selwyn) #mean BC between each station
bealey.stdev <-sd(mean.BC.selwyn[1:16,1]) #stdev of BC between each station
coalgate.stdev <-sd(mean.BC.selwyn[1:16,2])
coesford.stdev <-sd(mean.BC.selwyn[1:16,3])
gillanders.stdev <-sd(mean.BC.selwyn[1:16,4])
hawkins.stdev <-sd(mean.BC.selwyn[1:16,5])
highfield.stdev <-sd(mean.BC.selwyn[1:16,6])
hororata.stdev <-sd(mean.BC.selwyn[1:16,7])
mcgregor.stdev <-sd(mean.BC.selwyn[1:16,8])
oldbridgeroad.stdev <-sd(mean.BC.selwyn[1:16,9])
oldsouthroad.stdev <-sd(mean.BC.selwyn[1:16,10])
raywell.stdev <-sd(mean.BC.selwyn[1:16,11])
ridgens.stdev <-sd(mean.BC.selwyn[1:16,12])
scottsroad.stdev <-sd(mean.BC.selwyn[1:16,13])
waterford.stdev <-sd(mean.BC.selwyn[1:16,14])
westernras.stdev <-sd(mean.BC.selwyn[1:16,15])
withells.stdev <-sd(mean.BC.selwyn[1:16,16])

#plot of mean BC versus salinity anomaly
mean.BC.flow <-meandist(BC_selwyn, matrix_selwyn$intermittance) #BC matrix of sites by flow intermittance
flow.mean.selwyn <-colMeans(mean.BC.flow) #mean BC at each flow regime
plot(flow.mean.selwyn ~ intermit)
lm(flow.mean.selwyn)
intermit <-c(0,7,45,63,65,69,72,80,83,84,85,89)
intermit1 <-intermit/2.857
selwyn.flow.df <- data.frame(flow.mean.selwyn, intermit1)
model.selwyn <-lm(flow.mean.selwyn ~ intermit1, data=selwyn.flow.df)
lines(intermit1, predict(model.selwyn), col="#004D40") #line for this LM
summary(model.selwyn)  
#vector of all stdevs
stdev.selwyn <-c(bealey.stdev, coalgate.stdev, coesford.stdev, gillanders.stdev, hawkins.stdev, highfield.stdev,
                 hororata.stdev, mcgregor.stdev, oldbridgeroad.stdev, oldsouthroad.stdev, raywell.stdev, ridgens.stdev,
                 scottsroad.stdev, waterford.stdev, westernras.stdev, withells.stdev)

#null models - Datry 2007, 2014 - New Zealand data - one per site

#create null model for each est + station separately. Mean=Bexp (expected dissimilarity according to null model)
#Bdep =(Bobs - Bexp)/sd(Bexp)) -departure from null model


#Bealey
bealey_null <-nullmodel(matrix_selwyn[matrix_selwyn$site=="Bealey",1:89], method="swsh_both_r") #using algorithm from Huttunen et al., 2017 (quantitative swap and shuffle model)
sm_bealey <-simulate(bealey_null, nsim=1000) 
summary(sm_bealey) #mean of simulated null model. Mean Bexp =0.18634.
sd(sm_bealey) #sd(Bexp)
site.mean.selwyn[1] #Bobs
Bdep_bealey=(site.mean.selwyn[1]-mean(sm_bealey))/sd(sm_bealey)
Bdep_bealey #

#Coalgate
coalgate_null <-nullmodel(matrix_selwyn[matrix_selwyn$site=="Coalgate",1:89], method="swsh_both_r") #using algorithm from Huttunen et al., 2017 (quantitative swap and shuffle model)
sm_coalgate <-simulate(coalgate_null, nsim=1000) 
summary(sm_coalgate) #mean of simulated null model. Mean Bexp =0.18634.
sd(sm_coalgate) #sd(Bexp)
site.mean.selwyn[2] #Bobs
Bdep_coalgate=(site.mean.selwyn[2]-mean(sm_coalgate))/sd(sm_coalgate)
Bdep_coalgate #

#Coes Ford
coesford_null <-nullmodel(matrix_selwyn[matrix_selwyn$site=="Coes Ford",1:89], method="swsh_both_r") #using algorithm from Huttunen et al., 2017 (quantitative swap and shuffle model)
sm_coesford <-simulate(coesford_null, nsim=1000) 
summary(sm_coesford) #mean of simulated null model. Mean Bexp =0.18634.
sd(sm_coesford) #sd(Bexp)
site.mean.selwyn[3] #Bobs
Bdep_coesford=(site.mean.selwyn[3]-mean(sm_coesford))/sd(sm_coesford)
Bdep_coesford #

#Gillanders
gillanders_null <-nullmodel(matrix_selwyn[matrix_selwyn$site=="Gillanders",1:89], method="swsh_both_r") #using algorithm from Huttunen et al., 2017 (quantitative swap and shuffle model)
sm_gillanders <-simulate(gillanders_null, nsim=1000) 
summary(sm_gillanders) #mean of simulated null model. Mean Bexp =0.18634.
sd(sm_gillanders) #sd(Bexp)
site.mean.selwyn[4] #Bobs
Bdep_gillanders=(site.mean.selwyn[4]-mean(sm_gillanders))/sd(sm_gillanders)
Bdep_gillanders #

#Hawkins
hawkins_null <-nullmodel(matrix_selwyn[matrix_selwyn$site=="Hawkins",1:89], method="swsh_both_r") #using algorithm from Huttunen et al., 2017 (quantitative swap and shuffle model)
sm_hawkins <-simulate(hawkins_null, nsim=1000) 
summary(sm_hawkins) #mean of simulated null model. Mean Bexp =0.18634.
sd(sm_hawkins) #sd(Bexp)
site.mean.selwyn[5] #Bobs
Bdep_hawkins=(site.mean.selwyn[5]-mean(sm_hawkins))/sd(sm_hawkins)
Bdep_hawkins #

#Highfield
highfield_null <-nullmodel(matrix_selwyn[matrix_selwyn$site=="Highfield",1:89], method="swsh_both_r") #using algorithm from Huttunen et al., 2017 (quantitative swap and shuffle model)
sm_highfield <-simulate(highfield_null, nsim=1000) 
summary(sm_highfield) #mean of simulated null model. Mean Bexp =0.18634.
sd(sm_highfield) #sd(Bexp)
site.mean.selwyn[6] #Bobs
Bdep_highfield=(site.mean.selwyn[6]-mean(sm_highfield))/sd(sm_highfield)
Bdep_highfield #

#Hororata
hororata_null <-nullmodel(matrix_selwyn[matrix_selwyn$site=="Hororata",1:89], method="swsh_both_r") #using algorithm from Huttunen et al., 2017 (quantitative swap and shuffle model)
sm_hororata <-simulate(hororata_null, nsim=1000) 
summary(sm_hororata) #mean of simulated null model. Mean Bexp =0.18634.
sd(sm_hororata) #sd(Bexp)
site.mean.selwyn[7] #Bobs
Bdep_hororata=(site.mean.selwyn[7]-mean(sm_hororata))/sd(sm_hororata)
Bdep_hororata #

#McGregor
mcgregor_null <-nullmodel(matrix_selwyn[matrix_selwyn$site=="McGregor",1:89], method="swsh_both_r") #using algorithm from Huttunen et al., 2017 (quantitative swap and shuffle model)
sm_mcgregor <-simulate(mcgregor_null, nsim=1000) 
summary(sm_mcgregor) #mean of simulated null model. Mean Bexp =0.18634.
sd(sm_mcgregor) #sd(Bexp)
site.mean.selwyn[8] #Bobs
Bdep_mcgregor=(site.mean.selwyn[8]-mean(sm_mcgregor))/sd(sm_mcgregor)
Bdep_mcgregor #

#Old Bridge Road
oldbridgeroad_null <-nullmodel(matrix_selwyn[matrix_selwyn$site=="Old Bridge Road",1:89], method="swsh_both_r") #using algorithm from Huttunen et al., 2017 (quantitative swap and shuffle model)
sm_oldbridgeroad <-simulate(oldbridgeroad_null, nsim=1000) 
summary(sm_oldbridgeroad) #mean of simulated null model. Mean Bexp =0.18634.
sd(sm_oldbridgeroad) #sd(Bexp)
site.mean.selwyn[9] #Bobs
Bdep_oldbridgeroad=(site.mean.selwyn[9]-mean(sm_oldbridgeroad))/sd(sm_oldbridgeroad)
Bdep_oldbridgeroad #

#Old South Road
oldsouthroad_null <-nullmodel(matrix_selwyn[matrix_selwyn$site=="Old South Road",1:89], method="swsh_both_r") #using algorithm from Huttunen et al., 2017 (quantitative swap and shuffle model)
sm_oldsouthroad <-simulate(oldsouthroad_null, nsim=1000) 
summary(sm_oldsouthroad) #mean of simulated null model. Mean Bexp =0.18634.
sd(sm_oldsouthroad) #sd(Bexp)
site.mean.selwyn[10] #Bobs
Bdep_oldsouthroad=(site.mean.selwyn[10]-mean(sm_oldsouthroad))/sd(sm_oldsouthroad)
Bdep_oldsouthroad #

#Raywell
raywell_null <-nullmodel(matrix_selwyn[matrix_selwyn$site=="Raywell",1:89], method="swsh_both_r") #using algorithm from Huttunen et al., 2017 (quantitative swap and shuffle model)
sm_raywell <-simulate(raywell_null, nsim=1000) 
summary(sm_raywell) #mean of simulated null model. Mean Bexp =0.18634.
sd(sm_raywell) #sd(Bexp)
site.mean.selwyn[11] #Bobs
Bdep_raywell=(site.mean.selwyn[11]-mean(sm_raywell))/sd(sm_raywell)
Bdep_raywell #

#Ridgens
ridgens_null <-nullmodel(matrix_selwyn[matrix_selwyn$site=="Ridgens",1:89], method="swsh_both_r") #using algorithm from Huttunen et al., 2017 (quantitative swap and shuffle model)
sm_ridgens <-simulate(ridgens_null, nsim=1000) 
summary(sm_ridgens) #mean of simulated null model. Mean Bexp =0.18634.
sd(sm_ridgens) #sd(Bexp)
site.mean.selwyn[12] #Bobs
Bdep_ridgens=(site.mean.selwyn[12]-mean(sm_ridgens))/sd(sm_ridgens)
Bdep_ridgens #

#Scotts Road
scottsroad_null <-nullmodel(matrix_selwyn[matrix_selwyn$site=="Scotts Road",1:89], method="swsh_both_r") #using algorithm from Huttunen et al., 2017 (quantitative swap and shuffle model)
sm_scottsroad <-simulate(scottsroad_null, nsim=1000) 
summary(sm_scottsroad) #mean of simulated null model. Mean Bexp =0.18634.
sd(sm_scottsroad) #sd(Bexp)
site.mean.selwyn[13] #Bobs
Bdep_scottsroad=(site.mean.selwyn[13]-mean(sm_scottsroad))/sd(sm_scottsroad)
Bdep_scottsroad #

#Waterford
waterford_null <-nullmodel(matrix_selwyn[matrix_selwyn$site=="Waterford",1:89], method="swsh_both_r") #using algorithm from Huttunen et al., 2017 (quantitative swap and shuffle model)
sm_waterford <-simulate(waterford_null, nsim=1000) 
summary(sm_waterford) #mean of simulated null model. Mean Bexp =0.18634.
sd(sm_waterford) #sd(Bexp)
site.mean.selwyn[14] #Bobs
Bdep_waterford=(site.mean.selwyn[14]-mean(sm_waterford))/sd(sm_waterford)
Bdep_waterford #

#Westernras
westernras_null <-nullmodel(matrix_selwyn[matrix_selwyn$site=="Westernras",1:89], method="swsh_both_r") #using algorithm from Huttunen et al., 2017 (quantitative swap and shuffle model)
sm_westernras <-simulate(westernras_null, nsim=1000) 
summary(sm_westernras) #mean of simulated null model. Mean Bexp =0.18634.
sd(sm_westernras) #sd(Bexp)
site.mean.selwyn[15] #Bobs
Bdep_westernras=(site.mean.selwyn[15]-mean(sm_westernras))/sd(sm_westernras)
Bdep_westernras #

#Withells
withells_null <-nullmodel(matrix_selwyn[matrix_selwyn$site=="Withells",1:89], method="swsh_both_r") #using algorithm from Huttunen et al., 2017 (quantitative swap and shuffle model)
sm_withells <-simulate(withells_null, nsim=1000) 
summary(sm_withells) #mean of simulated null model. Mean Bexp =0.18634.
sd(sm_withells) #sd(Bexp)
site.mean.selwyn[16] #Bobs
Bdep_withells=(site.mean.selwyn[16]-mean(sm_withells))/sd(sm_withells)
Bdep_withells #

#SES (departure from null at each site)
Bdep_selwyn <-c(Bdep_bealey, Bdep_coalgate, Bdep_coesford, Bdep_gillanders, Bdep_hawkins, Bdep_highfield, Bdep_hororata, Bdep_mcgregor, Bdep_oldbridgeroad,
                Bdep_oldsouthroad, Bdep_raywell, Bdep_ridgens, Bdep_scottsroad, Bdep_waterford, Bdep_westernras, Bdep_withells)
mean.Bdep_selwyn <-mean(Bdep_selwyn)
sd.BDEP_S <-sd(Bdep_selwyn)

#######################################################################
###FIGURE 1: plot of BC dissimilarity and site - upper panel of figure 
plot_Bobs <-plot(site.meanBB ~ c(1:15), ylab = "observed Bray-Curtis dissimilarity", xlab = "site/station", ylim=c(0,1), pch=16, colour="lightblue1")
#text(x=15, y=0.92, "A")
z=1:15 #just the number of error bars
arrows(z, y0=site.meanBB-all.stdev, z, y1=site.meanBB+all.stdev, length=0.05, angle=90, code=3) #add stdev as error bar
points(site.mean, pch=19, col="gray4")
z=1:15 #just the number of error bars
arrows(z, y0=site.mean-stdev, z, y1=site.mean+stdev, length=0.05, angle=90, code=3) #add stdev as error bars
#z2=1:11
#arrows(z2, y0=site.mean.ell-stdev.ell, z2, y1=site.mean.ell+stdev.ell, length=0.05, angle=90, code=3) #add stdev as error bars
#points(site.mean.ell, pch=19, col="grey50")
z3=1:16
arrows(z3, y0=site.mean.selwyn-stdev.selwyn, z3, y1=site.mean.selwyn+stdev.selwyn, length=0.05, angle=90, code=3, xpd=FALSE)
points(site.mean.selwyn, pch=19, col="seagreen4", xpd=FALSE)
legend("bottomright", legend=c("lagoon", "river", "estuary"), cex=1.5, fill=c("gray4", "seagreen", "lightblue2"))

###########################################################
### FIGURE 2. Departure from null model vs. 2 SE line (marks significant difference from random chance)
#legend info
loc <- c("estuary", "lagoon", "river")
color <-c("lightblue", "gray4", "seagreen4")

#bdep df
bdep_df <-data.frame(Bdeps, Bdep_montagna, Bdep_selwyn[1:15], c(1:15))
colnames(bdep_df) <-c("barataria", "texas", "selwyn", "site")

#figure
barataria <- ggplot(data=bdep_df)+
  geom_point(mapping=aes(x=site, y=barataria, size=3), pch=21, fill="lightblue", colour="cornflowerblue")+
  geom_abline(slope=0,intercept=1.96*sd.BDEP,lwd=0.5, color="lightblue")+
  geom_area(aes(c(1:15), 1.96*sd.BDEP), fill = 'lightblue', alpha = 0.2)+
  theme_classic()+
  annotate("text", x=14, y=1.5, label= "estuary", color="cornflowerblue", size=10)+
  theme(legend.position="none")+
  theme(axis.title=element_text(size=15,face="bold"))+
  scale_y_continuous(name="", limits=c(0, 1.55))
barataria

texas <- ggplot(data=bdep_df)+
  geom_point(mapping=aes(x=site, y=abs(texas), size=3), color="gray4")+
  geom_abline(slope=0,intercept=1.96*sd.BDEP_M,lwd=0.5, color="gray4")+
  geom_area(aes(c(1:15), 1.96*sd.BDEP_M), fill = 'gray4', alpha = 0.2)+
  theme_classic()+
  annotate("text", x=13.5, y=1.05, label= "lagoons", color="gray4", size=10)+
  theme(legend.position="none")+
  theme(axis.title=element_text(size=15,face="bold"))+
  scale_y_continuous(name="departure from null model", limits=c(0, 1.1))
texas

selwyn <- ggplot(data=bdep_df)+
  geom_point(mapping=aes(x=site, y=abs(selwyn), size=3), color="seagreen4")+
  geom_abline(slope=0,intercept=1.96*sd.BDEP_S,lwd=0.5, color="seagreen4")+
  geom_area(aes(c(1:15), 1.96*sd.BDEP_S), fill = 'seagreen4', alpha = 0.2)+
  theme_classic()+
  annotate("text", x=14, y=0.55, label= "river", color="seagreen4", size=10)+
  theme(legend.position="none")+
  theme(axis.title=element_text(size=15,face="bold"))+
  scale_y_continuous(name="", limits=c(0, 0.6))

selwyn

#patchwork to put plots together!
library(patchwork)


(barataria /
    texas /
    selwyn)


############################################################
#FIGURE 3: conceptual diagram of flow extreme vs. BC dissimilarity  - all datasets included
library(yarrr)
plot(flow.mean ~ abs(sal_anomaly), data =flow.df1, col=alpha(0.4), ylim=c(0,1), xlim=c(0,35), ylab="mean Bray-Curtis dissimilarity", 
     xlab = "flow extreme") #baseline plot
#points
points(flow.mean ~abs(sal_anomaly), data=flow.df1, pch=16, col=transparent("gray4", trans.val=0.3), cex=0.5)
points(flow.mean.selwyn~ intermit1, data=selwyn.flow.df, pch=16, col=transparent("seagreen4", trans.val=0.1), cex=0.5)
points(flow.mean.BB ~ abs(sal.anom.bb), data=bb.flow.df1, pch=16, col=transparent("lightblue", trans.val=0.05), cex=0.5)
#lines (nls models)
lines(abs(flow.df1$sal_anomaly), predict(model.flow), col= "gray4", lwd=2) #Texas BC vs sal anom
lines(intermit1, predict(model.selwyn), col="seagreen4", lwd=2) #Selwyn BC vs intermittance
lines(abs(bb.flow.df1$sal.anom.bb), predict(model.bb), col="lightblue", lwd=2) #BB BC vs. sal anom
legend(28,0.20, legend=c("lagoon", "river", "estuary"), cex=2,fill=c("gray4", "seagreen", "lightblue"))


###########################################################
##### lm of env. parameters with dissimilarity (Bobs) as response variable

##check for correlation between variables and latitude - accounting for spatial autocorrelation
cor(all_comm$Salinity, all_comm$latitude)
#-0.73865, so do not include lat in lms

env_vars <-data.frame(all_comm$Salinity, all_comm$Site, all_comm$Season, all_comm$percent_clay,
                      all_comm$percent_LOI, all_comm$Temperature, all_comm$Bobs, all_comm$Bdep, simpson, all_comm$Corophium.louisianum, 
                      all_comm$Neanthes.succinea, all_comm$Mulinia.lateralis, all_comm$Tellina.versicolor, all_comm$Nemertean, 
                      all_comm$abundance, all_comm$richness, all_comm$Bdep_season, all_comm$Bobs_season)
colnames(env_vars) <-c("salinity", "site", "season", "percent_clay", "percent_LOI", "temperature", "Bobs", "Bdep", "D", "C.louisianum",
                       "N.succinea", "M.lateralis", "T.versicolor", "nemertean", "abundance", "richness", "Bdep_season", "Bobs_season")

#lms for Bobs with interactions terms
lm1 <-lm(Bobs ~salinity*season*percent_clay*percent_LOI, data=env_vars) #full interactive
lm2 <-lm(Bobs ~salinity+season+percent_clay+percent_LOI, data=env_vars) #all additive
lm3 <-lm(Bobs ~salinity*season*percent_LOI, data=env_vars)
lm4 <-lm(Bobs ~salinity*percent_LOI*season*D, data=env_vars) #R2=0.55
lm5 <-lm(Bobs~salinity*season, data=env_vars)
lm6 <-lm(Bobs ~salinity*percent_LOI*D, data=env_vars) #R2=0.41
lm7 <-lm(Bobs~salinity*season+percent_LOI, data=env_vars)
lm8 <-lm(Bobs~salinity*season+percent_clay, data=env_vars)
lm9 <-lm(Bobs ~season*C.louisianum*salinity, data=env_vars) #R2=0.32
lm10 <-lm(Bobs ~season*C.louisianum*salinity+percent_LOI, data=env_vars) #R2=0.46
lm11 <-lm(Bobs ~salinity*season*percent_LOI, data=env_vars)#R2=0.469
lm12 <-lm(Bobs~salinity*season*D, data=env_vars)
lm13 <-lm(Bobs~salinity*season+D, data=env_vars)
lm14 <-lm(Bobs~salinity*season+percent_LOI+D, data=env_vars)

#additive lms only for Bobs:
lm40 <-lm(Bobs ~salinity, data=env_vars) 
lm41 <-lm(Bobs ~salinity+season, data=env_vars) 
lm42 <-lm(Bobs ~salinity+season+percent_LOI, data=env_vars) 
lm43 <-lm(Bobs ~percent_LOI+season+salinity, data=env_vars) 
lm44 <-lm(Bobs ~salinity+season+percent_LOI+percent_clay, data=env_vars)
lm45 <-lm(Bobs ~salinity+season+percent_LOI+percent_clay+D, data=env_vars)
lm46 <-lm(Bobs ~percent_LOI, data=env_vars) 
lm47 <-lm(Bobs ~percent_LOI+salinity, data=env_vars) 
lm48 <-lm(Bobs ~percent_LOI+D, data=env_vars)
lm49 <-lm(Bobs ~salinity+season+percent_LOI+percent_clay+temperature+D, data=env_vars) 
lm50 <-lm(Bobs ~salinity*percent_LOI, data=env_vars) 

#additive lms only AIC for Bobs
Bobs_additive_AIC <-AICc(lm40, lm41, lm42, lm43, lm44, lm45, lm46, lm47, lm48, lm49, lm50)
Bobs_additive_AIC

#best additive model for Bobs
lm43 <-lm(Bobs ~percent_LOI+season+salinity, data=env_vars)#R2=0.4694
summary(lm43)

#use AICc for small samples - model selection - Bobs
Bobs_AICc <-AICc(lm1, lm2, lm3, lm4, lm5, lm6, lm7, lm8, lm9,lm10, lm11, lm12, lm13, lm14)
Bobs_AICc

#best model for Bobs:
lm1 <-lm(Bobs ~salinity*season*percent_clay*percent_LOI, data=env_vars) #full interactive, R2=0.6259,p<2.2^-16
summary(lm1)

#### lm of env. parameters with departure from null model as response variable
lm15 <-lm(Bdep ~salinity*season*percent_clay*percent_LOI, data=env_vars) #full interactive
lm16 <-lm(Bdep ~salinity+season+percent_clay+percent_LOI, data=env_vars) #all additive
lm17 <-lm(Bdep ~salinity*season*percent_LOI, data=env_vars)
lm18 <-lm(Bdep ~salinity*percent_LOI*season*D, data=env_vars) 
lm19 <-lm(Bdep~salinity*season, data=env_vars)
lm20 <-lm(Bdep ~salinity*percent_LOI*D, data=env_vars) 
lm21 <-lm(Bdep~salinity*season+percent_LOI, data=env_vars)
lm22 <-lm(Bdep~salinity*season+percent_clay, data=env_vars)
lm23 <-lm(Bdep ~season*C.louisianum*salinity, data=env_vars) 
lm24 <-lm(Bdep ~season*C.louisianum*salinity+percent_LOI, data=env_vars) 
lm25 <-lm(Bdep ~salinity*season*percent_LOI, data=env_vars)
lm26 <-lm(Bdep~salinity*season*D, data=env_vars)
lm27 <-lm(Bdep~salinity*season+D, data=env_vars)
lm28 <-lm(Bdep~salinity*season+percent_LOI+D, data=env_vars)

#additive lms only
lm29 <-lm(Bdep ~percent_LOI+season+salinity+percent_clay, data=env_vars) 
lm30 <-lm(Bdep ~percent_LOI+salinity, data=env_vars) #
lm31 <-lm(Bdep ~salinity+season+percent_LOI, data=env_vars)
lm32 <-lm(Bdep ~percent_LOI+percent_clay, data=env_vars) 
lm33 <-lm(Bdep ~salinity+season+percent_LOI+percent_clay, data=env_vars) 
lm34 <-lm(Bdep ~salinity+season+percent_LOI+percent_clay+D, data=env_vars)
lm35 <-lm(Bdep ~salinity+percent_LOI+D, data=env_vars) 
lm36 <-lm(Bdep ~salinity+season+percent_LOI+percent_clay+temperature, data=env_vars) 
lm37 <-lm(Bdep ~percent_LOI*salinity, data=env_vars)
lm38 <-lm(Bdep ~salinity+season+percent_LOI+percent_clay+temperature+D, data=env_vars) 
lm39 <-lm(Bdep ~percent_LOI+season, data=env_vars) 

#additive AICs with interaction between salinity and season
Bdep_additive_AIC <-AICc(lm29, lm30, lm31, lm32, lm33, lm34, lm35, lm36, lm37, lm38, lm39)
Bdep_additive_AIC

#use AICc for small samples - model selection- Bdep
Bdep_AICc <-AICc(lm15, lm16, lm17, lm18, lm19, lm20, lm21, lm22, lm23, lm24, lm25, lm26, lm27, lm28)
Bdep_AICc

#best model for Bdep is fully interactive one
lm15 <-lm(Bdep ~salinity*season*percent_clay*percent_LOI, data=env_vars) #full interactive, R2=0.6048, p<2.2e-16
summary(lm15)

#best additive model for Bdep:
lm29 <-lm(Bdep ~percent_LOI+season+salinity+percent_clay, data=env_vars) #R2=0.5213

##when accounting for stochasticity, the LMs don't change between Bobs and Bdeps
##now separate out ONLY the "stochastic" sites - site that were within 2SD of the null model
##sites, 3,4,8,9,10,15
#making sure that the conclusions do not change if we're only considering stochastic sites

#subset only select sites
library(dplyr)
stoch_matrix <-filter(all_comm, all_comm$Site==3|all_comm$Site==4|all_comm$Site==8|all_comm$Site==9|all_comm$Site==10|all_comm$Site==15)

#lm of Bobs for stochastic sites only
lma <-lm(Bobs~percent_LOI, data=stoch_matrix)
lmb <-lm(Bobs~percent_LOI+Season, data=stoch_matrix)
lmc <-lm(Bobs~percent_LOI+Salinity, data=stoch_matrix)
lmd <-lm(Bobs ~percent_LOI*Salinity, data=stoch_matrix)
lme <-lm(Bobs~percent_LOI+Salinity+percent_clay, data=stoch_matrix)
lmf <-lm(Bobs~percent_LOI+Salinity+percent_clay+Season, data=stoch_matrix)
lmg<-lm(Bobs~Salinity+percent_clay*percent_LOI, data=stoch_matrix)
lmh <-lm(Bobs ~percent_LOI+Salinity+Corophium.louisianum, data=stoch_matrix)
lmi <-lm(Bobs ~percent_LOI*percent_clay*Salinity, data=stoch_matrix)


#AICc of Bobs models - stoch only
Bobs_stoch_aic <-AICc <-AICc(lma, lmb, lmc, lmd, lme, lmf, lmg, lmh, lmi)
Bobs_stoch_aic

#best model overall
lmi <-lm(Bobs ~percent_LOI*percent_clay*Salinity, data=stoch_matrix) #R2=0.7493

#best additive model + partial regression plot
library(car)
lmc <-lm(Bobs ~percent_LOI+Salinity, data=stoch_matrix) #R2=0.3396
plot(Bobs~percent_LOI+Salinity, data=stoch_matrix)
summary(lmc)
par(mfrow=c(2,2))
avPlots(lmc, id=FALSE, main=" ")

par(mfrow=c(1,2))
#plot of only percent_LOI lm; R2=0.2893
ggplot(data=stoch_matrix,aes(x=percent_LOI, y=Bobs))+
  geom_point(alpha=0.5, pch=1)+ geom_jitter()+
  geom_smooth(method='lm', se=FALSE)+
  theme_classic()

#plot of only percent_LOI against Bdep; R2 =0.3359
lm <-lm(Bdep ~percent_LOI, data=stoch_matrix)
summary(lm)
ggplot(data=stoch_matrix,aes(x=percent_LOI, y=Bdep))+
  geom_point(alpha=0.5, pch=1)+ geom_jitter()+
  geom_smooth(method='lm', se=FALSE)+
  theme_classic()

#lm of Bdep for stochastic sites only
lmj <-lm(Bdep~percent_LOI, data=stoch_matrix)
lmk <-lm(Bdep~percent_LOI+Season, data=stoch_matrix)
lml <-lm(Bdep~percent_LOI+Salinity, data=stoch_matrix)
lmm <-lm(Bdep ~percent_LOI*Salinity, data=stoch_matrix)
lmn <-lm(Bdep~percent_LOI+Salinity+percent_clay, data=stoch_matrix)
lmo <-lm(Bdep~percent_LOI+Salinity+percent_clay+Season, data=stoch_matrix)
lmp<-lm(Bdep~Salinity+percent_clay*percent_LOI, data=stoch_matrix)
lmq <-lm(Bdep ~percent_LOI+Salinity+Corophium.louisianum, data=stoch_matrix)
lmr <-lm(Bdep ~percent_LOI*percent_clay*Salinity, data=stoch_matrix)

#AICc of Bdep models - stoch only
Bdeps_stoch_aic <-AICc(lmj, lmk, lml, lmm, lmn, lmo, lmp, lmq, lmr)
Bdeps_stoch_aic

#best model overall
lmr <-lm(Bdep ~percent_LOI*percent_clay*Salinity, data=stoch_matrix)#R2=0.5799

#best additive model
lml <-lm(Bdep~percent_LOI+Salinity, data=stoch_matrix)#0.3519

#partial regression plot
lm <-lm(Bobs ~percent_LOI, data=stoch_matrix) #R2=0.3396
plot(Bobs~percent_LOI, data=stoch_matrix)
summary(lm)
par(mfrow=c(2,2))
plot(Bobs ~percent_LOI, data=stoch_matrix)
abline()
avPlots(lm, id=FALSE, main=" ", grid=FALSE)

##accounting for spatial autocorrelation
#calculate distances between every site pair. Use Bobs for each pairwise site observation as response variable 
#distances between each point (lat/lon) as predictor variable. 

library(geosphere)

site_loc #dataset
site_dist <-site_loc[,1:3]
site_dist #subset with only site and lat/lon 
distm(c(lon1, lat1), c(lon2, lat2), fun = distHaversine)
site1.2 <-distm(c(site_dist[1,2], site_dist[1,3]), c(site_dist[2,2], site_dist[2,3]), fun = distHaversine) #use this to calculate distances between each site
site1.3 <-distm(c(site_dist[1,2], site_dist[1,3]), c(site_dist[3,2], site_dist[3,3]), fun = distHaversine)
site1.4 <-distm(c(site_dist[1,2], site_dist[1,3]), c(site_dist[4,2], site_dist[4,3]), fun = distHaversine)
site1.5 <-distm(c(site_dist[1,2], site_dist[1,3]), c(site_dist[5,2], site_dist[5,3]), fun = distHaversine)
site1.6 <-distm(c(site_dist[1,2], site_dist[1,3]), c(site_dist[6,2], site_dist[6,3]), fun = distHaversine)
site1.7 <-distm(c(site_dist[1,2], site_dist[1,3]), c(site_dist[7,2], site_dist[7,3]), fun = distHaversine)
site1.8 <-distm(c(site_dist[1,2], site_dist[1,3]), c(site_dist[8,2], site_dist[8,3]), fun = distHaversine)
site1.9 <-distm(c(site_dist[1,2], site_dist[1,3]), c(site_dist[9,2], site_dist[9,3]), fun = distHaversine)
site1.10 <-distm(c(site_dist[1,2], site_dist[1,3]), c(site_dist[10,2], site_dist[10,3]), fun = distHaversine)
site1.11 <-distm(c(site_dist[1,2], site_dist[1,3]), c(site_dist[11,2], site_dist[11,3]), fun = distHaversine)
site1.12 <-distm(c(site_dist[1,2], site_dist[1,3]), c(site_dist[12,2], site_dist[12,3]), fun = distHaversine)
site1.13 <-distm(c(site_dist[1,2], site_dist[1,3]), c(site_dist[13,2], site_dist[13,3]), fun = distHaversine)
site1.14 <-distm(c(site_dist[1,2], site_dist[1,3]), c(site_dist[14,2], site_dist[14,3]), fun = distHaversine)
site1.15 <-distm(c(site_dist[1,2], site_dist[1,3]), c(site_dist[15,2], site_dist[15,3]), fun = distHaversine)
site2.3 <-distm(c(site_dist[2,2], site_dist[2,3]), c(site_dist[3,2], site_dist[3,3]), fun = distHaversine)
site2.4 <-distm(c(site_dist[2,2], site_dist[2,3]), c(site_dist[4,2], site_dist[4,3]), fun = distHaversine)
site2.5 <-distm(c(site_dist[2,2], site_dist[2,3]), c(site_dist[5,2], site_dist[5,3]), fun = distHaversine)
site2.6 <-distm(c(site_dist[2,2], site_dist[2,3]), c(site_dist[6,2], site_dist[6,3]), fun = distHaversine)
site2.7 <-distm(c(site_dist[2,2], site_dist[2,3]), c(site_dist[7,2], site_dist[7,3]), fun = distHaversine)
site2.8 <-distm(c(site_dist[2,2], site_dist[2,3]), c(site_dist[8,2], site_dist[8,3]), fun = distHaversine)
site2.9 <-distm(c(site_dist[2,2], site_dist[2,3]), c(site_dist[9,2], site_dist[9,3]), fun = distHaversine)
site2.10 <-distm(c(site_dist[2,2], site_dist[2,3]), c(site_dist[10,2], site_dist[10,3]), fun = distHaversine)
site2.11 <-distm(c(site_dist[2,2], site_dist[2,3]), c(site_dist[11,2], site_dist[11,3]), fun = distHaversine)
site2.12 <-distm(c(site_dist[2,2], site_dist[2,3]), c(site_dist[12,2], site_dist[12,3]), fun = distHaversine)
site2.13 <-distm(c(site_dist[2,2], site_dist[2,3]), c(site_dist[13,2], site_dist[13,3]), fun = distHaversine)
site2.14 <-distm(c(site_dist[2,2], site_dist[2,3]), c(site_dist[14,2], site_dist[14,3]), fun = distHaversine)
site2.15 <-distm(c(site_dist[2,2], site_dist[2,3]), c(site_dist[15,2], site_dist[15,3]), fun = distHaversine)
site3.4 <-distm(c(site_dist[3,2], site_dist[3,3]), c(site_dist[4,2], site_dist[4,3]), fun = distHaversine)
site3.5 <-distm(c(site_dist[3,2], site_dist[3,3]), c(site_dist[5,2], site_dist[5,3]), fun = distHaversine)
site3.6 <-distm(c(site_dist[3,2], site_dist[3,3]), c(site_dist[6,2], site_dist[6,3]), fun = distHaversine)
site3.7 <-distm(c(site_dist[3,2], site_dist[3,3]), c(site_dist[7,2], site_dist[7,3]), fun = distHaversine)
site3.8 <-distm(c(site_dist[3,2], site_dist[3,3]), c(site_dist[8,2], site_dist[8,3]), fun = distHaversine)
site3.9 <-distm(c(site_dist[3,2], site_dist[3,3]), c(site_dist[9,2], site_dist[9,3]), fun = distHaversine)
site3.10 <-distm(c(site_dist[3,2], site_dist[3,3]), c(site_dist[10,2], site_dist[10,3]), fun = distHaversine)
site3.11 <-distm(c(site_dist[3,2], site_dist[3,3]), c(site_dist[11,2], site_dist[11,3]), fun = distHaversine)
site3.12 <-distm(c(site_dist[3,2], site_dist[3,3]), c(site_dist[12,2], site_dist[12,3]), fun = distHaversine)
site3.13 <-distm(c(site_dist[3,2], site_dist[3,3]), c(site_dist[13,2], site_dist[13,3]), fun = distHaversine)
site3.14 <-distm(c(site_dist[3,2], site_dist[3,3]), c(site_dist[14,2], site_dist[14,3]), fun = distHaversine)
site3.15 <-distm(c(site_dist[3,2], site_dist[3,3]), c(site_dist[15,2], site_dist[15,3]), fun = distHaversine)
site4.5 <-distm(c(site_dist[4,2], site_dist[4,3]), c(site_dist[5,2], site_dist[5,3]), fun = distHaversine)
site4.6 <-distm(c(site_dist[4,2], site_dist[4,3]), c(site_dist[6,2], site_dist[6,3]), fun = distHaversine)
site4.7 <-distm(c(site_dist[4,2], site_dist[4,3]), c(site_dist[7,2], site_dist[7,3]), fun = distHaversine)
site4.8 <-distm(c(site_dist[4,2], site_dist[4,3]), c(site_dist[8,2], site_dist[8,3]), fun = distHaversine)
site4.9 <-distm(c(site_dist[4,2], site_dist[4,3]), c(site_dist[9,2], site_dist[9,3]), fun = distHaversine)
site4.10 <-distm(c(site_dist[4,2], site_dist[4,3]), c(site_dist[10,2], site_dist[10,3]), fun = distHaversine)
site4.11 <-distm(c(site_dist[4,2], site_dist[4,3]), c(site_dist[11,2], site_dist[11,3]), fun = distHaversine)
site4.12 <-distm(c(site_dist[4,2], site_dist[4,3]), c(site_dist[12,2], site_dist[12,3]), fun = distHaversine)
site4.13 <-distm(c(site_dist[4,2], site_dist[4,3]), c(site_dist[13,2], site_dist[13,3]), fun = distHaversine)
site4.14 <-distm(c(site_dist[4,2], site_dist[4,3]), c(site_dist[14,2], site_dist[14,3]), fun = distHaversine)
site4.15 <-distm(c(site_dist[4,2], site_dist[4,3]), c(site_dist[15,2], site_dist[15,3]), fun = distHaversine)
site5.6 <-distm(c(site_dist[5,2], site_dist[5,3]), c(site_dist[6,2], site_dist[6,3]), fun = distHaversine)
site5.7 <-distm(c(site_dist[5,2], site_dist[5,3]), c(site_dist[7,2], site_dist[7,3]), fun = distHaversine)
site5.8 <-distm(c(site_dist[5,2], site_dist[5,3]), c(site_dist[8,2], site_dist[8,3]), fun = distHaversine)
site5.9 <-distm(c(site_dist[5,2], site_dist[5,3]), c(site_dist[9,2], site_dist[9,3]), fun = distHaversine)
site5.10 <-distm(c(site_dist[5,2], site_dist[5,3]), c(site_dist[10,2], site_dist[10,3]), fun = distHaversine)
site5.11 <-distm(c(site_dist[5,2], site_dist[5,3]), c(site_dist[11,2], site_dist[11,3]), fun = distHaversine)
site5.12 <-distm(c(site_dist[5,2], site_dist[5,3]), c(site_dist[12,2], site_dist[12,3]), fun = distHaversine)
site5.13 <-distm(c(site_dist[5,2], site_dist[5,3]), c(site_dist[13,2], site_dist[13,3]), fun = distHaversine)
site5.14 <-distm(c(site_dist[5,2], site_dist[5,3]), c(site_dist[14,2], site_dist[14,3]), fun = distHaversine)
site5.15 <-distm(c(site_dist[5,2], site_dist[5,3]), c(site_dist[15,2], site_dist[15,3]), fun = distHaversine)
site6.7 <-distm(c(site_dist[6,2], site_dist[6,3]), c(site_dist[7,2], site_dist[7,3]), fun = distHaversine)
site6.8 <-distm(c(site_dist[6,2], site_dist[6,3]), c(site_dist[8,2], site_dist[8,3]), fun = distHaversine)
site6.9 <-distm(c(site_dist[6,2], site_dist[6,3]), c(site_dist[9,2], site_dist[9,3]), fun = distHaversine)
site6.10 <-distm(c(site_dist[6,2], site_dist[6,3]), c(site_dist[10,2], site_dist[10,3]), fun = distHaversine)
site6.11 <-distm(c(site_dist[6,2], site_dist[6,3]), c(site_dist[11,2], site_dist[11,3]), fun = distHaversine)
site6.12 <-distm(c(site_dist[6,2], site_dist[6,3]), c(site_dist[12,2], site_dist[12,3]), fun = distHaversine)
site6.13 <-distm(c(site_dist[6,2], site_dist[6,3]), c(site_dist[13,2], site_dist[13,3]), fun = distHaversine)
site6.14 <-distm(c(site_dist[6,2], site_dist[6,3]), c(site_dist[14,2], site_dist[14,3]), fun = distHaversine)
site6.15 <-distm(c(site_dist[6,2], site_dist[6,3]), c(site_dist[15,2], site_dist[15,3]), fun = distHaversine)
site7.8 <-distm(c(site_dist[7,2], site_dist[7,3]), c(site_dist[8,2], site_dist[8,3]), fun = distHaversine)
site7.9 <-distm(c(site_dist[7,2], site_dist[7,3]), c(site_dist[9,2], site_dist[9,3]), fun = distHaversine)
site7.10 <-distm(c(site_dist[7,2], site_dist[7,3]), c(site_dist[10,2], site_dist[10,3]), fun = distHaversine)
site7.11 <-distm(c(site_dist[7,2], site_dist[7,3]), c(site_dist[11,2], site_dist[11,3]), fun = distHaversine)
site7.12 <-distm(c(site_dist[7,2], site_dist[7,3]), c(site_dist[12,2], site_dist[12,3]), fun = distHaversine)
site7.13 <-distm(c(site_dist[7,2], site_dist[7,3]), c(site_dist[13,2], site_dist[13,3]), fun = distHaversine)
site7.14 <-distm(c(site_dist[7,2], site_dist[7,3]), c(site_dist[14,2], site_dist[14,3]), fun = distHaversine)
site7.15 <-distm(c(site_dist[7,2], site_dist[7,3]), c(site_dist[15,2], site_dist[15,3]), fun = distHaversine)
site8.9 <-distm(c(site_dist[8,2], site_dist[8,3]), c(site_dist[9,2], site_dist[9,3]), fun = distHaversine)
site8.10 <-distm(c(site_dist[8,2], site_dist[8,3]), c(site_dist[10,2], site_dist[10,3]), fun = distHaversine)
site8.11 <-distm(c(site_dist[8,2], site_dist[8,3]), c(site_dist[11,2], site_dist[11,3]), fun = distHaversine)
site8.12 <-distm(c(site_dist[8,2], site_dist[8,3]), c(site_dist[12,2], site_dist[12,3]), fun = distHaversine)
site8.13 <-distm(c(site_dist[8,2], site_dist[8,3]), c(site_dist[13,2], site_dist[13,3]), fun = distHaversine)
site8.14 <-distm(c(site_dist[8,2], site_dist[8,3]), c(site_dist[14,2], site_dist[14,3]), fun = distHaversine)
site8.15 <-distm(c(site_dist[8,2], site_dist[8,3]), c(site_dist[15,2], site_dist[15,3]), fun = distHaversine)
site9.10 <-distm(c(site_dist[9,2], site_dist[9,3]), c(site_dist[10,2], site_dist[10,3]), fun = distHaversine)
site9.11 <-distm(c(site_dist[9,2], site_dist[9,3]), c(site_dist[11,2], site_dist[11,3]), fun = distHaversine)
site9.12 <-distm(c(site_dist[9,2], site_dist[9,3]), c(site_dist[12,2], site_dist[12,3]), fun = distHaversine)
site9.13 <-distm(c(site_dist[9,2], site_dist[9,3]), c(site_dist[13,2], site_dist[13,3]), fun = distHaversine)
site9.14 <-distm(c(site_dist[9,2], site_dist[9,3]), c(site_dist[14,2], site_dist[14,3]), fun = distHaversine)
site9.15 <-distm(c(site_dist[9,2], site_dist[9,3]), c(site_dist[15,2], site_dist[15,3]), fun = distHaversine)
site10.11 <-distm(c(site_dist[10,2], site_dist[10,3]), c(site_dist[11,2], site_dist[11,3]), fun = distHaversine)
site10.12 <-distm(c(site_dist[10,2], site_dist[10,3]), c(site_dist[12,2], site_dist[12,3]), fun = distHaversine)
site10.13 <-distm(c(site_dist[10,2], site_dist[10,3]), c(site_dist[13,2], site_dist[13,3]), fun = distHaversine)
site10.14 <-distm(c(site_dist[10,2], site_dist[10,3]), c(site_dist[14,2], site_dist[14,3]), fun = distHaversine)
site10.15 <-distm(c(site_dist[10,2], site_dist[10,3]), c(site_dist[15,2], site_dist[15,3]), fun = distHaversine)
site11.12 <-distm(c(site_dist[11,2], site_dist[11,3]), c(site_dist[12,2], site_dist[12,3]), fun = distHaversine)
site11.13 <-distm(c(site_dist[11,2], site_dist[11,3]), c(site_dist[13,2], site_dist[13,3]), fun = distHaversine)
site11.14 <-distm(c(site_dist[11,2], site_dist[11,3]), c(site_dist[14,2], site_dist[14,3]), fun = distHaversine)
site11.15 <-distm(c(site_dist[11,2], site_dist[11,3]), c(site_dist[15,2], site_dist[15,3]), fun = distHaversine)
site12.13 <-distm(c(site_dist[12,2], site_dist[12,3]), c(site_dist[13,2], site_dist[13,3]), fun = distHaversine)
site12.14 <-distm(c(site_dist[12,2], site_dist[12,3]), c(site_dist[14,2], site_dist[14,3]), fun = distHaversine)
site12.15 <-distm(c(site_dist[12,2], site_dist[12,3]), c(site_dist[15,2], site_dist[15,3]), fun = distHaversine)
site13.14 <-distm(c(site_dist[13,2], site_dist[13,3]), c(site_dist[14,2], site_dist[14,3]), fun = distHaversine)
site13.15 <-distm(c(site_dist[13,2], site_dist[13,3]), c(site_dist[15,2], site_dist[15,3]), fun = distHaversine)
site14.15 <-distm(c(site_dist[14,2], site_dist[14,3]), c(site_dist[15,2], site_dist[15,3]), fun = distHaversine)

#make a vector- 1:120
dist <-c(0,site1.2, site1.3, site1.4, site1.5, site1.6, site1.7, site1.8, site1.9, site1.10, site1.11, site1.12, site1.13, 
         site1.14, site1.15, 0, site2.3, site2.4, site2.5, site2.6, site2.7, site2.8, site2.9, site2.10, site2.11, site2.12, 
         site2.13, site2.14, site2.15, 0,site3.4, site3.5, site3.6, site3.7, site3.8, site3.9, site3.10, site3.11, site3.12,
         site3.13, site3.14, site3.15, 0, site4.5, site4.6, site4.7, site4.8, site4.9, site4.10, site4.11, site4.12, site4.13,
         site4.14, site4.15, 0, site5.6, site5.7, site5.8, site5.9, site5.10, site5.11, site5.12, site5.13, site5.14, site5.15, 
         0, site6.7, site6.8, site6.9, site6.10, site6.11, site6.12, site6.13, site6.14, site6.15, 0, site7.8, site7.9, site7.10,
         site7.11, site7.12, site7.13, site7.14, site7.15, 0, site8.9, site8.10, site8.11, site8.12, site8.13, site8.14, site8.15, 
         0, site9.10, site9.11, site9.12, site9.13, site9.14, site9.15, 0, site10.11, site10.12, site10.13, site10.14, site10.15,
         0, site11.12, site11.13,site11.14, site11.15, 0, site12.13, site12.14, site12.15, 0, site13.14, site13.15, 0, site14.15, 0)

mean.bray <-meandist(braycurtis, species_env_matrix$all_comm.Site) #from above

mean.bray.new<- mean.bray                           # Duplicate matrix
mean.bray.new[lower.tri(mean.bray.new)] <- NA       # Change diagonal of matrix
mean.bray.new
mean.bray.upper=c(mean.bray.new) #transform BC matrix into vector
mean.bray.upper=mean.bray.upper[!is.na(mean.bray.upper)]    #remove NAs from vector

#make new dataset
#there is no relationship between distance and Bobs - therefore no spatial autocorrelation. 
spatial.autocor <-data.frame(mean.bray.upper, dist)
colnames(spatial.autocor) <-c("Bobs", "dist")
plot(Bobs ~dist, data=spatial.autocor) #no obvious relationship based on visual
lm.spatial <-lm(Bobs~dist, data=spatial.autocor)
summary(lm.spatial)
plot(lm.spatial)



