# Here we will look at two ways of analysing a 2 x 2 design which has:
# 2 levels of word_type (1N，2B)
# 2 levels of familiarity(3L, 4H)
# We will start off with the same steps:

library(MASS)
library(lme4)
library(lattice)
library(plyr)
library(lmerTest)
library(Matrix)

rm(list=ls()) #clean

# Open "2x2 Data"
# 读取.csv的文件
datafile = read.csv(file.choose(), sep = ",", dec = ".")
colnames(datafile)
###################################################################################################################################################
datafile = Lex_RT_Clean_22
###################################################################################################################################################
# 第一列是被试，Subject被试，condition是条件，stim是刺激句子，crit删除的是正负三个标准差。
col.subject = 1      # specify which column your participant number is in
col.condition = 3   # specify which column contains your condition number (used for excluding outliers)
col.stim = 2       # specify which column your item number is in
crit = 3          # setting the criterium value you want, usually 3 or 2.5 standard deviations

# choose dependent variable - here I have several measures - I will use total reading times on word n
# Measure是要分析的指标，因变量。
measure = "Probe1.RT" 

# add to dataframe
# 所有的depvar就是因变量了。赋值的过程。
datafile$depvar = datafile[,measure]

# Work out which colums the fixed and random factors are in
# 继续命名的过程，被试命名成pp。
datafile$pp = datafile[,col.subject]
datafile$condition = datafile[,col.condition]
datafile$stim = datafile[,col.stim]

# Let's make sure all the variables are from the correct class
datafile$depvar = as.numeric(datafile$depvar)
datafile$pp = as.factor(datafile$pp)
datafile$stim = as.factor(datafile$stim)
datafile$condition = as.factor(datafile$condition)
str(datafile) # 数值型的命名。或者因素型的命名。

###################################################################################################################################################
# 删标准差的过程，正负三个标准差。
# Do we want to remove outliers
# If yes ...
nrcolumns = ncol(datafile)   # We'll be using this for cleanup later
# let's first have a look
head(datafile)
# make a new datafile with fixation durations > 0 ms
between1= datafile[datafile$depvar > 0, ]
head(between1)
# make a matrix with the means per subject and per condition for the dependent variable
mean.matrix = tapply(between1$depvar, list(between1$pp), mean, na.rm = T)
mean.matrix
# make a matrix with the standard deviations per subject and per condition for the dependent variable
sd.matrix = tapply(between1$depvar, list(between1$pp), sd, na.rm = T)
sd.matrix
# add this data to the actual dataframe
for(i in 1:nrow(between1)){between1$mean.sc[i] = mean.matrix[between1$pp[i]]}
for(i in 1:nrow(between1)){between1$sd.sc[i] = sd.matrix[between1$pp[i]]}
# calculate z-scores
for(i in 1:nrow(between1)) {between1$zscore[i] = (between1$depvar[i] - between1$mean.sc[i])/between1$sd.sc[i]}
# assign a zero value to cells with only 1 observation
between1$zscore[is.na(between1$zscore)] = 0
# make new matrix with only standard deviations below a certain criterium of z-score (in absolute value)
result = between1[abs(between1$zscore) < 3,]
# you might want to clean up this data file by removing the added columns
result =  result[,1:nrcolumns]
datafile <-result
nrow(datafile)

###################################################################################################################################################
# 中心化 launch_site 连续变量作为协变量，加入分析，还有很多可以这样分析，比如词频。
# centering  launch_site

#datafile$claunch_site = datafile$launch_site -mean(datafile$launch_site)
#head(datafile)
###################################################################################################################################################
# 平均数 标准差 标准误
# Code for calculating means between per pp and per condition
mean.tt = tapply(datafile$depvar, list(datafile$pp, datafile$condition), mean, na.rm = T)
mean.tt
sd.tt = tapply(datafile$depvar, list(datafile$pp, datafile$condition), sd, na.rm = T)
sd.tt
# then calculate the grand mean by condition
grand.mean=apply(mean.tt, 2, mean, na.rm = T)
grand.mean
grand.sd=apply(sd.tt,2,mean, na.rm = T)
grand.sd 
grand.se=grand.sd/sqrt(32) # this number is the total number of participants from your data, so it is changable.
grand.se

###################################################################################################################################################
# LOG之后的数据更线性，更符合预期，所以一般都要对注视时间上的数据进行log转换
# Do we need to log transform?
qqnorm(datafile$depvar) # log之前的图
qqnorm(log(datafile$depvar)) # log之后的图 看一下
# If yes ...
datafile$depvar = log(datafile$depvar) #这就是在log了
# How does it look at an individual level
qqmath(~depvar|pp, data = datafile) # log之后每个被试的图，如果趋势差的特别大，考虑删掉那个被试。

###################################################################################################################################################
# 相邻的比较。
contrasts(datafile$familiarity) <- contr.sdif(2)
contrasts(datafile$word_type) <- contr.sdif(2)

#full model, main effects
depvar.lmer1 = lmer(depvar ~ familiarity * word_type+ (1 + familiarity * word_type|pp) + (1 + familiarity * word_type|stim), datafile)
summary(depvar.lmer1, corr = FALSE)
#1代表的是截距，线性混合模型，贡献率是多少，preview*complex是对斜率的贡献（交互）

depvar.lmer2 = lmer(depvar ~ familiarity * word_type + (1 + familiarity + word_type|pp) + (1 + familiarity + word_type|stim), datafile)
summary(depvar.lmer2, corr = FALSE)
#*变成+，单独的贡献率是多少。

depvar.lmer3 = lmer(depvar ~ familiarity * word_type+ (1 + familiarity + word_type|pp) + (1|stim), datafile)
summary(depvar.lmer3, corr = FALSE)
#只考虑一个因素的情况下，越来越简单。

depvar.lmer4 = lmer(depvar ~  familiarity * word_type + (1|pp) + (1|stim), datafile)
summary(depvar.lmer4, corr = FALSE)
#最简单的模型。

depvar.lmer5 = lmer(depvar ~  familiarity * word_type + (1|pp) , datafile)
summary(depvar.lmer5, corr = FALSE)
#如果最简单的模型还跑不出来，就把stim也删掉，肯定能跑出来了。

###################################################################################################################################################
# 简单效应分析
depvar.lmer1 = lmer(depvar ~ familiarity+ (1+ familiarity*word_type|pp) + (1+ familiarity*word_type|stim),subset(datafile,word_type=='1N'))
summary(depvar.lmer1, corr = FALSE)

depvar.lmer2 = lmer(depvar ~ familiarity+ (1+ familiarity+word_type|pp) + (1+ familiarity+word_type|stim),subset(datafile,word_type=='1N'))
summary(depvar.lmer2, corr = FALSE)

depvar.lmer3 = lmer(depvar ~ familiarity+ (1+ familiarity+word_type|pp) + (1|stim),subset(datafile,word_type=='1N'))
summary(depvar.lmer3, corr = FALSE)

depvar.lmer4 = lmer(depvar ~ familiarity+ (1|pp) + (1|stim),subset(datafile,word_type=='1N'))
summary(depvar.lmer4, corr = FALSE)

depvar.lmer5 = lmer(depvar ~ familiarity+ (1|pp),subset(datafile,word_type=='1N'))
summary(depvar.lmer5, corr = FALSE)
# 同样情况下分析word_type=='2B'
# familiarity=='3L'
# familiarity=='4H'

###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
# sKip / IA_REGRESSION_OUT
# 3*2设计需要用到下述命令
contrast.matrix <- matrix(c(   
  -1, 0, +1,  # ID vs TP
  0, -1, +1)  # ID vs JC
  , 3, 2,      # number of conditions, number of comparisons
  dimnames=list(c("1P", "2U", "3N"), # label conditions
                c("1P vs 3N", " 2U vs 3N"))) # label contrasts

# Assign contrast matrix
inv.contrast.matrix <- matrix(t(ginv(contrast.matrix)), 3, 2,
                              dimnames=list(c("1P", "2U", "3N"), 
                                            c("1P vs 3N", " 2U vs 3N")))
(contrasts(datafile$emotion) <- fractions(inv.contrast.matrix))

###################################################################################################################################################
install.packages("languageR")

# The code for skipping
library(languageR)
library(lme4)
library(lattice)
library(MASS)
library(plyr)
library(lmerTest)

rm(list=ls()) # clean

# Open "2x2Data"
datafile = read.csv(file.choose(), sep = ",", dec = ".")
colnames(datafile)
###################################################################################################################################################
datafile = Lex_RT_Clean_22
###################################################################################################################################################
col.subject = 1      # specify which column your participant number is in
col.condition = 3  # specify which column contains your condition number (used for excluding outliers)
col.stim = 2     # specify which column your item number is in

# choose dependent variable
measure = "Probe1.ACC"

# add to dataframe
datafile$depvar = datafile[,measure]

# Work out which colums the fixed and random factors are in
datafile$pp = datafile[,col.subject]
datafile$condition = datafile[,col.condition]
datafile$stim = datafile[,col.stim]

# Let's make sure all the variables are from the correct class
datafile$depvar = as.numeric(datafile$depvar)
datafile$pp = as.factor(datafile$pp)
datafile$stim = as.factor(datafile$stim)
datafile$condition = as.factor(datafile$condition)

###################################################################################################################################################
# centering  launch_site
#datafile$claunch_site = datafile$launch_site -mean(datafile$launch_site)
#head(datafile)

###################################################################################################################################################
# Code for calculating means between per pp and per condition
mean.tt = tapply(datafile$depvar, list(datafile$pp, datafile$condition), mean, na.rm = T)
mean.tt
sd.tt = tapply(datafile$depvar, list(datafile$pp, datafile$condition), sd, na.rm = T)
sd.tt
# then calculate the grand mean by condition
grand.mean=apply(mean.tt, 2, mean, na.rm = T)
grand.mean
grand.sd=apply(sd.tt,2,mean, na.rm = T)
grand.sd 
grand.se=grand.sd/sqrt(32)  ## this number is the total number of participants from your data, so it is changable.
grand.se
# Add contrasts
contrasts(datafile$familiarity) <- contr.sdif(2)
contrasts(datafile$word_type) <- contr.sdif(2)

###################################################################################################################################################
#full model, main effects
depvar.lmer1 = glmer(depvar ~ familiarity * word_type+ (1 + familiarity * word_type|pp) + (1 + familiarity * word_type|stim), datafile, family = binomial)
summary(depvar.lmer1, corr = FALSE)
#1代表的是截距，线性混合模型，贡献率是多少，preview*complex是对斜率的贡献（交互）

depvar.lmer2 = glmer(depvar ~ familiarity * word_type + (1 + familiarity + word_type|pp) + (1 + familiarity + word_type|stim), datafile, family = binomial)
summary(depvar.lmer2, corr = FALSE)
#*变成+，单独的贡献率是多少。

depvar.lmer3 = glmer(depvar ~ familiarity * word_type+ (1 + familiarity + word_type|pp) + (1|stim), datafile, family = binomial)
summary(depvar.lmer3, corr = FALSE)
#只考虑一个因素的情况下，越来越简单。

depvar.lmer4 = glmer(depvar ~  familiarity * word_type + (1|pp) + (1|stim), datafile, family = binomial)
summary(depvar.lmer4, corr = FALSE)
#最简单的模型。

depvar.lmer5 = glmer(depvar ~  familiarity * word_type + (1|pp) , datafile, family = binomial)
summary(depvar.lmer5, corr = FALSE)
#如果最简单的模型还跑不出来，就把stim也删掉，肯定能跑出来了。

###################################################################################################################################################
# 简单效应分析
depvar.lmer1 = glmer(depvar ~ familiarity+ (1+ familiarity*word_type|pp) + (1+ familiarity*word_type|stim),subset(datafile,word_type=='1N'), family = binomial)
summary(depvar.lmer1, corr = FALSE)

depvar.lmer2 = glmer(depvar ~ familiarity+ (1+ familiarity+word_type|pp) + (1+ familiarity+word_type|stim),subset(datafile,word_type=='1N'), family = binomial)
summary(depvar.lmer2, corr = FALSE)

depvar.lmer3 = glmer(depvar ~ familiarity+ (1+ familiarity+word_type|pp) + (1|stim),subset(datafile,word_type=='1N'), family = binomial)
summary(depvar.lmer3, corr = FALSE)

depvar.lmer4 = glmer(depvar ~ familiarity+ (1|pp) + (1|stim),subset(datafile,word_type=='1N'), family = binomial)
summary(depvar.lmer4, corr = FALSE)

depvar.lmer4 = glmer(depvar ~ familiarity+ (1|pp) + (1|stim),subset(datafile,word_type=='2B'), family = binomial)
summary(depvar.lmer4, corr = FALSE)

depvar.lmer4 = glmer(depvar ~ word_type+ (1|pp) + (1|stim),subset(datafile,familiarity=='3L'), family = binomial)
summary(depvar.lmer4, corr = FALSE)

depvar.lmer4 = glmer(depvar ~ word_type+ (1|pp) + (1|stim),subset(datafile,familiarity=='4H'), family = binomial)
summary(depvar.lmer4, corr = FALSE)

depvar.lmer5 = glmer(depvar ~ familiarity+ (1|pp),subset(datafile,word_type=='1N'), family = binomial)
summary(depvar.lmer5, corr = FALSE)
# 同样情况下分析word_type=='2B'
# familiarity=='3L'
# familiarity=='4H'

###################################################################################################################################################
###################################################################################################################################################
# 原有命令
# 简单效应分析
m.1 = glmer(depvar ~ word_type + (1 |pp) + (1 |stim), subset(datafile,  familiarity=='3L'), family = binomial)
summary(m.1)

m.1 = glmer(depvar ~ word_type + (1 |pp) + (1 |stim), subset(datafile,  familiarity=='4H'), family = binomial)
summary(m.1)

m.1 = glmer(depvar ~ familiarity + (1 |pp) + (1 |stim), subset(datafile,  word_type=='1N'), family = binomial)
summary(m.1)

m.1 = glmer(depvar ~ familiarity + (1 |pp) + (1 |stim), subset(datafile,  word_type=='2B'), family = binomial)
summary(m.1)

# 2*2*3 简单简单效应分析
m.1 = glmer(depvar ~ word_type + (1 |pp) + (1 |stim), subset(datafile, nianling=='2O'& familiarity=='2P'), family = binomial)
summary(m.1)

###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
(table1 <- ddply(datafile, .(emotion, frequency, nianling), summarise, N=length(depvar[!is.na(depvar)]), M=mean(depvar, na.rm = TRUE), SD=sd(depvar, na.rm = TRUE), SE=SD/sqrt(N) ))
