## open file with data
#for this script, use HP > Box > New Computer > Auburn > Data > MammalMetaDataEdit
#the "meta_m.csv" data has better column titles and only the mammals. after class, more data will be added
setwd("C:/Users/HP/Box/New Computer/Auburn/Data") #set working directory in the Rstudio folder on my computer and box drive

meta = read.csv("meta_m.csv")  ##meta = read.csv(file.choose()) #use if you want to select file each time
meta  ##check data
##install metafor package if not already done
library(metafor) ##load metafor package

#add column for response ratio
"Will need to finalize the FST and COAL columns to make sure they are the correct value"
hist(meta$calcm)
hist(meta$m_rate)
#both calcm and m_rate are skewed to the right
#calcm -> final FST value to use. this was calculated as calcm ~ Nm ~ .25(1/Fst-1) AKA effective number of migrants
#m_rate -> final COAL value to use. this is the migration rate per generation, m, calculated with a coalescent approach
#RR -> response ratio
#RRlog -> log response ratio

#test for normality
shapiro.test(meta$calcm)
shapiro.test(meta$m_rate)
"Not normal, but that's not a surprise"

#visualize which ratio to go with -- Who's on top?
meta$ORR = meta$calcm/meta$m_rate
hist(meta$ORR) #skewed to the right

meta$ORR = meta$m_rate/meta$calcm
hist(meta$ORR) #also skewed to the right, but the scale is better
"m_rate is a top. calcm is a bottom"
#if <1, negative response. if = 1, no response. if >1, positive response

#check the log response ratio
meta$LOG = log(meta$m_rate+1) - log(meta$calcm) #add 1 for the low migration numbers
hist(meta$LOG) #looks good!! 
"Consider adding in more values of m = 0 when lit supports zero migration"
#if <0, negative response. if = 1, no response. if >1, positive response

#switching the variables doesnt change it, but here's the code if that seems to make more sense
#meta$LOG2 = log(meta$calcm) - log(meta$m_rate+1)
#hist(meta$LOG2)

#random effects model
results=rma(yi=LOG,vi=vi, data=meta, weighted=FALSE) #either way will be unweighted since the variance isnt available
results
forest(results)
funnel(results)


#measure=SMD "standardized mean difference"
#measure=SMDH "standardized mean difference with heteroscedastic population variances in the two groups (Bonett, 2008, 2009)"
#mrasure=ROM "log transformed ratio of means"
#FAILED::  datum2=escalc(measure="ROM", n1i=N_pop1, m1i=log(m_rate+1), sd1i=sd, n2i=N_pop2, m2i=log(calcm), sd2i=sd, data=meta)

#random effects model -- same output as above
#A fixed-effects model (with or without moderators) is fitted when using method="FE". Random/mixed-effects models are fitted by setting method equal to one of the following: "DL", "HE", "SJ", "ML", "REML", "EB", "HS", or "GENQ". Default is "REML".
#results2=rma(yi=LOG,vi=vi, data=meta, method="REML") #same output as method "HE"
#same results as above analysis

#Moderator analysis 1: use number of loci as a moderator
results3=rma(yi=LOG, vi=vi, mods=~ Loci, data=meta,method="DL") #can also do Loci-1 to get rid of the intercept. logically I don't know that we'd need to
#check here for info on intercepts https://www.metafor-project.org/doku.php/tips:models_with_or_without_intercept
#"When the model only includes continuous (i.e., numeric) predictors/moderators, then removing the intercept does just that: it removes the intercept. Hence, the model above forces the intercept to be 0, that is, it assumes that at the equator, the (average) log risk ratio is exactly zero. This seems like a rather strong assumption to make, so I would not recommend doing so. In fact, only in rare cases is it ever appropriate to remove the intercept term from a regression model (that involves only continuous predictors/moderators)."
results3
forest(results3)

#Moderator analysis 2: use study as a moderator
#this is the same as using species as moderator, as each study contained the same species
results4=rma(yi=LOG, vi=vi, mods=~ N, data=meta,method="DL")
forest(results4)
results4

#Moderator analysis 3: use comparison as a moderator (similar to random effect by id- i think?)
#this gave each row its own number, so this is basically just measuring ID
as.factor(meta$n)
results5=rma(yi=LOG, vi=vi, mods=~ n, data=meta,method="DL")
forest(results5)
results5
#FAILED: these values look weird. 

#Moderator analysis 4: use habitat as a moderator (terrestrial (t) and aquatic (a))
results6=rma(yi=LOG, vi=vi, mods=~ habitat, data=meta,method="DL")
forest(results6)
results6
#FAILED: these values look weird. 

results7=rma(yi=LOG, vi=vi, mods=~ Year, data=meta, method="DL")
forest(results7)
results7
#FAILED: these values look weird. 

#can test relationship of variables to see if they are independent. technically they SHOULDNT be
cor(meta$calcm, meta$m_rate)

#perform a log linear regression
#code found here https://www.scribbr.com/statistics/linear-regression-in-r/
reglog <- lm(log(calcm) ~ (log(m_rate+1)), data = meta)
summary(reglog)
logm = log(meta$m_rate+1)
logf = log(meta$calcm)
logm #check for problem values
NaN %in% meta$calcm
NaN %in% meta$m_rate
shapiro.test(residuals(reglog)) #Y

#perform a linear regression
reg <- lm(calcm ~ m_rate, data = meta)
summary(reg)
shapiro.test(residuals(reg)) #NOPE
"Log regression fits the data better!"

#website for ggplot2 https://www.guru99.com/r-scatter-plot-ggplot2.html#:~:text=Scatter%20Plot%20in%20R%20using%20ggplot2%20(with%20Example),he%20can%20start%20to%20explore%20the%20dataset.
#website for different smooths in ggplot2 https://stats.idre.ucla.edu/r/faq/how-can-i-explore-different-smooths-in-ggplot2/#:~:text=p%20+%20stat_smooth%20(method%20=%20%22lm%22,%20formula%20=,function%20can%20easily%20fit%20polynomials%20of%20arbitrary%20degree
library(ggplot2)
"The plot for a regular linear relationship is a bad fit"
scat <- ggplot(meta, aes(x = calcm, y = m_rate)) +
  geom_point() +
  stat_smooth(method = "lm", #linear regression
              col = "#C42126", #red color of line
              se = FALSE, #dont display standard error
              size = 1) #size of line

scat +
  labs(
    x = "Fst",
    y = "Coalescent",
    color = "Type",
    title = "Relationship of Fst and Coalescent values before log",
    caption = "Fst and coalescent values compared across available literature"
  )

#The log relationship looks better but still not a good fit
scatlog <- ggplot(meta, aes(x = logf, y = logm, size = 3)) +
  geom_point() +
  stat_smooth(method = "lm", #linear regression
              col = "#C42126", #red color of line
              se = FALSE, #dont display standard error
              size = 2) #size of line

scatlog +
  labs(
    x = "log calculated Fst",
    y = "log Coalescent",
    color = "Type",
    title = "Relationship of LOG Fst and Coalescent values",
    caption = "Fst and coalescent values compared across available literature"
  )

#locally weighted regression to see the pattern in the points... there isnt one
scatlog_weight <- ggplot(meta, aes(x = logf, y = logm, size = 3)) +
  geom_point() +
  stat_smooth(method = "loess", #linear regression
              col = "#C42126", #red color of line
              se = TRUE, #yes display standard error
              size = 2) #size of line

scatlog_weight +
  labs(
    x = "log calculated Fst",
    y = "log Coalescent",
    color = "Type",
    title = "Relationship of LOG Fst and Coalescent values",
    caption = "Fst and coalescent values compared across available literature"
  )

ggsave("Fst_v_Coal.png")
#one idea could be to classify as teh measurement (tajima, Fu, etc) titled "type" and then color by that on the chart. 
#ggplot(meta, aes(x = FST, y = COAL)) +
#geom_point(ae(color = factor(type)))

#~~~~~~~~~~~~~~~~~~~~~~~`
#OTHER ANALYSIS
#USE OTHER DATASET WITH MORE DATA
data = read.csv("MammalMetaDataNEW.csv")

"what is the relationship between heterozygosity/allelic richness and Fst?"
het <- lm(Ob_Het1 ~ Fstmicro, data = data)
summary(het)

ar <- lm(AllelicRichness1 ~ Fstmicro, data = data)
summary(ar)

"What is the relationship between mtDNA FST and microsat FST?"
mark <- lm(FstmtDNA ~ Fstmicro, data = data)
summary(mark)

#other data that can maybe added later for a GLMM: ICUN status, family/class, aquatic/terrestrial, clade/pylogenetic relationship analysis

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
meta$LOG

#forest plot
#random effects model
df = data.frame(LOG=meta$LOG, n=meta$n)
df = df[complete.cases(df),]
plot(-100,-100, xlim=c(-7, 2), ylim=c(1, nrow(df)), xlab="log response ratio", ylab="comparison", main="Random Effects model of Fst vs. migration rate")
points(y=seq(1,nrow(df)), x=df$LOG, pch=19, col="black")
abline(v=0)
