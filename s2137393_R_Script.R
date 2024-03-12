######## 1. HYPOTHESIS TESTING
##H0: L0 = L12 [Levels of bio-markers remains unchanged from the Inclusion Phase to a time of 12months]
##H0: L0 != L12 [Levels of bio-markers remains changes from the Inclusion Phase to a time of 12months]

###Setting the Environment
#attaching the bio-markers data of phase- At Inclusion vs 12 months
attach(Bio_0m_vs_12m_refined)
#the variables or bio-markers whose value are 0 it is showing N/A in the table so 
# replacing with 0 in the N/A s
Bio_0m_vs_12m_refined[is.na(Bio_0m_vs_12m_refined)] <- 0
#attaching the covariates table
attach(Cov)
#the variables or Covariates whose value are 0 it is showing N/A in the table so 
# replacing with 0 in the N/A s
Cov[is.na(Cov)] <- 0
#Merging both the tables on the ground of PatientID
mergeTable<- merge(Cov,Bio_0m_vs_12m_refined)
#Viewing the merged tables
View(mergeTable)

###Hypothesis Testing: Two Sided Test using Welch 2 Sample Test
#The hypothesis test is getting done for all the 9 bio-markers at a confidence 
#level of 95%, where alpha is 5%.
#After doing test for all the 9 bio-markers the t-test is also conducted for 
# the VAS At inclusion vs 12 months to establish the relevancy of the test result of
# the bio-markers.
t.test(`IL  8`,`IL  8_12m`,conf.level = 0.95)
t.test(`VEGF  A`,`VEGF  A_12m`,conf.level = 0.95)
t.test(OPG, OPG_12m,conf.level = 0.95)
t.test(`TGF  beta  1`,`TGF  beta  1_12m`,conf.level = 0.95)
t.test(`IL  6`,`IL  6_12m`,conf.level = 0.95)
t.test(CXCL9, CXCL9_12m,conf.level = 0.95)
t.test(CXCL1, CXCL1_12m,conf.level = 0.95)
t.test(`IL  18`,`IL  18_12m`,conf.level = 0.95)
t.test(`CSF  1`, `CSF  1_12m`,conf.level = 0.95)
t.test(`VAS-at-inclusion`,`Vas-12months`,conf.level = 0.95)

# Since multiple testing is done so the probability of making Type-I error increases.
#It is then called Family Wise Error rate (FWER)
#The probability of making at least one type I error Fwer is calculated
fwer<- 1-(1-0.05)^9
fwer

###t-test with "bonferroni" correction for the t-tests.
#producing the p-values and storing an relevant variable
p_IL8<- t.test(`IL  8`,`IL  8_12m`,conf.level = 0.95)$p.value
p_VEGF_A<- t.test(`VEGF  A`,`VEGF  A_12m`,conf.level = 0.95)$p.value
p_OPG<-t.test(OPG, OPG_12m,conf.level = 0.95)$p.value
p_TGF<-t.test(`TGF  beta  1`,`TGF  beta  1_12m`,conf.level = 0.95)$p.value
p_IL6<-t.test(`IL  6`,`IL  6_12m`,conf.level = 0.95)$p.value
p_CXCL9<-t.test(CXCL9, CXCL9_12m,conf.level = 0.95)$p.value
p_CXCL1<-t.test(CXCL1, CXCL1_12m,conf.level = 0.95)$p.value
p_IL18<-t.test(`IL  18`,`IL  18_12m`,conf.level = 0.95)$p.value
p_CSF1<-t.test(`CSF  1`, `CSF  1_12m`,conf.level = 0.95)$p.value

#creating a table with all the 9 p-values
p <- data.frame(p_IL8,p_VEGF_A,p_OPG,p_TGF,p_IL6,p_CXCL9,p_CXCL1,p_IL18,p_CSF1)

#adjusting the p-values with bonferroni correction
p<-  p.adjust(p, method = "bonferroni", n= length(p))
#viewing the adjusted p's
View(p)
p


##### 2.REGRESSION MODELLING
### a. Diving the merged table in 2 portions of 80% patient data and 20%
first_portion <- floor(0.8 * nrow(mergeTable))
## setting a seed to make the partition of the table reproducible [making sure 
# every time the same data/patient's data belongs to respective 80% and 20%]
set.seed(123)
#assigning the 80% data index
first_portion_idx <- sample(seq_len(nrow(mergeTable)), size = first_portion)
#creating 2 tables with the 80% and remainging 20% of data
first80 <- mergeTable[first_portion_idx, ]
remaining20 <- mergeTable[-first_portion_idx, ]
#Multiple linear regression model for 80% of the Data taking the VAS-12 months as
# Y (dependent variable) and all the covariates and bio-markers as the independent
#variables
regression80<- lm(`Vas-12months`~(Age+`Sex (1=male, 2=female)`+`Smoker (1=yes, 2=no)`+`IL  8`+`VEGF  A`+OPG+ `TGF  beta  1`+
                      `IL  6`+CXCL9+CXCL1+`IL  18`+`CSF  1`),first80)
#Doing a summary of the model
summary(regression80)
#Presenting the Fitted Parameter Values in a Table.
fittedTable<- data.frame(fitted(regression80))
View(fittedTable)


### c. Prediction using the date of remaining 20% data with the formerly created 
# regression model(80% data Model)
# Creating the independent variable of the remaining 20% data to use for the prediction.
remaining20_newX <-remaining20[c("Age","Sex (1=male, 2=female)","Smoker (1=yes, 2=no)",
                                 "IL  8","VEGF  A","OPG", "TGF  beta  1", "IL  6", 
                                 "CXCL9","CXCL1","IL  18","CSF  1")]
#Prediction using the remaining 20% data using the regression model of the 80% data
m2<- predict(regression80,remaining20_newX , level=0.95)
#Plotting the predicted data vs the observed or actual data of the 20% remaining pateints.
plot(m2,remaining20$`Vas-12months`
     ,xlab = "Predicted VAS_12m"
     ,ylab = "Observed VAS_12m")
#Adding a straight line in the plot
abline(a = 0,                                        
       b = 1,
       col = "red",
       lwd = 3)
#Comparing the predicted and actual VAS-12months of the remaining 20% patients
#in a table.
compareData <- data.frame(actual= remaining20$`Vas-12months`, predicted =m2)
View(compareData)
