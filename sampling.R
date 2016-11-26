#######################################
###    draw samples and do stats    ###
#######################################

#######################################
###     data prep and cleaning      ###
#######################################

# read in data, if not already in memory
load("station_crimes.RData")
data <- data.final

# split out into two data sets
crimes_Full <- data[data$CATEGORY == "property",]
crimes_Full <- crimes_Full[,c(1:7, 9)]
names(crimes_Full)[8] <- "PROPERTY_CRIMES"
crimes_Full$VIOLENT_CRIMES <- data$NUM_CRIMES[data$CATEGORY == "violent"]

#######################################
###         stratification          ###
#######################################

# pilot data is first six months of 2014
pilotData <- crimes_Full[as.Date(crimes_Full$WKND_START) < "2014-07-01",]
surveyData <- crimes_Full[as.Date(crimes_Full$WKND_START) >= "2014-07-01",]

# generate summary stratum data
stratumData <- 
  data.frame(Nh = table(surveyData$WARD),
             sh = tapply(pilotData$PROPERTY_CRIMES, pilotData$WARD, sd))
colnames(stratumData) <- c("ID", "Nh", "sh")

#######################################
###        Neyman allocation        ###
#######################################

# define desired probability and error 
alpha <- 0.05
e <- 0.05
z.alpha.2 <- qnorm(alpha / 2, lower.tail = F)

# gives a starting point for Exact Optimal allocation
n.neyman <- neyman.alloc(N.h = stratumData$Nh, S.h = stratumData$sh,
                         eps = e, alpha = alpha)

#######################################
###     Exact Optimal allocation    ###
#######################################

exopt.n <- exact.opt.alloc(N.h = stratumData$Nh, S.h = stratumData$sh, 
                           n = n.neyman$n, eps = e, alpha = alpha)

stratumData <- cbind(stratumData, n.h = exopt.n$n.h)

#######################################
###   statistics under stratified   ###
#######################################

# a single statistic 
set.seed(2)
sampleData <- sapply(stratumData$ID, FUN = function(x) {
  numUnits <- stratumData$n.h[x]
  sample(surveyData$PROPERTY_CRIMES[surveyData$WARD == x])[1:numUnits]
})

means <- sapply(sampleData, mean)
ybar <- sum(means * stratumData$Nh) / sum(stratumData$Nh)

# calculate a confidence interval for the population mean
sampleVars <- sapply(sampleData, FUN = var)
N_h_Over_N_Squared <- (stratumData$Nh / sum(stratumData$Nh)) ^ 2

# finite population correction
fpc <- (stratumData$Nh - stratumData$n.h) / stratumData$Nh
var_y_bars <- fpc * sampleVars / stratumData$n.h
var_y_bar_str <- sum(N_h_Over_N_Squared * fpc * sampleVars / stratumData$n.h)
print(c(confInt <- c(lowerBound = ybar - z.alpha.2 * sqrt(var_y_bar_str), 
                     upperBound = ybar + z.alpha.2 * sqrt(var_y_bar_str))))


# do random draws
ybars_Stratified <- sapply(1:500, FUN = sample_mean, data = surveyData, 
                           data.str = stratumData, stratify = T)
hist(ybars_Stratified, 
     main = "Stratified Random Sampling: Distribution of Estimates", 
     xlab = "Estimates")
abline(v = mean(surveyData$PROPERTY_CRIMES), col = "red", lwd = 3)
abline(v = mean(surveyData$PROPERTY_CRIMES) + 0.05, col = "purple", 
       lwd = 3, lty = 2)
abline(v = mean(surveyData$PROPERTY_CRIMES) - 0.05, col = "purple", 
       lwd = 3, lty  = 2)

# paste("Percent of samples within .05: ", 
#       length(which(abs(ybars_Stratified - 
#                          mean(surveyData$PROPERTY_CRIMES)) < 0.05)) / 5,
#       "%", sep = "")

#######################################
###       statistics under SRS      ###
#######################################

# do random draws
ybars_SRS <- sapply(1:500, FUN = sample_mean, data = surveyData, 
                    data.str = stratumData, stratify = F)
hist(ybars_SRS)
abline(v = mean(surveyData$PROPERTY_CRIMES), col = "red", lwd = 3)

plot(density(ybars_Stratified), lwd = 3, col = "blue", 
     main = "Comparison of Methods")
lines(density(ybars_SRS), lwd = 3, col = "green")
legend("topleft", legend = c("Stratified","Simple"), 
       col = c("blue", "green"), pch = 16)

#######################################
###  sample total under stratified  ###
#######################################

# a single statistic 
t.hat <- sum(means * stratumData$Nh)

# calculate a confidence interval, noting that the standard deviation of 
# \hat{t}_{str} is N * the standard deviation of \bar{y}_{str}
z <- z.alpha.2 * sum(stratumData$Nh) * sqrt(var_y_bar_str)
print(c(confInt <- c(lowerBound = t.hat - z, upperBound = t.hat + z)))
