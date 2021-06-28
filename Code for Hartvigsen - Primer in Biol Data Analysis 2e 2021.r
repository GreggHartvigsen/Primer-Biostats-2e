# Below is R code from the book:
# 
# A Primer in Biological Data Analysis and Visualization Using R
# Second Edition
# ISBN = 9780231202138 (paper)
# ISBN = 9780231202121 (hardcover)
# Gregg Hartvigsen
# Columbia University Press
# 2021
# http://cup.columbia.edu/book/a-primer-in-biological-data-analysis-and-visualization-using-r/9780231202138

#####################################################################
# CHAPTER 1
############

5-1 # Subtraction
2*3 # Multiplication
7/3 # Division
sqrt(9) # Use the sqrt() function to get the square root of 9
9^2 # 9 squared
log(3) # Natural logarithm of 3
log(3,10) # Log of 3 (base 10), or use log10(3)

1:5 # create an array of integers from 1 to 5

c(1, 2.5, 3, 4, 3.5) # combining five numbers into an array

dat = c(1, 2.5, 3, 4, 3.5) # store numbers the variable "dat  "

dat <- c(1, 2.5, 3, 4, 3.5) # store numbers the variable "dat  "

OoS = c("When", "on", "board", "H.M.S.", "Beagle....")

dog.mass = 25.2 # a 25.2 kg dog

sum(dat) # sum up all values in array "dat"
length(dat) # tells you how many numbers are in "dat"
sum(dat)/length(dat) # this calculates the mean
summary(dat) # more descriptive statistics for "dat"
dat/5 # divide each value in "dat" by 5
dat[5] # returns the fifth element in array "dat"
dat[1:3] # returns the first three elements of array "dat"

seq(1,10, by = 0.5) # sequence from 1 to 10 by 0.5

my.seq = seq(1,10, by = 0.5) # store result in variable "my.seq"

my.seq[12]

rep(c("A","B","C"), times = 2) # entire array twice
rep(c("A","B","C"), each = 2) # each element twice

p = c(1/2,1/4,1/4) # three proportions saved in an array

curve(2*x^2 + 4*x - 7,-10,10)

hist(rnorm(10000))

sum(1,2,3,4)

## install.packages("UsingR")

library(UsingR)

#####################################################################
# CHAPTER 2
############

male.hts = c(1.82, 1.79, 1.70, 1.76) # heights in meters
female.hts = c(1.55, 1.62, 1.58, 1.68)

Mean(male.hts) # calculates the arithmetic mean
median(female.hts) # calculates the median for continuous data

boxplot(male.hts,female.hts, names = c("Males","Females"), 
        xlab = "Sex",ylab = "Heights (m)", cex.lab = 1.5)

W = "ftp://aftp.cmdl.noaa.gov/products/trends/co2/co2_weekly_mlo.txt"
CO2 = read.table(W, skip = 49)[,4:5]
plot(CO2, ylim = c(300,420), pch = 16, cex = 0.5,
     xlab = "Year", ylab = "Atmospheric CO2 (ppm)")

#####################################################################
# CHAPTER 3
############

signif(mean(3.0 + 5.275),2)

a = 3.141592654 # pi
b = 3141592654 # a x 10^9
round(a,2) # control number of decimal places displayed
signif(a,2) # this controls significant digits
signif(b,2)


males = c(1.72, 1.59, 1.70)
females = c(1.55, 1.62, 1.58)

height.dat = data.frame(males,females)
height.dat

height.dat$males # this returns an array of data for males

summary(height.dat)

height.dat.stacked = stack(height.dat)
height.dat.stacked

names(height.dat.stacked) = c("height","sex") 

names(height.dat.stacked)

unstack(height.dat.stacked)

which(height.dat.stacked$sex == "females")

female.rows = which(height.dat.stacked$sex == "females")
height.dat.stacked$height[female.rows]

males = subset(height.dat.stacked, sex == "males")

males$height

height.dat.stacked$height # unsorted
sort(height.dat.stacked$height) # sorted low to high

order(height.dat.stacked$height)

height.dat.stacked[order(height.dat.stacked$height),]

#####################################################################
# CHAPTER 4
############

x = c(5,3,6,7,4)
mean(x)

vals = c(1,2,3)
weights = c(7,3,5) # 7 ones, 3 twos, and 5 threes
weighted.mean(vals,weights)

median(c(2,3,4,5,8)) # 4 is the middle value
median(c(2,3,4,8)) # the middle lies between 3 and 4

set.seed(100) # do this so you will have the same data
dat.raw = sample(1:10, 500, replace = T) # get 500 random 
# values from 1 to 10
dat.table = table(dat.raw) # cross tabulation for the data
dat.table # here are the data, ordered up by their frequency

as.numeric(names(dat.table)[which.max(dat.table)])

set.seed(100) # do this so you will have the same data
a = hist(rnorm(1000)) # store data from hist() in "a"
the.mode = a$mids[which.max(a$counts)]
cat("The mode is the bin centered on", the.mode,"\n")

set.seed(100) # do this so you will have the same data
x = rnorm(1000) # 1000 nums values from the "standard 
# normal distribution"

range(x) # returns the smallest and largest values

diff(range(x)) # for this, diff() gives us what we want

sd(x) # sd for the standard deviation

sd(x)^2
var(x)

SEM = sd(x)/sqrt(length(x)) 
SEM

n = length(x) # the size of a sample
CI95 = SEM * qt(0.975, df = n - 1)
CI95

IQR(x) # returns the interquartile range of the x array

library(UsingR)

simple.eda(rnorm(1000))

set.seed(10) # do this so you will have the same data
x1 = rnorm(1000) # 1000 values from standard normal dist.
x2 = runif(1000) # 1000 values from uniform dist.
x3 = rgamma(1000,2,1) # 1000 values from a gamma distribution
x4 = c(rep(3,10),1,2,4,5) # leptokurtic distribution

par(mfrow = c(2,2)) # create 2 x 2 graphics panel
hist(x1, main = "Normal Distribution", las = 1)
hist(x2, main = "Uniform Distribution", las = 1)
hist(x3, main = "Gamma Distribution", las = 1)
hist(x4, main = "Leptokurtic Distribution",breaks = 0:5, las = 1)

shapiro.test(x1)$p.value # large p -> fail to reject Ho
shapiro.test(x2)$p.value # small p -> reject Ho
shapiro.test(x3)$p.value # small p -> reject Ho
shapiro.test(x4)$p.value # small p -> reject Ho

library(e1071)

skewness(x1)
skewness(x2)
skewness(x3)
skewness(x4)

kurtosis(x1)
kurtosis(x2)
kurtosis(x3)
kurtosis(x4)

x = runif(100)+1 # creates 100 random values (1 <= x <= 2)

x[36] = 153 # makes the 36th value erroneously 153

which(x > 2) # find index of number(s) that are > 2

x[36]

x[36] = 1.53

x = x[-36] # -36 means removes the 36th value, then copies 
# the remaining values back into x

set.seed(10) # do this so you will have the same data
Data = exp(rnorm(100, mean = 5, sd = 0.5))
log.Data = log(Data) # the natural logarithm of the x1 data
par(mfrow = c(1,2))
hist(Data, xlab = "Data", main = "",ylim = c(0,30), las = 1,
     cex.lab = 1.5)
hist(log.Data, xlab = "ln(Data)", main = "", cex.lab = 1.5,
     xlim = c(3.5,6.5),ylim = c(0,20), las = 1)

shapiro.test(Data) # these are NOT normally distributed
shapiro.test(log.Data) # transformed -> normally distributed

#####################################################################
# CHAPTER 5
############

set.seed(100) # do this so you will have the same data
x = rnorm(100, mean = 5) # 100 random nums, mean = 5, sd = 1
y = rnorm(100, mean = 6) # 100 random nums, mean = 6, sd = 1 
z = rnorm(100, mean = 8) # 100 random nums, mean = 8, sd = 1

par(mfrow = c(1,2)) # make graph window 1 row and 2 columns
plot(x[1:20]) # left plot with no graphics parameters
plot(x[1:20], xlim = c(0,20),ylim = c(0,10), las = 1,
     xlab = "My X-Axis Label", ylab = "My Y-Axis Label", 
     main = "My Custom Title", pch = 16, cex = 0.75,
     cex.lab = 1.5) # a professional-looking graph!

par(mfrow = c(1,2))
hist(x, cex.lab = 1.5, main = "", las = 1)
abline(v = mean(x),lwd = 5)
hist(x, cex.lab = 1.5, main = "", breaks = 34, 
     las = 1, ylim = c(0,10))
abline(v = mean(x),lwd = 5)
par(mfrow = c(1,1))

boxplot(x,y,z,names = c("Control","Low","High"), 
        xlab = "Fertilizer Treatment", las = 1,
        ylab = "Tree Height (m)", 
        cex.lab = 1.5, ylim = c(0, max(c(x,y,z))))

Ht = c(mean(x),mean(y),mean(z))
barplot(Ht,
        xlab = "Fertilizer Treatment",
        ylab = "Tree Height (m)", las = 1,
        names = c("Low","Medium","High"), cex.lab = 1.5)
abline(h=0)

library(plotrix)
set.seed(10)
a = rnorm(1000,mean = 100,sd = 14)
# b = rnorm(1000, mean = 100, sd = 10)
b = rgamma(1000, shape = .1, scale = 5)*10 + 95
par(mfrow = c(1,2))
pos = barplot(c(mean(a),mean(b)),ylim = c(0,200), 
              xlab = "Treatments", ylab = "Biomass (g)", 
              names = c("A","B"),cex.lab = 1.5, las = 1)
abline(h=0)
M = c(mean(a),mean(b))
S = c(sd(a),sd(b))
plotCI(pos,M,S,pch = NA, add = T)
boxplot(a,b, ylab = "Biomass (g)",names = c("A","B"),
        ylim = c(0, 200), cex.lab = 1.5, xlab = "Treatments", 
        pch = 16, cex = .5,las = 1)
par(mfrow = c(1,1))

females = c(36.3,29.5)
males = c(36.2,28.9) # first is biology, then non-bio majors
MCAT = matrix(c(females,males),byrow = T, nrow = 2)
MCAT

par(mfrow = c(1,2)) # graphics panel with 1 row, 2 columns
MCAT = matrix(c(females,males),byrow = T, nrow = 2)
barplot(MCAT,beside = T, names = c("Biology","Nonbiology"),
        ylim = c(0,60), xlab = "Major", cex.lab = 1.5, las = 1,
        ylab = "Mean MCAT Score", legend.text = c("Males","Females"))
abline(h=0)
# Alternatively, separate by major
bio = c(36.2,36.3) # first males, then females
non.bio = c(28.9,29.5)
MCAT = matrix(c(bio,non.bio),byrow = T, nrow = 2)
barplot(MCAT, beside = T, names = c("Males","Females"),
        ylim = c(0,60), xlab = "Sex", las = 1,
        cex.lab = 1.5, ylab = "Mean MCAT Score",
        legend.text = c("Biology","Nonbiology"))
abline(h=0)

par(mar = c(5.1,5.1,2.1,1.1))
M.artiodactyl = c(37800, 196500,69100, 325000,21500,58600,20500)
BMR.artiodactyl = c(9318,41242,19120,51419,8308,25609,5945)
plot(M.artiodactyl, BMR.artiodactyl,
     xlim = c(1000,300000), ylim = c(400,100000),
     log = "xy", cex.lab = 1.5, xaxt = "n", yaxt = "n",
     xlab = "Mass (g)",
     ylab = expression(paste("BMR (ml O"[2]," hr"^-1,")")))
axis(2,at = c(1e3,1e4,1e5), 
     labels = expression("10"^3,"10"^4,"10"^5),las = 1)
axis(1, at = c(1e3,1e4,1e5), 
     labels = expression("10"^3,"10"^4,"10"^5))
M.canids = c(3600,10000,7720,5444,1215,1769,4440)
BMR.canids = c(1374,2687, 3860,1524,583,887,2442)
points(M.canids, BMR.canids,#log = "xy",
       pch = 16, las = 1)
legend("topleft", legend = c("Artiodactyla","Canidae"), pch = c(1,16))

N = 20 # number of individuals in sample
time1 = rnorm(N, mean = 3)
time2 = rnorm(N, mean = 7)
plot(0, xlim = c(0.75,2.25),ylim = c(0,10),
     type = "n",xaxt = "n", xlab = "Time", 
     ylab = "Measurement", cex.lab = 1.5, las = 1)
axis(1,at = c(1,2),labels = c("Before","After"), 
     cex.axis = 1.25)
for (i in 1:length(time1)) { # add the lines one at a time
  lines(c(1,2),c(time1[i],time2[i]),lty = 1) 
}
points(rep(1,length(time1)),time1, pch = 16, cex = 1.5)
points(rep(2,length(time2)),time2, pch = 17, cex = 1.5)

x1 = c(10,20,30,40)
x2 = c(5,25,35,35)
my.dat = matrix(c(x1,x2),nrow = 2,byrow = T)
my.col = gray(seq(0,1.0,length=4))
layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE)) # set layout
pie(my.dat[1,], col = my.col,radius = 1,
    main = "Pie Chart 1")
pie(my.dat[2,], col = my.col,radius = 1,
    main = "Pie Chart 2")
leg.txt = c("Pie 1 data","Pie 2 data")
barplot(my.dat,beside = T, ylim = c(0,50), names = 1:4, las = 1,
        col = gray.colors(2), ylab = "Height (cm)",
        xlab = "Groups", cex.lab = 1.5, main = "Barplot")
legend("topleft",leg.txt, fill = gray.colors(2))
abline(h=0)

# min() gets smallest value; floor() rounds down
xmin = floor(min(c(x,y,z))) 
# max gets highest value; ceiling() rounds up
xmax = ceiling(max(c(x,y,z))) 
par(mfrow = c(3,1)) # set graphic window: 3 rows, 1 col
hist(x,xlim = c(xmin,xmax), las = 1)
hist(y,xlim = c(xmin,xmax), las = 1)
hist(z,xlim = c(xmin,xmax), las = 1)
par(mfrow = c(1,1)) # reset the graphics window

dat = data.frame(x,y,z)
pairs(dat)

#####################################################################
# CHAPTER 6
############

pot = 1:6

trmt = rep(c("control","fertilizer"), each = 3)

set.seed(11) # so you get the same order as I do (don't do if 
#              doing this for an experiment!)
pot = sample(pot) # randomize the array pot numbers
design = data.frame(trmt,pot)
design

power.t.test(delta = 30, sd = 25, power = 0.9)

sample(1:10,5) # from the array 1-10, randomly choose 5

#####################################################################
# CHAPTER 7
#####################

SAT.scores = c(1130, 1090, 1190, 1110, 1160, 1120, 1160, 1110,
               1080, 1160, 1050, 1120, 1030, 1010, 1080, 1090,
               1170, 1110, 1090, 1140)
your.SAT = 1120

par(mfrow = c(1,2))
hist(SAT.scores, main = "", las = 1, xlab = "SAT Scores")
abline(v = your.SAT, lwd = 3)
boxplot(SAT.scores, las = 1, ylab = "SAT Scores")
abline(h = your.SAT, lwd = 3)

shapiro.test(SAT.scores)$p.value

par(mfrow = c(1,2))
library(plotrix) # install this package, if necessary
boxplot(SAT.scores, ylab = "SAT Scores", cex.lab = 1.5, 
        ylim = c(1000,1200), las = 1)
M = mean(SAT.scores)
s = sd(SAT.scores)
SEM = s/sqrt(length(SAT.scores))
CI95 = qt(0.975, df = length(SAT.scores) - 1)*SEM
a = barplot(M, ylim = c(1000,1200), xpd = F, cex.lab = 1.5,
            ylab = "Mean SAT Scores", las = 1)
abline(h=1000)
plotCI(a, M, CI95, pch = NA, add = T, lwd = 2)

t.test(SAT.scores, mu = your.SAT) # mu value tested against sample

SEM = sd(SAT.scores)/sqrt(20)
(mean(SAT.scores) - your.SAT)/SEM

breath = c(33, 79, 41, 41, 25, 51, 50, 46, 53, 61)

par(mfrow = c(1,2))
hist(breath, xlab = "Breath (sec)", las = 1, main = "")
abline(v=65, lwd = 3)
boxplot(breath, ylab = "Breath (sec)", las = 1)
abline(h=65, lwd = 3)

shapiro.test(breath)

t.test(breath, mu = 65, alt = "l")

jumpers = c(48, 118, 136, 134, 129, 123, 119, 87, 106, 119)
par(mfrow = c(1,2))
hist(jumpers)
abline(v = 137, lwd = 3)
boxplot(jumpers)
abline(h=137, lwd = 3)

shapiro.test(jumpers)$p.value

the.frog = 137 # distance in centimeters
wilcox.test(jumpers, mu = the.frog, alt = "g")

first = c(106.70, 72.00, 58.80, 48.00, 53.53, 35.93, 39.91,
          31.00, 45.85, 78.50)
second = c(129.00, 101.00, 64.20, 58.00, 64.78, 50.92, 48.50, 
           42.09, 70.00, 124.60)

breath.diff = second - first # create a new variable

par(mfrow = c(1,2))
boxplot(first, second, names = c("First","Second"),
        xlab = "Attempt", ylab = "Time (sec)", 
        cex.lab = 1.5, las = 1)
boxplot(breath.diff, ylim = c(-10,50),cex.lab = 1.5,
        ylab = "Difference (sec)", las = 1)
abline(h = 0, lty = 3, lwd = 3)

shapiro.test(breath.diff) # test this single sample for normality

t.test(breath.diff, alt = "g")

t.test(second, first, paired = T, alt = "g") # T for TRUE

before = c(75, 53, 78, 89, 96, 86, 69, 87, 73, 87)
after = c(124, 103, 126, 193, 181, 122, 120, 197, 127, 146)

LDL.diff = after-before
shapiro.test(LDL.diff)$p.value 

wilcox.test(LDL.diff)

cont = c(64.7, 86.6, 67.1, 62.5, 75.1, 83.8, 
         71.7, 83.4, 90.3, 82.7) # control plants
fert = c(110.3, 130.4, 114.0, 135.7, 129.9,  
         98.2, 109.4, 131.4, 127.9, 125.7) # fertilized treatment plants

boxplot(cont, fert, names = c("Control","Fertilizer"), 
        xlab = "Treatment", ylab = "Plant Height (cm)", 
        ylim = c(60,140), cex.lab = 1.5, las = 1)

shapiro.test(cont)$p.value # normal?
shapiro.test(fert)$p.value # normal?

var.test(cont,fert)

t.test(cont, fert, alt = "l",var.equal = TRUE)

t.test(fert, cont, alt = "g", var.equal = TRUE)

A = c(5800, 5000, 6500, 5000, 5200, 7900, 5300, 5900, 7600, 5400)
B = c(7100, 8600, 8900, 7500, 7300, 7100, 7400, 8300, 7000, 7000)

boxplot(A,B, names = c("Treatment A","Treatment B"), las = 1,
        ylab = "Milk Production (kg)", cex.lab = 1.5)

shapiro.test(A)$p.value # normal?
shapiro.test(B)$p.value # normal?

wilcox.test(A,B)

#####################################################################
# CHAPTER 8
#####################

Low = c(52.3, 48.0, 42.3, 50.8, 53.3, 45.1)
Med = c(50.4, 53.8, 53.4, 58.1, 56.7, 61.2)
High = c(66.3, 59.9, 57.1, 61.3, 58.3, 55.4)
fish.dat = data.frame(Low,Med,High)
fish.stacked = stack(fish.dat)

names(fish.stacked) = c("Mass","Trmt")

boxplot(fish.stacked$Mass ~ fish.stacked$Trmt)

tapply(fish.stacked$Mass, fish.stacked$Trmt, shapiro.test)

bartlett.test(fish.stacked$Mass, fish.stacked$Trmt)

class(fish.stacked$Trmt)

fish.stacked$Trmt = as.factor(fish.stacked$Trmt)

fish.aov = aov(fish.stacked$Mass ~ fish.stacked$Trmt)

summary(fish.aov)

summary(lm(fish.stacked$Mass ~ fish.stacked$Trmt))

TukeyHSD(fish.aov)

contrast.labels = c("a","b","b")

par(mfrow = c(1,2))
library(plotrix)
boxplot(fish.stacked$Mass ~ fish.stacked$Trmt, # that's a tilde symbol
        names = c("L","M","H"), ylim = c(0,80),
        xlab = "Feeding Rate", ylab = "Fish Mass (g)", 
        las = 1, cex.lab = 1.5)
data.tops = tapply(fish.stacked$Mass, fish.stacked$Trmt, max)
text(1:3, data.tops, labels = contrast.labels, pos = 3)
#------------------------------------------------------
M = tapply(fish.stacked$Mass,fish.stacked$Trmt,mean)
S = tapply(fish.stacked$Mass,fish.stacked$Trmt,sd)
SEM = S/sqrt(6) # because S is an array, SEM is an array
CI95 = qt(1-0.05/2,5) * SEM
a = barplot(M, ylim = c(0,80), las = 1, cex.lab = 1.5,
            names.arg = c("L", "M", "H"), 
            xlab = "Feeding Rate", ylab = "Fish Mass (g)") 
abline(h = 0)
plotCI(a,M,CI95, add = T, pch = NA)
text(a, M+CI95, labels = contrast.labels, pos = 3)

sp.A = c(7.03, 7.22, 7.16, 7.15, 7.33, 7.22, 7.77, 
         7.75, 7.33, 7.17)
sp.B = c(5.84, 5.59, 5.38, 5.38, 5.35, 5.87, 5.29, 
         5.75, 5.85, 5.33)
sp.C = c(7.09, 7.49, 7.06, 7.19, 7.10, 7.22, 7.45, 
         6.20, 6.96, 6.02)

tail.dat = data.frame(sp.A, sp.B, sp.C)
tail.dat = stack(tail.dat)
names(tail.dat) = c("length","species")

boxplot(tail.dat$length ~ tail.dat$species,names = c("A","B","C"), 
        ylab = "Tail Length (cm)", las = 1,
        xlab = "Species",cex.lab = 1.5)

tapply(tail.dat$length, tail.dat$species, shapiro.test)

tapply(log(tail.dat$length), tail.dat$species, shapiro.test)

kruskal.test(tail.dat$length ~ tail.dat$species)

LW.LN = c(3.84,4.47,4.45,4.17,5.41,3.82,3.83,4.67)
LW.HN = c(8.84,6.54,7.60,7.22,7.63,9.24,7.89,8.17)
HW.LN = c(7.57,8.67,9.13,10.02,8.74,8.70,10.62,8.23)
HW.HN = c(16.42,14.45,15.48,15.72,17.01,15.53,16.30,15.58)

water = rep(c("LW","HW"), each = 16)
nutr = rep(c("LN","HN"), each = 8, times = 2)
plant.mass = c(LW.LN,LW.HN,HW.LN,HW.HN) # combine the data

water
nutr

plant.dat = data.frame(water,nutr,plant.mass)

head(plant.dat) # view the first 6 rows

boxplot(plant.dat$plant.mass ~ plant.dat$water * plant.dat$nutr, 
        ylab = "Biomass (g)", xlab = "Treatment", 
        las = 1, cex.lab = 1.5)

plant.dat$plant.mass[plant.dat$water == "LW" & 
                       plant.dat$nutr == "LN"]

shapiro.test(plant.dat$plant.mass[water == "LW" & 
                                    nutr == "LN"])$p.value
shapiro.test(plant.dat$plant.mass[water == "LW" & 
                                    nutr == "HN"])$p.value
shapiro.test(plant.dat$plant.mass[water == "HW" & 
                                    nutr == "LN"])$p.value
shapiro.test(plant.dat$plant.mass[water == "HW" & 
                                    nutr == "HN"])$p.value

plant.aov = aov(plant.dat$plant.mass ~ 
                  plant.dat$water * plant.dat$nutr)

summary(plant.aov)

interaction.plot(plant.dat$water, trace.label = "nutr", las = 1,
                 trace.factor = plant.dat$nutr, plant.dat$plant.mass)

TukeyHSD(plant.aov)

contrast.labels = c("c","b","b","a")

M = tapply(plant.dat$plant.mass,
           list(plant.dat$water, plant.dat$nutr), mean)
SD = tapply(plant.dat$plant.mass,
            list(plant.dat$water, plant.dat$nutr), sd)
SEM = SD/sqrt(8) # there are 8 observations per sample
CI95 = qt(0.975,7)*SEM
a = barplot(M, beside = TRUE, ylim = c(0,20), 
            xlab = "Nutrient Level", ylab = "Mass (g)",
            legend = c("High Water","Low Water"),
            cex.lab = 1.5, las = 1)
abline(h = 0)
library(plotrix)
plotCI(a, M, CI95, pch = NA, add = T)
text(a, M + CI95, labels = contrast.labels, pos = 3)


#####################################################################
# CHAPTER 9
#####################

Low = c(52.3, 48.0, 42.3, 50.8, 53.3, 45.1)
Med = c(50.4, 53.8, 53.4, 58.1, 56.7, 61.2)
High = c(66.3, 59.9, 57.1, 61.3, 58.3, 55.4)
fish.dat = data.frame(Low,Med,High)
fish.stacked = stack(fish.dat)

names(fish.stacked) = c("Mass","Trmt")

boxplot(fish.stacked$Mass ~ fish.stacked$Trmt)

tapply(fish.stacked$Mass, fish.stacked$Trmt, shapiro.test)

bartlett.test(fish.stacked$Mass, fish.stacked$Trmt)

class(fish.stacked$Trmt)

fish.stacked$Trmt = as.factor(fish.stacked$Trmt)

fish.aov = aov(fish.stacked$Mass ~ fish.stacked$Trmt)

summary(fish.aov)

summary(lm(fish.stacked$Mass ~ fish.stacked$Trmt))

TukeyHSD(fish.aov)

contrast.labels = c("a","b","b")

par(mfrow = c(1,2))
library(plotrix)
boxplot(fish.stacked$Mass ~ fish.stacked$Trmt, # that's a tilde symbol
        names = c("L","M","H"), ylim = c(0,80),
        xlab = "Feeding Rate", ylab = "Fish Mass (g)", 
        las = 1, cex.lab = 1.5)
data.tops = tapply(fish.stacked$Mass, fish.stacked$Trmt, max)
text(1:3, data.tops, labels = contrast.labels, pos = 3)
#------------------------------------------------------
M = tapply(fish.stacked$Mass,fish.stacked$Trmt,mean)
S = tapply(fish.stacked$Mass,fish.stacked$Trmt,sd)
SEM = S/sqrt(6) # because S is an array, SEM is an array
CI95 = qt(1-0.05/2,5) * SEM
a = barplot(M, ylim = c(0,80), las = 1, cex.lab = 1.5,
            names.arg = c("L", "M", "H"), 
            xlab = "Feeding Rate", ylab = "Fish Mass (g)") 
abline(h = 0)
plotCI(a,M,CI95, add = T, pch = NA)
text(a, M+CI95, labels = contrast.labels, pos = 3)

sp.A = c(7.03, 7.22, 7.16, 7.15, 7.33, 7.22, 7.77, 
         7.75, 7.33, 7.17)
sp.B = c(5.84, 5.59, 5.38, 5.38, 5.35, 5.87, 5.29, 
         5.75, 5.85, 5.33)
sp.C = c(7.09, 7.49, 7.06, 7.19, 7.10, 7.22, 7.45, 
         6.20, 6.96, 6.02)

tail.dat = data.frame(sp.A, sp.B, sp.C)
tail.dat = stack(tail.dat)
names(tail.dat) = c("length","species")

boxplot(tail.dat$length ~ tail.dat$species,names = c("A","B","C"), 
        ylab = "Tail Length (cm)", las = 1,
        xlab = "Species",cex.lab = 1.5)

tapply(tail.dat$length, tail.dat$species, shapiro.test)

tapply(log(tail.dat$length), tail.dat$species, shapiro.test)

kruskal.test(tail.dat$length ~ tail.dat$species)

LW.LN = c(3.84,4.47,4.45,4.17,5.41,3.82,3.83,4.67)
LW.HN = c(8.84,6.54,7.60,7.22,7.63,9.24,7.89,8.17)
HW.LN = c(7.57,8.67,9.13,10.02,8.74,8.70,10.62,8.23)
HW.HN = c(16.42,14.45,15.48,15.72,17.01,15.53,16.30,15.58)

water = rep(c("LW","HW"), each = 16)
nutr = rep(c("LN","HN"), each = 8, times = 2)
plant.mass = c(LW.LN,LW.HN,HW.LN,HW.HN) # combine the data

water
nutr

plant.dat = data.frame(water,nutr,plant.mass)

head(plant.dat) # view the first 6 rows

boxplot(plant.dat$plant.mass ~ plant.dat$water * plant.dat$nutr, 
        ylab = "Biomass (g)", xlab = "Treatment", 
        las = 1, cex.lab = 1.5)

plant.dat$plant.mass[plant.dat$water == "LW" & 
                       plant.dat$nutr == "LN"]

shapiro.test(plant.dat$plant.mass[water == "LW" & 
                                    nutr == "LN"])$p.value
shapiro.test(plant.dat$plant.mass[water == "LW" & 
                                    nutr == "HN"])$p.value
shapiro.test(plant.dat$plant.mass[water == "HW" & 
                                    nutr == "LN"])$p.value
shapiro.test(plant.dat$plant.mass[water == "HW" & 
                                    nutr == "HN"])$p.value

plant.aov = aov(plant.dat$plant.mass ~ 
                  plant.dat$water * plant.dat$nutr)

summary(plant.aov)

interaction.plot(plant.dat$water, trace.label = "nutr", las = 1,
                 trace.factor = plant.dat$nutr, plant.dat$plant.mass)

TukeyHSD(plant.aov)

contrast.labels = c("c","b","b","a")

M = tapply(plant.dat$plant.mass,
           list(plant.dat$water, plant.dat$nutr), mean)
SD = tapply(plant.dat$plant.mass,
            list(plant.dat$water, plant.dat$nutr), sd)
SEM = SD/sqrt(8) # there are 8 observations per sample
CI95 = qt(0.975,7)*SEM
a = barplot(M, beside = TRUE, ylim = c(0,20), 
            xlab = "Nutrient Level", ylab = "Mass (g)",
            legend = c("High Water","Low Water"),
            cex.lab = 1.5, las = 1)
abline(h = 0)
library(plotrix)
plotCI(a, M, CI95, pch = NA, add = T)
text(a, M + CI95, labels = contrast.labels, pos = 3)

#####################################################################
# CHAPTER 10
#####################

obs = c(50,20) # observed counts
expP = c(0.75,0.25) # expected probabilities

chisq.test(obs, p = expP)

16*18/36

mat = matrix(c(10,6,8,12), byrow = TRUE, nrow = 2)
mat

chisq.test(mat)

M = matrix(c(22,4,15,10),byrow = TRUE, nrow = 2)
M

barplot(M, beside = TRUE, ylim = c(0,30),
        xlab = "Treatment", ylab = "Number of Patients", 
        names = c("Resolved","Not Resolved"), 
        legend = c("DT","Cryo"), cex.lab = 1.5)
abline(h=0)

chisq.test(M)

chisq.test(M, correct = F) # not using Yates' correction

M = matrix(c(22,4,15,10),byrow = TRUE, nrow = 2)
fisher.test(M)

#####################################################################
# CHAPTER 11
#####################

dat = c(6,3,4,5,3,1)

my.mean = function (x) {
  ans = sum(x)/length(x)
  return (ans)
}

mean(dat)
my.mean(dat)

mean.default
lm.line = function (x, y) {
  lines(x,fitted(lm(y~x)))
}
set.seed(100) # do this so you will have the same data
x = 5:20
y = 0.4*x + 2 + 7.5*runif(length(x))
mod = lm(y~x)
par(mfrow = c(1,2)) # create a two-panel graphics window
plot(x,y, xlim = c(0,25), ylim = c(0,20), pch=16, 
     cex.lab = 1.5, las = 1)
abline(mod)
plot(x,y, xlim = c(0,25), ylim = c(0,20), pch=16, 
     cex.lab = 1.5, las = 1)
lm.line(x,y)

lm.line = function (x, y) {
  lines(x,fitted(lm(y~x)))
}

set.seed(100) # do this so you will have the same data
x = 5:20
y = 0.4*x + 2 + 7.5*runif(length(x))
mod = lm(y~x)

par(mfrow = c(1,2)) # create a two-panel graphics window
plot(x,y, xlim = c(0,25), ylim = c(0,20), pch=16, 
     cex.lab = 1.5, las = 1)
abline(mod)
plot(x,y, xlim = c(0,25), ylim = c(0,20), pch=16, 
     cex.lab = 1.5, las = 1)
lm.line(x,y)

plot(BOD$Time, BOD$demand, pch = 16, ylim = c(0,30),
     ylab = "BOD (mg/l)", xlab = "Time (Day)", las = 1, cex.lab = 1.5)
mod = lm(BOD$demand ~ BOD$Time)
lm.line(BOD$Time, BOD$demand) # line function from earlier
newx = BOD$Time # create a new variable for Time
prd = predict(mod, interval = c("confidence"), level = 0.95,
              type="response")
lines(newx, prd[,2], lty=2, lwd = 2)
lines(newx, prd[,3], lty=2, lwd = 2)

set.seed(101) # do this so you will have the same data
x = 1:10
y = round(5*x + 25 + 14*runif(length(x)),1)
y

mod1 = lm(y~x)
mod2 = lm(y ~ poly(x,9))
par(mfrow = c(1,2))
plot(x,y, xlim = c(0,11), ylim = c(0,100), 
     ylab = "Dolphin Mass (kg)", xlab = "Time (weeks)", 
     pch = 16, cex.lab = 1.5, las = 1, xaxt = "n")
axis(1, at = seq(0,11,by = 5))
lines(x,fitted(mod1),lwd = 2)
plot(x,y, xlim = c(0,11), ylim = c(0,100), 
     ylab = "Dolphin Mass (kg)", xlab = "Time (weeks)", 
     pch = 16, cex.lab = 1.5, las = 1, xaxt = "n")
axis(1, at = seq(0,11,by = 5))
xvals = seq(1,10,by=0.1)
yvals = predict(mod2,list(x = xvals))
lines(xvals,yvals,lwd = 2)

summary(mod1)
summary(lm(y~x))

plant.mass = c(4.0, 8.5, 12.2, 13.4, 15.0, 17.8, 19.3, 19.4, 21.2, 
               21.7, 23.4, 23.8, 24.1, 24.7, 24.9, 25.5) 
time = 1:16 # weeks
mod = lm(plant.mass ~ time)
plot(time, plant.mass, xlim = c(0,16), ylim = c(0,30), 
     ylab = "Plant Mass (g)", xlab = "Time (weeks)", 
     pch = 16, cex.lab = 1.5, las = 1)
lines(time,fitted(mod), lwd = 2)

plot(time, plant.mass, xlim = c(0,16), ylim = c(0,30), 
     ylab = "Plant Mass (g)", xlab = "Time (weeks)", 
     pch = 16, cex.lab = 1.5, las = 1)
curve(-0.12*x^2+3.1*x+3,from = 1,to = 16, add = TRUE)

plot(time, plant.mass, xlim = c(0,16), ylim = c(0,30), 
     ylab = "Plant Mass (g)", xlab = "Time (weeks)", 
     pch = 16, cex.lab = 1.5, las = 1)
curve(25 * (1 - exp(-(1*x))),from=0,to=16,add = TRUE)

mod = nls(plant.mass ~ a*(1 - exp(-b * time)), 
          start = list(a = 25, b = 1))

summary(mod)

plot(time, plant.mass, xlim = c(0,16), ylim = c(0,30), pch = 16,
     ylab = "Plant Mass (g)", xlab = "Time (weeks)", 
     cex.lab = 1.5, las = 1)
xv = seq(0,16,0.1)
yv = predict(mod, list(time = xv))
lines(xv, yv)
asymptote = sum(coef(mod)[1]) # asymptote is the 
# sum of these two coefficients
abline(h = asymptote, lwd = 2, lty = 2)

S = c(0.2,0.5,1,2,5) # x-axis data ([S])
v = c(0.576, 0.762, 0.846, 0.860, 0.911) # velocity

plot(S,v, xlab = "[S]", ylab = "Velocity", xlim = c(0,max(S)),
     ylim = c(0,1), cex.lab = 1.5, pch=16, las = 1)

mod = nls(v ~ (Vmax * S)/(Km + S),
          start = list(Vmax = 1, Km = .1))
summary(mod)

plot(S,v, xlab = "[S]", ylab = "Velocity", xlim = c(0,max(S)),
     ylim = c(0,1),cex.lab = 1.5, pch = 16, las = 1)
x = seq(0,max(S),by = 0.01) # an array seq. of x values
y = predict(mod,list(S = x)) # the predicted y values
lines(x,y, lwd = 2) # draw the line
asymptote = coef(mod)[1] # value of the asymptote
abline(h = asymptote)

Vmax = coef(mod)[1] # this is Vmax
Km = coef(mod)[2] # this is Km

my.exp = expression(Vmax * my.S/(Km+my.S))
my.deriv = D(my.exp,"my.S")
my.deriv # look a the derivative returned by D()

my.S = 1 # choose where on [S] to get the rate (tangent line)

plot(S,v, xlab = "[S]", ylab = "Velocity", xlim = c(0,max(S)),
     ylim = c(0,1), cex.lab = 1.5, pch = 16, las = 1)
der.slope = eval(my.deriv) # get slow using deriv. at my.S
der.y = eval(my.exp) # Ht of the best-fit function at my.S
der.int = der.y - der.slope*my.S # The intercept of tangent
lines(x,y, lwd = 2)
asymptote = Vmax
abline(der.int,der.slope,lwd = 2,lty = 2) # draw tangent
abline(v=my.S) # place a vert line at tangent 

cat("The velocity of the reaction at",
    my.S," is ",der.slope,"\n")

# install.packages("deSolve")
library(deSolve)  # library for solving diff. equations
num_yrs = 10
r = 0.2 # the growth rate parameter
N0 = 100 # starting population
xstart = c(N=N0) # create a "list" of starting values
parms = r # here's our one parameter
mod = function(t,x,parms) {
  N = x[1]
  with(as.list(parms) , {
    dN.dt = r*N
    res=c(dN.dt)
    return(list(res))
  })
}

time=seq(0,num_yrs, by = 0.1)  # set number of time steps
# RUN THE MODEL in the next line!
output = as.data.frame(ode(xstart,time,mod,parms)) 

plot(output$time,output$N, xlab = "Time", ylab = "N", 
     type = "l", cex.lab = 1.5, las = 1)

Num_Days = 20  # number of days to run simulation
B = 0.006      # transmission rate (0.006)
v = 0.3        # recovery rate (0.3)
So = 499	     # initial susceptible pop (299)
Io = 1	       # initial infectious pop (1)
Ro = 0	       # initial recovered pop (0)
xstart = c(S=So,I = Io, R = Ro)
parms = c(B, v)
times=seq(0,Num_Days,length=200)  # set up time steps

mod = function(t,x,parms) {
  S = x[1] # init num of susceptibles
  I = x[2] # init num of infectious
  R = x[3] # init num of recovered
  with(as.list(parms) , {
    dS = -B*S*I       # dS/dt
    dI = B*S*I - v*I  # dI/dt
    dR = v*I          # dR/dt
    res=c(dS,dI,dR)
    list(res)
  })}

output = as.data.frame(lsoda(xstart,times,mod,parms))  

plot(output$time, output$S, type="n",ylab="Abundance",
     xlab="Time", cex.lab = 1.5,
     ylim = c(0,So*1.6), las = 1)
lines(output$time,output$S,lty = 1,lwd=3)
lines(output$time,output$I,lty = 2,lwd=3)
lines(output$time,output$R,lty = 3,lwd=3)
leg.txt = c("S","I","R")
legend("topright",leg.txt,lwd=2,lty = 1:3, horiz = T)

#####################################################################
# CHAPTER 12
#####################

set.seed(20) # do this so you will have the same data

num.values = 100 # the number of observations per sample
num.samples = 5000 # the number of samples to draw

means.norm = numeric(num.samples) # will hold 5000 values
means.unif = numeric(num.samples)

for (i in 1:num.samples) {
  means.norm[i] = mean(rnorm(num.values))
  means.unif[i] = mean(runif(num.values))
}

par(mfrow = c(2,2)) 
hist(rnorm(num.samples),xlab = "X",main = "Normal Sample", 
     ylim = c(0,1000), las = 1)
hist(runif(num.samples),xlab = "X",main = "Uniform Sample", 
     ylim = c(0,1000), las = 1)
hist(means.norm, xlab = "X",main = "Means from Normal Dist", 
     ylim = c(0,1000), las = 1)
hist(means.unif,xlab = "X",main = "Means from Uniform Dist", 
     ylim = c(0,1000), las = 1)
par(mfrow = c(1,1)) 

shapiro.test(means.norm)$p # from normal pop
shapiro.test(means.unif)$p # from unif pop

set.seed(10) # do this so you will have the same data
# Step 1: Define variables
n.time.steps = 100 # how long to run simulation
pop.size = 50 # how many individuals there are
n.reps = 20 # how many populations to simulate
# Step 2: Create an empty plot for simulation
plot(0,ylim = c(0,1),xlim = c(1,n.time.steps), type = "n", 
     ylab = "P(Allele 1)",
     xlab = "Time Step", cex.lab = 1.5, las = 1)
abline(h=0.5,lty = 2, lwd = 3)
# Step 3: Run the simulation n.reps times
for (i in 1:n.reps) {
  # Step 3a: Create population of 0 and 1 alleles
  pop = c(rep(0,pop.size/2),rep(1,pop.size/2))
  # Step 3b: Store proportion of 1s in variable
  P = sum(pop)/pop.size
  # Step 3c: Run the simulation for this replicate
  for (j in 2:n.time.steps) {
    pop = sample(pop,pop.size, replace = T)
    P[j] = sum(pop)/pop.size
  }
  # Step 3d: Add line to graph for this replicate
  lines(P)
}

normality.test = function(x) {
  ans = shapiro.test(x)
  if(ans$p.value > 0.05) {
    cat("The data are normally distributed: p = ",ans$p.value)
  } else {
    cat("The data are not normally distributed: p = ",ans$p.value)
  }
}

x = c(1,2,5,3,2,1) # enter your data into this variable
normality.test(x)

rgamma(100, shape = .5)

#####################################################################
# CHAPTER 13
#####################

# No code for this chapter

