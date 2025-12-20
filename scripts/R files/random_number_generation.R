### to play with random numbers (RN)
#
# by Alexander Skupin 2023/10/25
#

### pseudo random number generation (PRNG)

nn <- 4000 # number of sequence

# starting point
xs = 1/123 #.5 #1/9# .51245422624
#multiplier 
mul = 2.132 # 2
#drift 
dri =  0.212134 # 0

#vector for random numbers and initiation 
rannum <- c(rep(0,nn))
rannum[1] <- xs

#loop to generate RN
for (i in c(2:nn)){
  rannum[i] <- (mul * rannum[i-1] + dri) %% 1
}
plot (rannum)

plot (rannum[c(1:nn-1)],rannum[c(2:nn)]) # plot dependencies

hist(rannum, breaks = 60) # histogram

disr <- hist(rannum, breaks=250) # properties of histogram

# to check for properties
mean(rannum) # average (to be 0.5 - see lecture)

sd(rannum)^2 # variance (to be 1/12 - see lecture)
# calculate entropy and average
ent = 0
danv = 0
for (i in c(1:length(disr$counts))){
  ptemp <- disr$counts[i]/nn
  ent = ent - ptemp  * log2(ptemp)
  danv = danv + ptemp *disr$mids[i] # the mids are the x
}

print(ent) # depends on binning, i.e breaks of the R function; 
print(danv)
#How does this relate to the analytical result for uniform distribution S = log(b-a) = 0?

# differential enrtopy (continous description) has to consider bin withs

# calculate differential entropy
nbis = 500 # 4500
ddisr <- hist(rannum, breaks=nbis) # properties of histogram
dent = 0
for (i in c(1:length(ddisr$counts))){
  ptemp <- ddisr$counts[i]/nn * 1/nbis
  if(ptemp >0){dent = dent - ptemp * log2(ptemp)}
}

print(dent) # depends on binning, i.e breaks of the R function; 


# exp-distribution by inverse sampling approach
lambda <- 1 #0.1
expran <- - lambda * log(rannum) # inverse function of exponential
expden <- hist(expran, breaks=500, freq=F)
sum(expden$counts)/nn # check if normailzed

sum ( expden$counts * expden$mids)/(nn) # average by distribution
mean(expran) # average of numbers
sd(expran) # standard deviation

#differential entropy - should be 1 - ln(lambda)?
# nebis = 100 # 4500
# dedisr <- hist(expran, breaks=nebis, freq=F) # properties of histogram
# deent = 0
# for (i in c(1:length(dedisr$counts))){
#   ptemp <- dedisr$counts[i]/nn  * 1/nebis
#   if(ptemp >0){deent = deent - ptemp * log2(ptemp )}
# }

# print(deent) # dep

#generate uniform random number by R FUNCTIONS
runif(10)

#plot random numbers
plot(runif(1000))

mean(runif(100000))

# in vector
rn <- runif(1000)

plot( rn[c(1:999)] , rn[c(2:1000)])

rnhist <- hist(rn,breaks=30)

#normalized histogram
sum((rnhist$counts/1000)) # check if normalized
barplot(rnhist$counts/1000) # plot normalized histogram

#density plot
plot(density(rn))

#binning depends on number of replicates - N(bin) = sqrt N(repl)
rnhist <- hist(rn, breaks=100)
rnhist <- hist(rn, breaks=30)

rn2 <- 2*runif(1000)
plot(rn[c(1:100)],ylim=c(0,2))
points(rn2, col="red")

plot(density(rn),xlim=c(0,3))
lines(density(rn2), col="red")

# sum of two uniform random numbers (like 2 dices)
rnsum <- rn+rn2
hist(rnsum)
plot(density(rnsum),ylim=c(0,2))

plot(density(rn),ylim=c(0,2),xlim=c(0,4))
lines(density(rn2), col="red")
lines(density(rnsum), col="blue")

#for generating exponential RN (inverse sampling approach)
hist(-log(runif(100000)))
plot(density(-log(runif(1000000))))


# as vector
rne <- -log(runif(1000))

rnehist <- hist(rne)
sum((rnehist$counts/1000)) # check if normalized
barplot(rnehist$counts/1000) # plot normalized histogram


# Gaussian random numbers

rng <- rnorm(10000, mean=10, sd=2) # with mean 10 and standard deviation 2

plot(rng[c(1:1000)])
rnghist <- hist(rng)
sum((rnghist$counts/10000)) # check if normalized
barplot(rnghist$counts/10000) # plot normalized histogram

# for significance tests

wilcox.test(rng,rn)
t.test(rng,rne) # wrong since requires normal distribution - to be checked by chi^2

hist(rn)
hist(rng)
hist(rne)
#set starting point (seed) for random number generator) - good for reproducability
set.seed(2)
rng2 <- rnorm(1000, mean=9.5, sd=2)

# to plot 
plot(rng2[c(1:100)])
points(rng, col="red")

#for density plots
plot(density(rng))
lines(density(rng2), col="red")

wilcox.test(rng,rng2)
wilcox.test(rng[1:120],rng2[1:120]) # significance depends on sample size!

t.test(rng,rng2)
t.test(rng[1:120],rng2[1:120])