#######################
# In this r scripts, we provide all the codes for generating the figures and the tables.
#######################
library(rempsyc)
library(rootSolve)
library(Rcpp)
library(BH)
library(bfpwr)
library(tikzDevice)  
library(colorspace)
sourceCpp("boost_noncentralt.cpp", cacheDir = "tmp_cache")
sourceCpp("pt.cpp", cacheDir = "tmp_cache")               
source("onesample.r")
source("twosample.r") 
####################### Figure 1 
# panel 1
# Specification
location = 0            #location parameter
scale = .707            #scaling parameter


# prior density for delta
delta = seq(-4,4,.001)
d_tstude = tstude(delta,location,scale,1)
d_norm = dnorm(delta,location,scale)
d_nlp = dnlp(delta,location,scale)


# Code for bottom panel of Figure 1:
t = seq(-20,20,.01)
df = 150
predict_cauchy = NA
predict_normal = NA
predict_nlp = NA

# Equation (7)  for marginalized likelihood of the predictive distributions
for (i in 1:length(t)){
  
  int <-function(delta){
    dnct(t[i],df,ncp = delta *sqrt(df+1))*tstude(delta,location, scale,1)}
  predict_cauchy[i]= integrate(int,lower = -Inf,upper = Inf)$value
  
  int <-function(delta){
    dnct(t[i],df,ncp = delta *sqrt(df+1))*dnorm(delta,location, scale)}
  predict_normal[i]= integrate(int,lower = -Inf,upper = Inf)$value
  
  int <-function(delta){
    dnct(t[i],df,ncp = delta *sqrt(df+1))*dnlp(delta,location,scale)}
  predict_nlp[i]= integrate(int,lower = -Inf,upper = Inf)$value
  
}



# Figure 1:
plot.name <- "Figure1"
# Part 1/2 - Save to .tex:
# Top panel:
tikz(paste0(plot.name, ".tex"), standAlone=TRUE, width = 6, height = 6) 
par(mfrow = c(2, 1), mar = c(4, 4, 1.5, 0))
plot(delta, d_tstude, 
     xlab = "", ylab = "", main = "", xaxt = "n", yaxt = "n", frame.plot = FALSE, 
     ylim = c(0, .6), 
     type = "l", lty = 1, lwd = 2)
axis(1, seq(-4, 4, 2), paste0("$", seq(-4, 4, 2), "$")) 
axis(2, seq(0, .6, .2), las = 1)
mtext("Effect size $\\delta$", 1, 2) 
mtext("Probability density", 2, 3)
mtext("\\itshape Prior distributions for $\\delta$ under $\\mathcal{H}_1$", 3, 0)
lines(delta, d_norm, lty = 2, lwd = 2)
lines(delta, d_nlp,  lty = 9, lwd = 2)
legend("topleft", 
       legend = c("Cauchy $(0, r = .707)$", "$\\mathcal{N}(0, \\sigma_\\delta^2 = 0.707^2)$", "Non-local $(0, \\tau = 0.707)$"), 
       lty = c(1, 2, 9), lwd = 2, y.intersp=1.1,
       bty = "n", cex = 1) 

# Bottom panel:
plot(t, predict_cauchy, 
     xlab = "", ylab= "", main = "", xaxt = "n", yaxt = "n", frame.plot = FALSE, 
     ylim = c(0, .05), # ylim=c(0,max(predict_normal)), 
     type = "l", lty = 1, lwd = 2)
axis(1, seq(-20, 20, 10), paste0("$", seq(-20, 20, 10), "$")) 
axis(2, seq(0, .05, .01), las = 1)
mtext("Test statistic $t$", 1, 2) 
mtext("Marginalized density", 2, 3)
mtext("\\itshape Prior predictive distribution of $t$-values", 3, 0)
lines(t, predict_normal, lty = 2, lwd = 2)
lines(t, predict_nlp,    lty = 9, lwd = 2)
dev.off()

# Part 2/2 - Create PNG (.tex -> .pdf -> .png -> clean up):
system(paste0("pdflatex ", plot.name, ".tex; 
       magick -density 300 ", plot.name, ".pdf ", plot.name, ".png; 
       rm *.aux; rm *.log; rm *.tex; rm *.pdf")
)


####################### Figure 2 
D = 10            # bound of compelling evidence 
df = 150          # degree of freedom
location = 0      # location parameter 
scale = .707      # scaling parameter
hypothesis = "!=" # the direction of the alternative 
model = "Cauchy"  # distribution 
target = .8       # targeted power 
#######
# t-values
t= seq(from = -5,to = 5, .01)
# BFs based on different t-values
BF10 = BF10_one_sample(t,df,model ,location,scale,dff,hypothesis)
#######
# density of predictive distribution under the null and alternative
Dh0 = dt(t,df=df)
Dh1 = NA
# Equation 7 
for (i in 1:length(t)){
  int <-function(delta){
    dnct(t[i],df,ncp = delta *sqrt(df+1))*tstude(delta,location,scale,1)}
  Dh1[i]= integrate(int,lower = -Inf,upper = Inf)$value
  
}

# finding t-values when BF10 = D or BF01 = D 
BF_D1 = BF_bound_10(D = D, df =df,model = model,location =location,scale=scale,dff = dff, hypothesis)
BF_D2 = BF_bound_01(D = D, df =df,model = model,location =location,scale=scale,dff = dff, hypothesis)

######################### shading the curve
# Find the indices for shading areas outside BF_D1 and inside BF_D2
outside_left_BF_D1  <- t <= BF_D1[1]
outside_right_BF_D1 <- t >= BF_D1[2]
inside_BF_D2        <- t >= BF_D2[1] & t <= BF_D2[2]



# Figure 2:
plot.name <- "Figure2"
# Part 1/2 - Save to .tex:
tikz(paste0(plot.name, ".tex"), standAlone=TRUE, width = 6, height = 6)
par(mfrow = c(2, 2),               # setting the plots to be 2 by 2 
    mar   = c(4, 4.5, 1.5, 0))     ### mar   = c(5.1, 4.1, 4.1, 2.1)) # setting the margin of the plot

# Top-left plot:
plot(t, log(1/BF10), 
     xlab = "", ylab = "", main = "", xaxt = "n", yaxt = "n", frame.plot = FALSE, 
     ylim = c(-9, 3), 
     type = "l", lty = 1, lwd = 2)
axis(1, c(-5, BF_D2, 5), c("$-5$", "", "", "5")) 
mtext(paste0("$", round(BF_D2[1], 3), "$"), 1, 1, at = -1, cex = .85) # Cheating on the position of labels, to separate them
mtext(paste0("$", round(BF_D2[2], 3), "$"), 1, 1, at = .8, cex = .85) 
axis(2, seq(-9, 3, 3), paste0("$", seq(-9, 3, 3), "$"), las = 1)
# mtext("Test statistic $t$", 1, 2.5, cex = .8) 
mtext("Logarithm of $BF_{01}$", 2, 3.5, cex = .8)
mtext(paste0("\\itshape $BF_{01}=10$ when $t=\\pm $", round(BF_D2[2], 3)), 3, 0, cex = .8)
segments(BF_D2, c(-10, -10), BF_D2, c(3, 3), lty = 3, col = "gray60")

# Top-right plot:
plot(t, Dh0, 
     xlab = "", ylab= "", main = "", xaxt = "n", yaxt = "n", frame.plot = FALSE, 
     ylim = c(0, .4), 
     type = "l", lty = 1, lwd = 2)
axis(1, c(-5, BF_D1[1], BF_D2, BF_D1[2], 5), c("$-5$", paste0("$", round(BF_D1[1], 3), "$"), "", "", paste0("$", round(BF_D1[2], 3), "$"), "5")) 
mtext(paste0("$", round(BF_D2[1], 3), "$"), 1, 1, at = -1, cex = .85) # Cheating on the position of labels, to separate them
mtext(paste0("$", round(BF_D2[2], 3), "$"), 1, 1, at = .8, cex = .85) 
axis(2, seq(0, .4, .2), las = 1)
# mtext("Test statistic $t$", 1, 2.5, cex = .8) 
mtext("Probability density", 2, 3.5, cex = .8)
mtext("\\itshape Prior predictive distribution under $\\mathcal{H}_0$", 3, 0, cex = .8)
# Shade the areas outside BF_D1 in red:
polygon(c(t[outside_left_BF_D1], rev(t[outside_left_BF_D1])), c(Dh0[outside_left_BF_D1], rep(0, sum(outside_left_BF_D1))), 
        col = "#FF00004C", border = NA)
polygon(c(t[outside_right_BF_D1], rev(t[outside_right_BF_D1])), c(Dh0[outside_right_BF_D1], rep(0, sum(outside_right_BF_D1))), 
        col = "#FF00004C", border = NA)
# Shade the area inside BF_D2 in blue:
polygon(c(t[inside_BF_D2], rev(t[inside_BF_D2])), c(Dh0[inside_BF_D2], rep(0, sum(inside_BF_D2))), 
        col = "#0000FF4C", border = NA)
# Draw vertical lines at the boundaries of BF_D1 and BF_D2:
segments(BF_D1, c(-1, -1), BF_D1, c(.4, .4), lty = 2, col = "gray60")   # Vertical lines for BF_D1
segments(BF_D2, c(-1, -1), BF_D2, c(.4, .4), lty = 3, col = "gray60")   # Vertical lines for BF_D2

# Bottom-left plot:
plot(t, log(BF10), 
     xlab = "", ylab = "", main = "", xaxt = "n", yaxt = "n", frame.plot = FALSE, 
     ylim = c(-3, 9), 
     type = "l", lty = 1, lwd = 2)
axis(1, c(-5, BF_D1, 5), c("$-5$", paste0("$", round(BF_D1, 3), "$"), "5")) 
axis(2, seq(-3, 9, 3), paste0("$", seq(-3, 9, 3), "$"), las = 1)
mtext("Test statistic $t$", 1, 2.5, cex = .8) 
mtext("Logarithm of $BF_{10}$", 2, 3.5, cex = .8)
mtext(paste0("\\itshape $BF_{10}=10$ when $t=\\pm $", round(BF_D1[2], 3)), 3, 0, cex = .8)
segments(BF_D1, c(-4, -4), BF_D1, c(9, 9), lty = 2, col = "gray60")

# Bottom-right:
plot(t, Dh1, 
     xlab = "", ylab= "", main = "", xaxt = "n", yaxt = "n", frame.plot = FALSE, 
     ylim = c(.020, .036), 
     type = "l", lty = 1, lwd = 2)
axis(1, c(-5, BF_D1[1], BF_D2, BF_D1[2], 5), c("$-5$", paste0("$", round(BF_D1[1], 3), "$"), "", "", paste0("$", round(BF_D1[2], 3), "$"), "5")) 
mtext(paste0("$", round(BF_D2[1], 3), "$"), 1, 1, at = -1, cex = .85) # Cheating on the position of labels, to separate them
mtext(paste0("$", round(BF_D2[2], 3), "$"), 1, 1, at = .8, cex = .85) 
axis(2, seq(.020, .036, .004), las = 1)
mtext("Test statistic $t$", 1, 2.5, cex = .8) 
mtext("Marginalized density", 2, 3.5, cex = .8)
mtext("\\itshape Prior predictive distribution under $\\mathcal{H}_1$", 3, 0, cex = .8)
# Shade the areas outside BF_D1 in red:
polygon(c(t[outside_left_BF_D1], rev(t[outside_left_BF_D1])), c(Dh1[outside_left_BF_D1], rep(0, sum(outside_left_BF_D1))), 
        col = "#FF00004C", border = NA)
polygon(c(t[outside_right_BF_D1], rev(t[outside_right_BF_D1])), c(Dh1[outside_right_BF_D1], rep(0, sum(outside_right_BF_D1))), 
        col = "#FF00004C", border = NA)
# Shade the area inside BF_D2 in blue:
polygon(c(t[inside_BF_D2], rev(t[inside_BF_D2])), c(Dh1[inside_BF_D2], rep(0, sum(inside_BF_D2))), 
        col = "#0000FF4C", border = NA)
# Draw vertical lines at the boundaries of BF_D1 and BF_D2:
segments(BF_D1, c(.015, .015), BF_D1, c(.036, .036), lty = 2, col = "gray60")   # Vertical lines for BF_D1
segments(BF_D2, c(.015, .015), BF_D2, c(.036, .036), lty = 3, col = "gray60")   # Vertical lines for BF_D2

dev.off()

# Part 2/2 - Create PNG (.tex -> .pdf -> .png -> clean up):
system(paste0("pdflatex ", plot.name, ".tex; 
       magick -density 300 ", plot.name, ".pdf ", plot.name, ".png; 
       rm *.aux; rm *.log; rm *.tex; rm *.pdf")
)



####################### Figure 4 
# loading the data for the plot
load("figure4.RData")
# Colors:
colors <- qualitative_hcl(3, palette = "Dark 3")

# New Figure:
plot.name <- "figure4"
# Part 1/2 - Save to .tex:
tikz(paste0(plot.name, ".tex"), standAlone = TRUE, width = 1.2*8, height = 1.2*5)
par(mfrow = c(1, 2), mar = c(4, 5, 1.5, .5))
# Left panel:
plot(NULL, xlim = c(0, 150), ylim = c(0, .061),
     xlab = "", ylab = "", main = "", xaxt = "n", yaxt = "n", frame.plot = FALSE)
abline(h = .05, lty = 3, col = "gray60")
lines(n1, alpha[, 1,1], col = colors[1], lwd = 2)
lines(n1, alpha[, 2,1], col = colors[2],  lwd = 2)
lines(n1, alpha[, 3,1], col = colors[3],   lwd = 2)
lines(n1, alpha[, 1,2], col = colors[1], lwd = 2, lty = 2)
lines(n1, alpha[, 2,2], col = colors[2],  lwd = 2, lty = 2)
lines(n1, alpha[, 3,2], col = colors[3],   lwd = 2, lty = 2)
axis(1, seq(0, 150, 30), paste0("$", seq(0, 150, 30), "$"), cex.axis = 1.4)
axis(2, c(.05, seq(0, .06, .02)), paste0("$", c(.05, seq(0, .06, .02)), "$"),
     las = 1,cex.axis = 1.4)
mtext("False Positive Evidence", 2, 3.5, cex = 1.4)
mtext("Sample size per group $N$", 1, 2.5, cex = 1.4)
mtext("$\\delta=0$", 3, 0, cex = 1.4)

# Right panel:
plot(NULL, xlim = c(0, 150), ylim = c(0, 1),
     xlab = "", ylab = "", main = "", xaxt = "n", yaxt = "n", frame.plot = FALSE)
lines(n1, f_power, lty = 3, col = "gray60",lwd = 2)
lines(n1, power[, 1,1], col = colors[1], lwd = 2)
lines(n1, power[, 2,1], col = colors[2], lwd = 2)
lines(n1, power[, 3,1], col = colors[3], lwd = 2)
lines(n1, power[, 1,2], col = colors[1],  lty = 2, lwd = 2)
lines(n1, power[, 2,2], col = colors[2], lty = 2, lwd = 2)
lines(n1, power[, 3,2], col = colors[3],  lty = 2, lwd = 2)
axis(1, seq(0, 150, 30), paste0("$", seq(0, 150, 30), "$"), cex.axis = 1.4)
axis(2, seq(0, 1, .2), paste0("$", seq(0, 1, .2), "$"), las = 1,cex.axis = 1.4)
mtext("True Positive Evidence", 2, 3, cex = 1.4)
mtext("Sample size per group $N$", 1, 2.5, cex = 1.4)
mtext("$\\delta=.64$", 3, 0, cex = 1.4)

# Legend:
legend(76, .4, title = "$BF_b$", legend = c(3, 10),
       lty = c(1, 2),  cex = 1.2, bty = "n", lwd = 2, seg.len = 3, title.adj = .2)
legend(76, .2, title = "Analysis prior", legend = c("Cauchy(0, .707)", "Normal(0, 1)", "NLP(0, .45)"),
       col = c(colors[1], colors[2], colors[3]), lty = 1, cex = 1.2, bty = "n", lwd = 2, seg.len = 3, title.adj = .14)
dev.off()

# Part 2/2 - Create PNG (.tex -> .pdf -> .png -> clean up):
system(paste0("pdflatex ", plot.name, ".tex;
       magick -density 300 ", plot.name, ".pdf ", plot.name, ".png;
       rm *.aux; rm *.log; rm *.tex; rm *.pdf")
)

################################ generating the data for the plots
# input
D = c(3,10)        # decision bound 
location = 0            # location parameter
hypothesis = ">"       # the direction of the hypotheses
scale = c(.707,.1,.45) # scaling parameter
n1 = seq(2,150,by = 2)  
r = 1
model = c("Cauchy","Normal","NLP")
alpha = array(NA, dim = c(length(n1),length(model),length(D)))
for (iii in 1:length(D)){
  for (ii in 1:length(model)){
    for (i in 1:length(n1)){
      
      t = BF_bound_10_two(D[iii] ,n1[i] ,r,model = model[ii] ,location ,scale = scale[ii] ,dff ,hypothesis )
      alpha[i,ii,iii] = false_positive_evidence_two (t,n1[i],r,model = "Point" ,location = .64 ,scale,dff, hypothesis)
      
    }}}

power = array(NA, dim = c(length(n1),length(model),length(D)))
for (iii in 1:length(D)){
  for (ii in 1:length(model)){
    for (i in 1:length(n1)){
      
      t = BF_bound_10_two(D[iii] ,n1[i] ,r,model = model[ii] ,location ,scale = scale[ii] ,dff ,hypothesis )
      power[i,ii,iii] = pro_compelling_BF_two(t,n1[i],r,model = "Point" ,location = .64 ,scale,dff, hypothesis)
      
    }}}

f_power = power.t.test(n = n1, delta = .64,sig.level = .05,type = "two.sample", alternative = "one.sided")$power
save(alpha,power,f_power,n1,  file ="new.RData")
####################### Figure 5  
##  loading the needed information for making the plot
load("probabilities.RData")

# Figure 5:
plot.name <- "Figure5"
# Part 1/2 - Save to .tex:
tikz(paste0(plot.name, ".tex"), standAlone=TRUE, width = 8, height = 5) 
par(mfrow = c(2, 3),               
    mar   = c(4, 5, 1.5, .5))
# i = 1 Cauchy prior
# i = 2 Normal prior
# i = 3 Non-local prior 

# Top row:
for (i in 1:3){
  plot(n1, alpha[ , 4, 1, i], 
       xlab = "", ylab = "", main = "", xaxt = "n", yaxt = "n", frame.plot = FALSE, 
       ylim = c(0, .06), 
       type = "l", lty = 1, lwd = 2)
  axis(1, seq(0, 500, 100),cex.axis = 1.4) 
  axis(2, seq(0, .06, .02), las = 1,cex.axis = 1.4)
  axis(2,.05,, las = 1,cex.axis = 1.4)
  abline( h = .05, lty = 2, col = "gray60")
  # mtext("Sample size per group $N$", 1, 2.5, cex = .8) 
  if (i == 1) mtext("False Positive Evidence", 2, 3, cex = .8)
  if (i == 1) mtext("\\itshape $\\mathcal{H}_1:\\delta\\sim$ Cauchy $(r)$", 3, 0, cex = 1)
  if (i == 2) mtext("\\itshape $\\mathcal{H}_1:\\delta\\sim\\mathcal{N}(0, \\sigma^2_\\delta)$", 3, 0, cex = 1)
  if (i == 3) mtext("\\itshape $\\mathcal{H}_1:\\delta\\sim$ Non-local $(0, \\tau)$", 3, 0, cex = 1)
  lines(n1, alpha[, 3,1,i], col = "blue")
  lines(n1, alpha[, 2,1,i], col = "red")
  lines(n1, alpha[, 1,1,i], col = "green")
  lines(n1, alpha[, 4,2,i], col = "black", lty = 2)
  lines(n1, alpha[, 3,2,i], col = "blue",  lty = 2)
  lines(n1, alpha[, 2,2,i], col = "red",   lty = 2)
  lines(n1, alpha[, 1,2,i], col = "green", lty = 2)
  if (i == 1) {
    legend(290, .05, title = "Scaling parameter", legend = c(1,.707,.5,.1),
           col = c("black","blue","red","green"), lty = 1, cex = 1.2, bty = "n", 
           lwd = 2, seg.len = 3, title.adj = -0.4,
           y.intersp = .7,x.intersp = 1)
    legend(305, .025, title = "$BF_b$", legend = c(3, 10),
           lty = c(1, 2),  cex = 1.2, bty = "n", lwd = 2, seg.len = 3, title.adj = .4,
           y.intersp = .7,x.intersp = 1)
  }
}

# Bottom row:
for (i in 1:3){
  plot(n1, power[ , 4, 1, i], 
       xlab = "", ylab = "", main = "", xaxt = "n", yaxt = "n", frame.plot = FALSE, 
       ylim = c(0, 1), 
       type = "l", lty = 1, lwd = 2)
  axis(1, seq(0, 500, 100),cex.axis = 1.4) 
  axis(2, seq(0, 1, .2), las = 1,cex.axis = 1.4)
  abline( h = .8, lty = 2, col = "gray60")
  mtext("Sample size per group $N$", 1, 2.5, cex = 1) 
  if (i == 1) mtext("True Positive Evidence", 2, 3, cex = .8)
  # if (i == 1) mtext("\\itshape $\\mathcal{H}_1:\\delta\\sim$ Cauchy $(r)$", 3, 0, cex = .8)
  # if (i == 2) mtext("\\itshape $\\mathcal{H}_1:\\delta\\sim\\mathcal{N}(0, \\sigma^2_\\delta)$", 3, 0, cex = .8)
  # if (i == 3) mtext("\\itshape $\\mathcal{H}_1:\\delta\\sim$ Non-local $(0, \\tau)$", 3, 0, cex = .8)
  lines(n1, power[, 3,1,i], col = "blue")
  lines(n1, power[, 2,1,i], col = "red")
  lines(n1, power[, 1,1,i], col = "green")
  lines(n1, power[, 4,2,i], col = "black", lty = 2)
  lines(n1, power[, 3,2,i], col = "blue",  lty = 2)
  lines(n1, power[, 2,2,i], col = "red",   lty = 2)
  lines(n1, power[, 1,2,i], col = "green", lty = 2)

}

dev.off()

# Part 2/2 - Create PNG (.tex -> .pdf -> .png -> clean up):
system(paste0("pdflatex ", plot.name, ".tex; 
       magick -density 300 ", plot.name, ".pdf ", plot.name, ".png; 
       rm *.aux; rm *.log; rm *.tex; rm *.pdf")
)

#####################colorblind friendly version
# Figure 5:
plot.name <- "Figure5colorfriendly"

colors <- qualitative_hcl(4, palette = "Dark 2")

# Part 1/2 - Save to .tex:
tikz(paste0(plot.name, ".tex"), standAlone=TRUE, width = 8, height = 5) 
par(mfrow = c(2, 3),               
    mar   = c(4, 5, 1.5, .5))
# i = 1 Cauchy prior
# i = 2 Normal prior
# i = 3 Non-local prior 

# Top row:
for (i in 1:3){
  plot(n1, alpha[ , 4, 1, i], 
       xlab = "", ylab = "", main = "", xaxt = "n", yaxt = "n", frame.plot = FALSE, 
       ylim = c(0, .06), 
       type = "l", lty = 1, lwd = 2,col = colors[1])
  axis(1, seq(0, 500, 100),cex.axis = 1.4) 
  axis(2, seq(0, .06, .02), las = 1,cex.axis = 1.4)
  axis(2,.05,, las = 1,cex.axis = 1.4)
  abline( h = .05, lty = 2, col = "gray60")
  # mtext("Sample size per group $N$", 1, 2.5, cex = .8) 
  if (i == 1) mtext("False Positive Evidence", 2, 3, cex = .8)
  if (i == 1) mtext("\\itshape $\\mathcal{H}_1:\\delta\\sim$ Cauchy $(r)$", 3, 0, cex = 1)
  if (i == 2) mtext("\\itshape $\\mathcal{H}_1:\\delta\\sim\\mathcal{N}(0, \\sigma^2_\\delta)$", 3, 0, cex = 1)
  if (i == 3) mtext("\\itshape $\\mathcal{H}_1:\\delta\\sim$ Non-local $(0, \\tau)$", 3, 0, cex = 1)
  lines(n1, alpha[, 3,1,i], col = colors[2])
  lines(n1, alpha[, 2,1,i], col = colors[3])
  lines(n1, alpha[, 1,1,i], col = colors[4])
  lines(n1, alpha[, 4,2,i], col = colors[1], lty = 2)
  lines(n1, alpha[, 3,2,i], col = colors[2],  lty = 2)
  lines(n1, alpha[, 2,2,i], col = colors[3],   lty = 2)
  lines(n1, alpha[, 1,2,i], col = colors[4], lty = 2)
  if (i == 1) {
    legend(290, .05, title = "Scaling parameter", legend = c(1,.707,.5,.1),
           col = colors, lty = 1, cex = 1.2, bty = "n", 
           lwd = 2, seg.len = 3, title.adj = -0.4,
           y.intersp = .7,x.intersp = 1)
    legend(305, .025, title = "$BF_b$", legend = c(3, 10),
           lty = c(1, 2),  cex = 1.2, bty = "n", lwd = 2, seg.len = 3, title.adj = .4,
           y.intersp = .7,x.intersp = 1)
  }
  

  
}


# Bottom row:
for (i in 1:3){
  plot(n1, power[ , 4, 1, i], 
       xlab = "", ylab = "", main = "", xaxt = "n", yaxt = "n", frame.plot = FALSE, 
       ylim = c(0, 1), 
       type = "l", lty = 1, lwd = 2,col = colors[1])
  axis(1, seq(0, 500, 100),cex.axis = 1.4) 
  axis(2, seq(0, 1, .2), las = 1,cex.axis = 1.4)
  abline( h = .8, lty = 2, col = "gray60")
  mtext("Sample size per group $N$", 1, 2.5, cex = 1) 
  if (i == 1) mtext("True Positive Evidence", 2, 3, cex = .8)
  # if (i == 1) mtext("\\itshape $\\mathcal{H}_1:\\delta\\sim$ Cauchy $(r)$", 3, 0, cex = .8)
  # if (i == 2) mtext("\\itshape $\\mathcal{H}_1:\\delta\\sim\\mathcal{N}(0, \\sigma^2_\\delta)$", 3, 0, cex = .8)
  # if (i == 3) mtext("\\itshape $\\mathcal{H}_1:\\delta\\sim$ Non-local $(0, \\tau)$", 3, 0, cex = .8)
  lines(n1, power[, 3,1,i], col = colors[2])
  lines(n1, power[, 2,1,i], col = colors[3])
  lines(n1, power[, 1,1,i], col = colors[4])
  lines(n1, power[, 4,2,i], col = colors[1], lty = 2)
  lines(n1, power[, 3,2,i], col = colors[2],  lty = 2)
  lines(n1, power[, 2,2,i], col = colors[3],   lty = 2)
  lines(n1, power[, 1,2,i], col = colors[4], lty = 2)
}

dev.off()

# Part 2/2 - Create PNG (.tex -> .pdf -> .png -> clean up):
system(paste0("pdflatex ", plot.name, ".tex; 
       magick -density 300 ", plot.name, ".pdf ", plot.name, ".png; 
       rm *.aux; rm *.log; rm *.tex; rm *.pdf")
)

## the information from "probabilities.RData" is generated using the following codes:
# input
D = c(10,3)             # decision bound 
location = 0            # location parameter
hypothesis = ">"       # the direction of the hypotheses
scale = c(.1,.5,.707,1) # scaling parameter
n1 = seq(2,500,by = 5)  
r = 1
df = n1+n2*r -2
#############################################
dff = 1
model = c("Cauchy","Normal","NLP")
title = c(bquote("H"[1]~":"~delta~"~ Cauchy(r)"),
          bquote("H"[1]~":"~delta~"~ Normal(0,"~sigma[delta]^2~")"),
          bquote("H"[1]~":"~delta~"~ Non-local(0,"~tau~")"))
alpha = array(NA, dim = c(length(n1),length(scale),length(D),length(model)))
for(iv in 1:length(model)){
  for(iii in 1:length(D)){
    for(ii in 1:length(scale)){
      for ( i in 1:length(n1)){
        t = BF_bound_10_two(D[iii] ,n1[i] ,r,model[iv] ,location ,scale[ii] ,dff ,hypothesis )
        alpha[i,ii,iii,iv] = false_positive_evidence_two (t,n1[i],r,model[iv] ,location ,scale[ii],dff, hypothesis)
      }}}}

power = array(NA, dim = c(length(n1),length(scale),length(D),length(model)))
for(iv in 1:length(model)){
  for(iii in 1:length(D)){
    for(ii in 1:length(scale)){
      for ( i in 1:length(n1)){
        bound = BF_bound_10_two(D[iii] ,n1[i] ,r,model[iv] ,location ,scale[ii] ,dff ,hypothesis )
        power[i,ii,iii,iv] = pro_compelling_BF_two(bound,n1[i],r,model[iv] ,location ,scale[ii],dff , hypothesis )
      }}}}
save(alpha,power,title,n1,  file ="probabilities.RData")


####################### Table 1 & 2
# Specification 
# Mode:
Mode = 1 # 1: A priori BFDA   0: post hoc BFDA

# Specification of design and analysis priors

# Analysis prior 
BF_bound =  10
Power = seq(.5,.95,.05)        # desired probability of compelling evidence for the alternative given a design prior
Model = "t-distribution"    # Other prior distributions: "t-distribution"  , "Normal"  , "NLP"
location = 0      # location parameter
scale = .707       # scaling parameter or the standard deviation for the normal distribution
df = 3          # degree of freedom for t-distribution
hypothesis = hypothesis_d =  "!="  #"!=": two-sided   ">": effect size greater than zero "<": effect size smaller than zero
r = 1             # Ratio of sample size in group 2 over the one in group 1 for independent t-test
Design_prior = 0  # 1: Analysis prior and design prior are the same
# 0: Analysis prior and design prior are different

# Design priors
Model_d = "Normal"   # Alternative models: "t-distribution" "Normal" "NLP"
location_d = 0   # location parameter
scale_d = 1  # scaling parameter or the standard deviation for normal model
df_d = 1          # degree of freedom for t-distribution
####################### Table 1 
Required_N = data.frame(matrix(ncol = length(Power), nrow = 2))
colnames(Required_N) = Power
for (i in 1:length(Power)){
  Required_N[1,i] = ntbf01(k = 1/BF_bound, power = Power[i], pscale = scale,pdf =df,dpm = location_d, dpsd = scale_d,type = "two.sample",  alternative = "two.sided")
  Required_N[2,i] =  ceiling(N1_finder_two(BF_bound,r,Power[i],Model,location,scale,df, hypothesis,
                                           Model_d,location_d,scale_d,dff_d, hypothesis_d,Design_prior))
}

Required_N
nice_table(round(Required_N))
####################### Table 2
n  = seq(4,20,2)
n = c(n,1000)
Power_methods = array(NA,c(2,length(n)))
colnames(Power_methods)=n
for (i in 1:length(n)){  
  Power_methods[1,i] = powertbf01(k = 1/BF_bound, n[i],plocation= location, pscale = scale, pdf = df,dpm = location_d, dps = scale_d, alternative = "two.sided")$power
  t = BF_bound_10_two(BF_bound,n[i],r,Model ,location ,scale ,df ,hypothesis )
  Power_methods[2,i]=pro_compelling_BF_two(t,n[i],r ,Model_d ,location_d ,scale_d ,df_d ,hypothesis_d )
}

nice_table(Power_methods*100)
