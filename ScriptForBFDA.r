###########################################################
#In this r script, we provide the codes for conducting BFDA.
###########################################################
# loading packages and needed functions
library(rootSolve)
library(Rcpp)
library(BH)

sourceCpp("boost_noncentralt.cpp")
sourceCpp("pt.cpp")
source("onesample.r")
source("twosample.r")
#################################################
# Specification 
# Mode:
Mode = 1 # 1: A priori BFDA   0: BFDA for  a fixed sample size

#for post hoc BFDA
N1 = 51  # sample size for group 1 or total sample size for one sample and paired t-test
N2 = 51  # sample size for group 2


# Specification of design and analysis priors

# Analysis prior 
BF_bound =  3
Power = .8        # desired probability of compelling evidence for the alternative given a design prior
Models = "NLP"    # Prior distributions: "t-distribution"  , "Normal"  , "NLP" , "Cauchy"
location = 0      # location parameter
scale = .45       # scaling parameter or the standard deviation for the normal distribution
df = 1            # degree of freedom for t-distribution
hypothesis = ">"  #"!=": two-sided   ">": effect size greater than zero "<": effect size smaller than zero
r = 1             # Ratio of sample size in group 2 over the one in group 1 for independent t-test
Design_prior = 1  # 1: Analysis prior and design prior are the same
                  # 0: Analysis prior and design prior are different

# Design priors
Model_d = "NLP"   # Prior distributions: "t-distribution"  , "Normal"  , "NLP" , "Cauchy"
location_d = 0    # location parameter
scale_d = .3      # scaling parameter or the standard deviation for normal distribution
df_d = 1          # degree of freedom for t-distribution
hypothesis_d = hypothesis 

# For independent t-test, the BFDA can be done by:
print (Table_two (BF_bound,r,Power,
                         Models,location,scale,dff, hypothesis,Model_d,
                         location_d,scale_d,dff_d, hypothesis_d,Design_prior,N1,N2,Mode))


# For one-sample or paired t test, the BFDA can be done by:
print(Table(BF_bound,Power,
                    Models,location,scale,dff, hypothesis,Model_d,
                    location_d,scale_d,dff_d, hypothesis_d,Design_prior,
                    N1,Mode))

