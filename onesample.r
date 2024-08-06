# likelihood of non-local prior
dnlp <-function(delta,mu,ta){
  ((delta-mu)^2)/(sqrt(2*pi)*ta^3)*exp(-((delta-mu)^2)/(2*ta^2))
}
# likelihood of informed t prior
tstude <- function(t, location = 0, scale = sqrt(2)/2, df = 1) {
  gamma((df+1)/2) * ((df+((t-location)/scale)^2)/df)^(-((df+1)/2)) / (scale*sqrt(df*pi)*gamma(df/2))
  
  #dnct((t-location)/scale,df,ncp = 0)/scale
}


# likelihood of t under the null 
ml_H0 <-function(t,df){
  dnct(t,df,ncp = 0)
}

# marginal likelihood under the alternative
ml_H1_one_sample <-function(t,df,model ,location,scale,dff , hypothesis ){

  bound  <- switch(hypothesis,
                   ">" = c(a = 0, b = Inf),
                   "<" = c(a = -Inf, b = 0),
                   "!=" = c(a = -Inf, b = Inf)
  )
  x <- numeric(length(t))
  
  for(i in 1:length(t)){
    
    
    normalization  <- switch(model,
                             "Cauchy" = integrate(function(delta) tstude(delta,location, scale,1),lower = bound[1],upper = bound[2])$value,
                             "Normal" = integrate(function(delta)dnorm(delta,location,scale),lower = bound[1],upper = bound[2])$value,
                             "NLP" = integrate(function(delta)dnlp(delta,location,scale),lower = bound[1],upper = bound[2])$value,
                             "t-distribution" = integrate(function(delta)tstude(delta,location,scale,dff),lower = bound[1],upper = bound[2])$value
    )
    int  <- switch(model,
                   "Cauchy" = function(delta){
                     dnct(t[i],df,ncp = delta *sqrt(df+1))* tstude(delta,location, scale,1)/normalization},
                   
                   "Normal" = function(delta){
                     dnct(t[i],df,ncp = delta *sqrt(df+1))*dnorm(delta,location,scale)/normalization},
                   
                   "NLP" = int <-function(delta){
                     dnct(t[i],df,ncp = delta *sqrt(df+1))*dnlp(delta,location,scale)/normalization},
                   
                   "t-distribution" = int <-function(delta){
                     dnct(t[i],df,ncp = delta *sqrt(df+1))*tstude(delta,location, scale,dff)/normalization}
    )

    error = 1e-8
    if (model == "NLP" & scale <.3 ){
      error = 1e-14
    }
    x[i]= integrate(int,lower = bound[1],upper = bound[2], rel.tol=error,stop.on.error = F)$value
    
  }
  return(x)
}



# the Bayes Factor
BF10_one_sample <- function(t,df,model ,location,scale,dff , hypothesis ){
  ml_H1_one_sample(t=t,df=df,model=model ,location=location ,scale=scale,dff=dff, hypothesis) /ml_H0(t,df)
}

# finding the t that correspond to BF10=D
BF_bound_10 <-function(D, df,model ,location ,scale,dff , hypothesis){
  y <- numeric(0)
  Bound_finding <-function(t){
    BF10_one_sample(t,df=df,model=model,location=location,scale=scale,dff=dff, hypothesis =hypothesis )- D
  }
  
  if (hypothesis=="!="){
    x <- tryCatch({
      uniroot(Bound_finding, lower = -8, upper = 0)$root
    }, error = function(e) {
      if (grepl("f\\(\\) values at end points not of opposite sign", e$message)) {
        return(numeric(0)) # Return numeric(0) if the specific error is encountered
      } else {
        stop(e) # Rethrow the error if it's a different error
      }
    })
    y <- tryCatch({
      uniroot(Bound_finding, lower = 0, upper = 8)$root
    }, error = function(e) {
      if (grepl("f\\(\\) values at end points not of opposite sign", e$message)) {
        return(numeric(0)) # Return numeric(0) if the specific error is encountered
      } else {
        stop(e) # Rethrow the error if it's a different error
      }
    })
    
    
  } else {
    x <- tryCatch({
      uniroot.all(Bound_finding, lower = -8, upper = 8)
    }, error = function(e) {
      if (grepl("f\\(\\) values at end points not of opposite sign", e$message)) {
        return(numeric(0)) # Return numeric(0) if the specific error is encountered
      } else {
        stop(e) # Rethrow the error if it's a different error
      }
    })
    
    
    
  }
  
  if (length(x) ==0){
    x = "no bound is found"
    return(x)
  } 
  if (length(y) == 1){
    x = cbind(x,y)
  }
  BF = BF10_one_sample(x,df,model,location,scale,dff,hypothesis )
  x = x[round(BF,2)==round(D,2)]
  if (length(x) ==0){
    x = "no bound is found"
    return(x)
  } 
  return(x)
}


# finding the t that correspond to BF01=D
BF_bound_01 <-function(D , df,model ,location ,scale,dff , hypothesis){
  
  Bound_finding <-function(t){
    1/BF10_one_sample(t,df=df,model=model,location=location,scale=scale,dff=dff, hypothesis =hypothesis )-D
  }
  
  x = uniroot.all(Bound_finding, lower = -8,upper = 8)
  if (length(x) == 0 ){
    x = "bound cannot be found"
    return(x)
  }
  BF = 1/BF10_one_sample(x ,df=df,model=model,location=location,scale=scale,dff=dff, hypothesis =hypothesis )
  x = x[round(BF,1)== round(D,1)]
  
  return(x)
}

# p(BF01>D|H0)
pro_compelling_BF_null <- function(t , df,model,location ,scale,dff , hypothesis){

  if (any(t == "bound cannot be found")){
    return(t)
  }
  if (length(t)>1){
    pro = pt(max(t),df) - pt(min(t),df)
  }else{
    if (t >0){
      pro = pt(t,df)
    } else {
      pro = 1-pt(t,df)
    }
  }
  return(pro)
}

# p(BF10>D|H1)
pro_compelling_BF<-function(t,df,model ,location ,scale,dff , hypothesis ){
  
  if (any(t =="no bound is found" | length(t)==0)){
    t=0
    return(t)
  }
  
  bound  <- switch(hypothesis,
                   ">" = c(a = 0, b = Inf),
                   "<" = c(a = -Inf, b = 0),
                   "!=" = c(a = -Inf, b = Inf)
  )
  x = NULL

  normalization  <- switch(model,
                           "Cauchy" = integrate(function(delta) tstude(delta,location,scale,1),lower = bound[1],upper = bound[2])$value,
                           "Normal" = integrate(function(delta)dnorm(delta,location,scale),lower = bound[1],upper = bound[2])$value,
                           "NLP" = integrate(function(delta)dnlp(delta,location,scale),lower = bound[1],upper = bound[2])$value,
                           "t-distribution" = integrate(function(delta)tstude(delta,location,scale,dff),lower = bound[1],upper = bound[2])$value
  )
  if (length(t) > 1){
    
    int  <- switch(model,
                   "Cauchy" = function(delta){
                     pro1 = pnct(t[t>0],df,ncp = delta *sqrt(df+1),lower  = F)
                     pro2 = pnct(t[t<0],df,ncp = delta *sqrt(df+1),lower  = T)
                     (pro1+pro2)* tstude(delta,location,scale,1)/normalization },
                   
                   "Normal" = function(delta){
                     pro1 = pnct(t[t>0],df,ncp = delta *sqrt(df+1),lower  = F)
                     pro2 = pnct(t[t<0],df,ncp = delta *sqrt(df+1),lower  = T)
                     
                     (pro1+pro2)* dnorm(delta,location,scale)/normalization   },
                   
                   "NLP" = function(delta){
                     pro1 = pnct(t[t>0],df,ncp = delta *sqrt(df+1),lower  = F)
                     pro2 = pnct(t[t<0],df,ncp = delta *sqrt(df+1),lower  = T)
                     
                     (pro1+pro2)* dnlp(delta,location,scale)/normalization   } ,
                   
                   "t-distribution" = function(delta){
                     pro1 = pnct(t[t>0],df,ncp = delta *sqrt(df+1),lower  = F)
                     pro2 = pnct(t[t<0],df,ncp = delta *sqrt(df+1),lower  = T)
                     
                     (pro1+pro2)* tstude(delta,location,scale,dff)/normalization   }
    )
    
    
    
  }  else if (t >= 0) { 
    int  <- switch(model,
                   "Cauchy" = function(delta){
                     pro1 = pnct(t,df,ncp = delta *sqrt(df+1),lower  = F)
                     (pro1)*tstude(delta,location,scale,1)/normalization  },
                   
                   "Normal" = function(delta){
                     pro1 = pnct(t,df,ncp = delta *sqrt(df+1),lower  = F)
                     (pro1)* dnorm(delta,location,scale)/normalization   },
                   
                   "NLP" = function(delta){
                     pro1 = pnct(t,df,ncp = delta *sqrt(df+1),lower  = F)
                     (pro1)* dnlp(delta,location,scale)/normalization   },
                   
                   "t-distribution" = function(delta){
                     pro1 = pnct(t,df,ncp = delta *sqrt(df+1),lower  = F)
                     (pro1)* tstude(delta,location,scale,dff)/normalization}
    )
    
  } else if (t < 0){
    int  <- switch(model,
                   "Cauchy" = function(delta){
                     pro2 = pnct(t,df,ncp = delta *sqrt(df+1),lower  = T)
                     (pro2)* tstude(delta,location,scale,1)/normalization},
                   
                   
                   "Normal" = function(delta){
                     pro2 = pnct(t,df,ncp = delta *sqrt(df+1),lower  = T)
                     (pro2)* dnorm(delta,location,scale)/normalization   },
                   
                   "NLP" = function(delta){
                     pro2 = pnct(t,df,ncp = delta *sqrt(df+1),lower  = T)
                     (pro2)* dnlp(delta,location,scale)/normalization},
                   
                   "t-distribution" = int <-function(delta){
                     pro2 = pnct(t,df,ncp = delta *sqrt(df+1),lower  = T)
                     (pro2)* tstude(delta,location,scale,dff)/normalization   }
    )
    
  }
  if (scale >.3){
    error = .Machine$double.eps^0.25
  }else {
    error = 1e-8
  }
  
  if (model == "NLP" & scale <.3 ){
    error = 1e-14
  }
  x = integrate(int,lower = bound[1],upper = bound[2], rel.tol = error,stop.on.error=FALSE)$value
  
  return(x) 
  
} 

# p(BF01>D|H1)
false_negative_BF<-function(t,df,model ,location ,scale,dff , hypothesis ){
  
  if (any(t == "bound cannot be found")){
    return(t)
  }
  
  bound  <- switch(hypothesis,
                   ">" = c(a = 0, b = Inf),
                   "<" = c(a = -Inf, b = 0),
                   "!=" = c(a = -Inf, b = Inf)
  )
  x = NULL
  t <- as.numeric(t)
  if (length(t) == 4) {  # Corrected condition check
    t = t[2:3]
  }
  
  normalization  <- switch(model,
                           "Cauchy" = integrate(function(delta) tstude(delta,location,scale,1),lower = bound[1],upper = bound[2])$value,
                           "Normal" = integrate(function(delta)dnorm(delta,location,scale),lower = bound[1],upper = bound[2])$value,
                           "NLP" = integrate(function(delta)dnlp(delta,location,scale),lower = bound[1],upper = bound[2])$value,
                           "t-distribution" = integrate(function(delta)tstude(delta,location,scale,dff),lower = bound[1],upper = bound[2])$value
  )
  if (length(t)>1){
    
    int  <- switch(model,
                   "Cauchy" = function(delta){
                     pro1 = pnct(max(t),df,ncp = delta *sqrt(df+1),lower  = T)
                     pro2 = pnct(min(t),df,ncp = delta *sqrt(df+1),lower  = T)
                     (pro1-pro2)* tstude(delta,location,scale,1)/normalization },
                   
                   "Normal" = function(delta){
                     pro1 = pnct(max(t),df,ncp = delta *sqrt(df+1),lower  = T)
                     pro2 = pnct(min(t),df,ncp = delta *sqrt(df+1),lower  = T)
                     
                     (pro1-pro2)* dnorm(delta,location,scale)/normalization   },
                   
                   "NLP" = function(delta){
                     pro1 = pnct(max(t),df,ncp = delta *sqrt(df+1),lower  = T)
                     pro2 = pnct(min(t),df,ncp = delta *sqrt(df+1),lower  = T)
                     
                     (pro1-pro2)* dnlp(delta,location,scale)/normalization   } ,
                   
                   "t-distribution" = function(delta){
                     pro1 = pnct(max(t),df,ncp = delta *sqrt(df+1),lower  = T)
                     pro2 = pnct(min(t),df,ncp = delta *sqrt(df+1),lower  = T)
                     
                     (pro1-pro2)* tstude(delta,location,scale,dff)/normalization   }
    )
    
  } else if (t>0) { 
    int  <- switch(model,
                   "Cauchy" = function(delta){
                     pro1 = pnct(t,df,ncp = delta *sqrt(df+1),lower  = T)
                     (pro1)*tstude(delta,location,scale,1)/normalization  },
                   
                   "Normal" = function(delta){
                     pro1 = pnct(t,df,ncp = delta *sqrt(df+1),lower  = T)
                     (pro1)* dnorm(delta,location,scale)/normalization   },
                   
                   "NLP" = function(delta){
                     pro1 = pnct(t,df,ncp = delta *sqrt(df+1),lower  = T)
                     (pro1)* dnlp(delta,location,scale)/normalization   },
                   
                   "t-distribution" = function(delta){
                     pro1 = pnct(t,df,ncp = delta *sqrt(df+1),lower  = T)
                     (pro1)* tstude(delta,location,scale,dff)/normalization}
    )
    
  } else if (t<0){
    int  <- switch(model,
                   "Cauchy" = function(delta){
                     pro2 = pnct(t,df,ncp = delta *sqrt(df+1),lower  = F)
                     (pro2)* tstude(delta,location,scale,1)/normalization},
                   
                   
                   "Normal" = function(delta){
                     pro2 = pnct(t,df,ncp = delta *sqrt(df+1),lower  = F)
                     (pro2)* dnorm(delta,location,scale)/normalization   },
                   
                   "NLP" = function(delta){
                     pro2 = pnct(t,df,ncp = delta *sqrt(df+1),lower  = F)
                     (pro2)* dnlp(delta,location,scale)/normalization},
                   
                   "t-distribution" = int <-function(delta){
                     pro2 = pnct(t,df,ncp = delta *sqrt(df+1),lower  = F)
                     (pro2)* tstude(delta,location,scale,dff)/normalization   }
    )
    
  }
  x = integrate(int,lower = bound[1],upper = bound[2],stop.on.error = FALSE)$value
  
  if (length(t) == 1 ){
    if (t>0 & hypothesis == "<"){
      x = 1-x
    }
    if (t<0 & hypothesis == ">") {
      x = 1-x
    }}
  return(x) 
} 

# p(BF10>D|H0)
false_positive_evidence <- function(t,df,model ,location ,scale,dff, hypothesis){
  if (hypothesis=="!="){
    
    if (any(t[1] == "no bound is found"|t[2] == "no bound is found")){
    return(0)
  }
  } else if (t == "no bound is found"){
    return(0)
  }
  
  if (length(t) == 4) {  # Corrected condition check
    t = t[2:3]
  }
  if (length(t)>1){
    pro = pt(max(t),df=df,lower.tail = F) +pt(min(t),df=df,lower.tail = T)
  }else if (t>0){
    pro = pt(t[t>0],df=df,lower.tail = F)
  }else if (t<0){
    pro = pt(t[t<0],df=df,lower.tail = T)
  }
  return(pro)
  
}

# Finding the degree of freedom that ensure p(BF10>D|H1) > targeted probability

N_finder<-function(D,target,model,location,scale,dff=1, hypothesis ,
                   model_d,location_d,scale_d,dff_d=1, hypothesis_d,de_an_prior ){
 
  lo = 2
  t = BF_bound_10(D=D,df=lo,model=model,location=location ,scale=scale ,dff = dff, hypothesis =hypothesis)
  pro = pro_compelling_BF(t , df=lo , model=model , location=location ,scale=scale,dff = dff, hypothesis=hypothesis)
  
  if (pro>target){
    N = 2
    return(N)
  }
  
  up= 100000
  if (de_an_prior == 1){
  Power_root <- function(df){
    
    t = BF_bound_10(D=D,df=df,model=model,location=location ,scale=scale ,dff=dff, hypothesis =hypothesis)
    pro = pro_compelling_BF(t , df=df , model=model , location=location ,scale=scale,dff=dff, hypothesis=hypothesis)
    return(target- pro)}
  }else {

    Power_root <- function(df){
      
      t = BF_bound_10(D=D,df=df,model=model,location=location ,scale=scale ,dff=dff, hypothesis =hypothesis)
      pro = pro_compelling_BF(t , df=df , model=model_d , location=location_d ,scale=scale_d,dff=dff_d, hypothesis=hypothesis_d)
      return(target- pro)}
    
    
  }
  N = uniroot(Power_root,lower = 2,upper =  up)$root
  return(N)}


# probability table
Table <- function(D,target,model,location,scale,dff, hypothesis,
                  model_d,location_d,scale_d,dff_d, hypothesis_d,de_an_prior,N, mode  ){
  hypothesis_d= hypothesis 
  if (mode == 0){
  df = N-1
 }else {
  
  df = N_finder(D,target,model,location,scale,dff, hypothesis,
                model_d,location_d,scale_d,dff_d=1, hypothesis_d,de_an_prior  )
  }
  
  t = BF_bound_10(D ,df ,model ,location ,scale ,dff ,hypothesis )
  
  if (de_an_prior == 1 ){
  p_BF10_D_H11 = pro_compelling_BF(t,df ,model ,location ,scale ,dff ,hypothesis )
  p_BF10_D_H01 = false_positive_evidence(t,df,model ,location ,scale,dff, hypothesis)
  
  max_BF = 1/BF10_one_sample(0,df=df,model=model,location=location,scale=scale,dff=dff, hypothesis =hypothesis )
  BF_D = BF_bound_01(D = D, df =df,model = model,location =location,scale=scale,dff = dff, hypothesis)
  
  if (any(hypothesis == "!=" & max_BF<D |BF_D == "bound cannot be found"  )) {
    p_BF10_D_H10 = 0
    p_BF10_D_H00 = 0
  }else{
    t = BF_bound_01(D,df,model ,location ,scale ,dff ,hypothesis )
    p_BF10_D_H10 = false_negative_BF(t,df,model ,location ,scale,dff, hypothesis)
    p_BF10_D_H00 = pro_compelling_BF_null(t , df,model ,location ,scale,dff , hypothesis)
  }
  table <- data.frame(
    H11 = p_BF10_D_H11,
    H10 = p_BF10_D_H10,
    H00 = p_BF10_D_H00,
    H01 = p_BF10_D_H01,
    N =  ceiling(df+1),
    df = df )
  colnames(table) <- c(sprintf("p(BF10> %0.f|H1)",D), sprintf("p(BF01> %0.f|H1)",D), sprintf("p(BF01> %0.f|H0)",D), sprintf("p(BF10> %0.f|H0)",D), " Reqiured N", "exact df")
  
  } else {
    p_BF10_D_H11 = pro_compelling_BF(t,df,model_d ,location_d ,scale_d ,dff_d ,hypothesis_d )
    p_BF10_D_H01 = false_positive_evidence(t,df,model_d ,location_d ,scale_d,dff_d, hypothesis_d)
    
    max_BF = 1/BF10_one_sample(0,df=df,model=model,location=location,scale=scale,dff=dff, hypothesis =hypothesis )
    BF_D = BF_bound_01(D = D, df =df,model = model,location =location,scale=scale,dff = dff, hypothesis)
    
    if (any(hypothesis == "!=" & max_BF<D |BF_D == "bound cannot be found"  )) {
      p_BF10_D_H10 = 0
      p_BF10_D_H00 = 0
    }else{
      t = BF_bound_01(D,df,model ,location ,scale ,dff ,hypothesis )
      p_BF10_D_H10 = false_negative_BF(t,df,model_d,location_d ,scale_d,dff_d, hypothesis_d)
      p_BF10_D_H00 = pro_compelling_BF_null(t , df,model_d ,location_d ,scale_d,dff_d , hypothesis_d)
    }
    table <- data.frame(
      H11 = p_BF10_D_H11,
      H10 = p_BF10_D_H10,
      H00 = p_BF10_D_H00,
      H01 = p_BF10_D_H01,
      N =  ceiling(df+1),
      df = df )
    colnames(table) <- c(sprintf("p(BF10> %0.f|H1)",D), sprintf("p(BF01> %0.f|H1)",D), sprintf("p(BF01> %0.f|H0)",D), sprintf("p(BF10> %0.f|H0)",D), " Reqiured N", "exact df")
    
    }
  return(table)
}

# plot for the selected prior 
prior_plot <-function(D =3,target,model = "NA",location =0,scale=.707,dff = 1, hypothesis,model_d,location_d,scale_d,dff_d=1, hypothesis_d,de_an_prior){
  par(mfrow = c(1, 1))
  bound  <- switch(hypothesis,
                   ">" = c(a = 0, b = 5),
                   "<" = c(a = -5, b = 0),
                   "!=" = c(a = -5, b = 5)
  )
  tt= seq(bound[1],bound[2],.01)
  
  
  
    bound  <- switch(hypothesis,
                   ">" = c(a = 0, b = Inf),
                   "<" = c(a = -Inf, b = 0),
                   "!=" = c(a = -Inf, b = Inf)
  )
  normalization  <- switch(model,
                           "Cauchy" = integrate(function(delta) tstude(delta,location,scale,1),lower = bound[1],upper = bound[2])$value,
                           "Normal" = integrate(function(delta)dnorm(delta,location,scale),lower = bound[1],upper = bound[2])$value,
                           "NLP" = integrate(function(delta)dnlp(delta,location,scale),lower = bound[1],upper = bound[2])$value,
                           "t-distribution" = integrate(function(delta)tstude(delta,location,scale,dff),lower = bound[1],upper = bound[2])$value
  )
  
  
  
  
  
  
  prior_DELTA = NA
  prior_DELTA = switch(model,
                       "Cauchy" = tstude(tt,location,scale,1)/normalization,
                       "Normal" = dnorm(tt,location,scale)/normalization,
                       "t-distribution" = tstude(tt,location,scale,dff)/normalization,
                       "NLP" = dnlp(tt,location,scale)/normalization
  )
  
  plot(tt,prior_DELTA,xlab= bquote(bold(delta)),ylab= "density",type = "l",main  = bquote(bold("prior distribution on "~delta~" under the alternative hypothesis")),frame.plot = FALSE)
if (de_an_prior ==0){
  
  normalization_d  <- switch(model_d,
                           "Cauchy" = integrate(function(delta) tstude(delta,location_d,scale_d,1),lower = bound[1],upper = bound[2])$value,
                           "Normal" = integrate(function(delta)dnorm(delta,location_d,scale_d),lower = bound[1],upper = bound[2])$value,
                           "NLP" = integrate(function(delta)dnlp(delta,location_d,scale_d),lower = bound[1],upper = bound[2])$value,
                           "t-distribution" = integrate(function(delta)tstude(delta,location_d,scale_d,dff_d),lower = bound[1],upper = bound[2])$value
  )
  prior_DELTA_D = switch(model_d,
                       "Cauchy" = tstude(tt,location_d,scale_d,1)/normalization_d,
                       "Normal" = dnorm(tt,location_d,scale_d)/normalization_d,
                       "t-distribution" = tstude(tt,location_d,scale_d,dff_d)/normalization_d,
                       "NLP" = dnlp(tt,location_d,scale_d)/normalization_d
  )
  plot(tt,prior_DELTA,xlab= bquote(bold(delta)),ylim = c(0,max(max(prior_DELTA_D),max(prior_DELTA))),ylab= "density",type = "l",main  = bquote(bold("prior distribution on "~delta~" under the alternative hypothesis")),frame.plot = FALSE)

  lines(tt,prior_DELTA_D ,lty = 2)
  legend("topright", 
         legend = c("Analysis prior", "Design prior"), 
         lty = c(1, 2), 
         col = c("black", "black"),
         bty = "n") 
  
}
  
  
  
  
  
  
  }

# plots for showing the relationship between BF and t-values 

bf10_t <-function(D =3,df, target,model = "NA",location =0,scale=.707,dff = 1, hypothesis ){
  
  tt= seq(from = -5,to = 5,.2)
  BF_D = BF_bound_10(D = D, df =df,model = model,location =location,scale=scale,dff = dff, hypothesis)
  BF10 = BF10_one_sample(tt,df,model ,location,scale,dff,hypothesis)
  
  
  if (length(BF_D) == 1){
    main =  bquote(bold("BF"[10]~"="~.(D) ~"when t = "~.(format(BF_D, digits = 4))))
    #sprintf("BF10 = %.0f when t = %.3f ",D,BF_D)
  } else {
    main =  bquote(bold("BF"[10]~"="~.(D) ~"when t = "~.(format(BF_D[1], digits = 4))~"or"~.(format(BF_D[2], digits = 4))))
    #sprintf("BF10 = %.0f when t = %.3f or %.3f ",D,BF_D[1],BF_D[2])
  }
  par(mfrow = c(1, 2))
  plot(tt,log10(BF10),xlab= "t-value",type="l", ylab = expression("logarithm of BF"[10]),main =   main,frame.plot = FALSE,xaxt = "n")
  abline(v = BF_D)
  axis(1, c(-5,5))
  if (length(BF_D) != 0 ){
  axis(1, round(BF_D,2))}
  
  max_BF = 1/BF10_one_sample(0,df=df,model=model,location=location,scale=scale,dff=dff, hypothesis ="!=" )
  BF_D = BF_bound_01(D = D, df =df,model = model,location =location,scale=scale,dff = dff, hypothesis)

  
  
  plot(tt,log10(1/BF10),xlab= "t-value",type="l",main = "",frame.plot = FALSE,ylab = bquote("logarithm of BF"[0][1]),xaxt = "n")
  axis(1, c(-5,5))
  if (any(hypothesis == "!=" & max_BF<D |BF_D == "bound cannot be found" ) ) {
    main = bquote(bold("It is impossible to have BF"[0][1]~"="~.(D)))
    title(main = main)
    #sprintf("It is impossible to have BF01 = %.3f ",D)
  } else      {
    abline(v = BF_D)
    axis(1, round(BF_D,2))
    if (length(BF_D) == 1){
      main =  bquote(bold("BF"[0][1]~"="~.(D) ~"when t = "~.(format(BF_D, digits = 4))))
      title(main = main)
    } else {
      main =  bquote(bold("BF"[0][1]~"="~.(D) ~"when t = "~.(format(BF_D[1], digits = 4))~"or"~.(format(BF_D[2], digits = 4))))
      title(main = main)
    }}


  

  
}
