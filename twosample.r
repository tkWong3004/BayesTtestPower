
BF10_two_sample <-function(t,n1,r,model ,location,scale,dff , hypothesis ){
  n2 = n1*r
  df = n1+n2-2
  constant = sqrt((n1*n2)/(n1+n2))
  
  bound  <- switch(hypothesis,
                   ">" = c(a = 0, b = Inf),
                   "<" = c(a = -Inf, b = 0),
                   "!=" = c(a = -Inf, b = Inf)
  )
  x <- numeric(length(t))
  
  for(i in 1:length(t)){
    
    
    normalization  <- switch(model,
                             "Cauchy" = integrate(function(delta) dcauchy(delta,location,scale),lower = bound[1],upper = bound[2])$value,
                             "Normal" = integrate(function(delta)dnorm(delta,location,scale),lower = bound[1],upper = bound[2])$value,
                             "NLP" = integrate(function(delta)dnlp(delta,location,scale),lower = bound[1],upper = bound[2])$value,
                             "t-distribution" = integrate(function(delta)tstude(delta,location,scale,dff),lower = bound[1],upper = bound[2])$value
    )
    int  <- switch(model,
                   "Cauchy" = function(delta){
                     dnct(t[i],df,ncp = delta *constant)* dcauchy(delta,location,scale)/normalization},
                   
                   "Normal" = function(delta){
                     dnct(t[i],df,ncp = delta *constant)*dnorm(delta,location,scale)/normalization},
                   
                   "NLP" = int <-function(delta){
                     dnct(t[i],df,ncp = delta *constant)*dnlp(delta,location,scale)/normalization},
                   
                   "t-distribution" = int <-function(delta){
                     dnct(t[i],df,ncp = delta *constant)*tstude(delta,location, scale,dff)/normalization}
    )
    if (scale >.3){
      error = .Machine$double.eps^0.25
    }else {
      error = 1e-8
    }
    if (model == "NLP" & scale <.3 ){
      error = 1e-14
    }
    x[i]= integrate(int,lower = bound[1],upper = bound[2], rel.tol=error,stop.on.error=FALSE)$value
    
  }
  x = x/dnct(t,df,ncp = 0) 
  return(x)
}

BF_bound_10_two <-function(D, n1,r,model ,location ,scale,dff , hypothesis){
  y <- numeric(0)
  Bound_finding <-function(t){
    BF10_two_sample(t,n1,r,model=model,location=location,scale=scale,dff=dff, hypothesis =hypothesis )- D
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
  BF = BF10_two_sample(x,n1,r,model,location,scale,dff,hypothesis )
  x = x[round(BF,1)==round(D,1)]
  if (length(x) ==0){
    x = "no bound is found"
    return(x)
  } 
  return(x)
}

BF_bound_01_two <-function(D , n1,r,model ,location ,scale,dff , hypothesis){
  
  Bound_finding <-function(t){
    1/BF10_two_sample(t, n1,r,model=model,location=location,scale=scale,dff=dff, hypothesis =hypothesis )-D
  }
  
  x = uniroot.all(Bound_finding, lower = -20,upper = 20)
  if (length(x) == 0 ){
    x = "bound cannot be found"
    return(x)
  }
  BF = 1/BF10_two_sample(x , n1,r,model=model,location=location,scale=scale,dff=dff, hypothesis =hypothesis )
  x = x[round(BF,1)== round(D,1)]
  
  return(x)
}

pro_compelling_BF_two<-function(t,n1,r,model ,location ,scale,dff , hypothesis ){
  n2 = n1*r
  df = n1+n2-2
  constant =  sqrt((n1*n2)/(n1+n2))
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
                           "Cauchy" = integrate(function(delta) tstude(delta,location,scale,df=1),lower = bound[1],upper = bound[2])$value,
                           "Normal" = integrate(function(delta)dnorm(delta,location,scale),lower = bound[1],upper = bound[2])$value,
                           "NLP" = integrate(function(delta)dnlp(delta,location,scale),lower = bound[1],upper = bound[2])$value,
                           "t-distribution" = integrate(function(delta)tstude(delta,location,scale,dff),lower = bound[1],upper = bound[2])$value
  )
  if (length(t) > 1){
    
    int  <- switch(model,
                   "Cauchy" = function(delta){
                     pro1 = pnct(t[t>0],df,ncp = delta *constant,lower = F)
                     pro2 = pnct(t[t<0],df,ncp = delta *constant,lower = T)
                     (pro1+pro2)* tstude(delta,location,scale,df=1)/normalization },
                   
                   "Normal" = function(delta){
                     pro1 = pnct(t[t>0],df,ncp = delta *constant,lower = F)
                     pro2 = pnct(t[t<0],df,ncp = delta *constant,lower = T)
                     
                     (pro1+pro2)* dnorm(delta,location,scale)/normalization   },
                   
                   "NLP" = function(delta){
                     pro1 = pnct(t[t>0],df,ncp = delta *constant,lower = F)
                     pro2 = pnct(t[t<0],df,ncp = delta *constant,lower = T)
                     
                     (pro1+pro2)* dnlp(delta,location,scale)/normalization   } ,
                   
                   "t-distribution" = function(delta){
                     pro1 = pnct(t[t>0],df,ncp = delta *constant,lower = F)
                     pro2 = pnct(t[t<0],df,ncp = delta *constant,lower = T)
                     
                     (pro1+pro2)* tstude(delta,location,scale,dff)/normalization   }
    )
    
    
    
  }  else if (t >= 0) { 
    int  <- switch(model,
                   "Cauchy" = function(delta){
                     pro1 = pnct(t,df,ncp = delta *constant,lower = F)
                     (pro1)*dcauchy(delta,location,scale)/normalization  },
                   
                   "Normal" = function(delta){
                     pro1 = pnct(t,df,ncp = delta *constant,lower = F)
                     (pro1)* dnorm(delta,location,scale)/normalization   },
                   
                   "NLP" = function(delta){
                     pro1 = pnct(t,df,ncp = delta *constant,lower = F)
                     (pro1)* dnlp(delta,location,scale)/normalization   },
                   
                   "t-distribution" = function(delta){
                     pro1 = pnct(t,df,ncp = delta *constant,lower = F)
                     (pro1)* tstude(delta,location,scale,dff)/normalization}
    )
    
  } else if (t < 0){
    int  <- switch(model,
                   "Cauchy" = function(delta){
                     pro2 = pnct(t,df,ncp = delta *constant,lower = T)
                     (pro2)* dcauchy(delta,location,scale)/normalization},
                   
                   
                   "Normal" = function(delta){
                     pro2 = pnct(t,df,ncp = delta *constant,lower = T)
                     (pro2)* dnorm(delta,location,scale)/normalization   },
                   
                   "NLP" = function(delta){
                     pro2 = pnct(t,df,ncp = delta *constant,lower = T)
                     (pro2)* dnlp(delta,location,scale)/normalization},
                   
                   "t-distribution" = int <-function(delta){
                     pro2 = pnct(t,df,ncp = delta *constant,lower = T)
                     (pro2)* tstude(delta,location,scale,dff)/normalization   }
    )
    
  }
  if (scale >.3){
error = .Machine$double.eps^0.25
  }else {
   error = 1e-8
 }
  
  if (model == "NLP" & scale <.3 ){
    error = 1e-10
  }
  x = integrate(int,lower = bound[1],upper = bound[2], rel.tol = error)$value
  
  return(x) 
  
}

pro_compelling_BF_null_two <- function(t , n1,r,model,location ,scale,dff , hypothesis){
  if (any(t =="no bound is found" | length(t)==0)){
    t=0
    return(t)
  }
  
  df = n1+n1*r-2

  if (length(t)>1){
    pro = pnct(max(t),df) - pnct(min(t),df)
  }else if (hypothesis ==">"){
    pro = pnct(t,df=df,lower = T)
  }else if (hypothesis =="<"){
    pro = pnct(t,df=df,lower = F)
  }
  return(pro)
}

false_positive_evidence_two <- function(t_cv,n1,r,model ,location ,scale,dff, hypothesis){
  if (any(t_cv =="no bound is found" | length(t_cv)==0)){
    t_cv=0
    return(t_cv)
  }
  
  df = n1+n1*r-2
  if (length(t_cv) == 4) {  # Corrected condition check
    t_cv = t_cv[2:3]
  }
  if (length(t_cv)>1){
    pro = pnct(max(t_cv),df=df,lower = F) +pnct(min(t_cv),df=df,lower = T)
  }else if (hypothesis ==">"){
    pro = pnct(t_cv,df=df,lower = F)
  }else if (hypothesis =="<"){
    pro = pnct(t_cv,df=df,lower = T)
  }
  return(pro)
  
}

false_negative_BF_two<-function(t,n1,r,model ,location ,scale,dff , hypothesis ){
  
  n2 = n1*r
  df= n1+n2-2
  constant =  sqrt((n1*n2)/(n1+n2))
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
                           "Cauchy" = integrate(function(delta) dcauchy(delta,location,scale),lower = bound[1],upper = bound[2])$value,
                           "Normal" = integrate(function(delta)dnorm(delta,location,scale),lower = bound[1],upper = bound[2])$value,
                           "NLP" = integrate(function(delta)dnlp(delta,location,scale),lower = bound[1],upper = bound[2])$value,
                           "t-distribution" = integrate(function(delta)tstude(delta,location,scale,dff),lower = bound[1],upper = bound[2])$value
  )
  if (length(t)>1){
    
    int  <- switch(model,
                   "Cauchy" = function(delta){
                     pro1 = pnct(max(t),df,ncp = delta *constant,lower = T)
                     pro2 = pnct(min(t),df,ncp = delta *constant,lower = T)
                     (pro1-pro2)* dcauchy(delta,location,scale)/normalization },
                   
                   "Normal" = function(delta){
                     pro1 = pnct(max(t),df,ncp = delta *constant,lower = T)
                     pro2 = pnct(min(t),df,ncp = delta *constant,lower = T)
                     
                     (pro1-pro2)* dnorm(delta,location,scale)/normalization   },
                   
                   "NLP" = function(delta){
                     pro1 = pnct(max(t),df,ncp = delta *constant,lower = T)
                     pro2 = pnct(min(t),df,ncp = delta *constant,lower = T)
                     
                     (pro1-pro2)* dnlp(delta,location,scale)/normalization   } ,
                   
                   "t-distribution" = function(delta){
                     pro1 = pnct(max(t),df,ncp = delta *constant,lower = T)
                     pro2 = pnct(min(t),df,ncp = delta *constant,lower = T)
                     
                     (pro1-pro2)* tstude(delta,location,scale,dff)/normalization   }
    )
    
  } else if (t>0) { 
    int  <- switch(model,
                   "Cauchy" = function(delta){
                     pro1 = pnct(t,df,ncp = delta *constant,lower = T)
                     (pro1)*dcauchy(delta,location,scale)/normalization  },
                   
                   "Normal" = function(delta){
                     pro1 = pnct(t,df,ncp = delta *constant,lower = T)
                     (pro1)* dnorm(delta,location,scale)/normalization   },
                   
                   "NLP" = function(delta){
                     pro1 = pnct(t,df,ncp = delta *constant,lower = T)
                     (pro1)* dnlp(delta,location,scale)/normalization   },
                   
                   "t-distribution" = function(delta){
                     pro1 = pnct(t,df,ncp = delta *constant,lower = T)
                     (pro1)* tstude(delta,location,scale,dff)/normalization}
    )
    
  } else if (t<0){
    int  <- switch(model,
                   "Cauchy" = function(delta){
                     pro2 = pnct(t,df,ncp = delta *constant,lower = F)
                     (pro2)* dcauchy(delta,location,scale)/normalization},
                   
                   
                   "Normal" = function(delta){
                     pro2 = pnct(t,df,ncp = delta *constant,lower = F)
                     (pro2)* dnorm(delta,location,scale)/normalization   },
                   
                   "NLP" = function(delta){
                     pro2 = pnct(t,df,ncp = delta *constant,lower = F)
                     (pro2)* dnlp(delta,location,scale)/normalization},
                   
                   "t-distribution" = int <-function(delta){
                     pro2 = pnct(t,df,ncp = delta *constant,lower = F)
                     (pro2)* tstude(delta,location,scale,dff)/normalization   }
    )
    
  }
  x = integrate(int,lower = bound[1],upper = bound[2])$value
  
  
  
  if (length(t) == 1 ){
    if (t>0 & hypothesis == "<"){
      x = 1-x
    }
    if (t<0 & hypothesis == ">") {
      x = 1-x
    }}
  
  return(x) 
} 

N1_finder_two<-function(D,r,target,model,location,scale,dff, hypothesis,
                        model_d,location_d,scale_d,dff_d, hypothesis_d,de_an_prior){

  lo = 2
  t = BF_bound_10_two(D=D,n1=lo,r=r,model=model,location=location ,scale=scale ,dff = dff, hypothesis =hypothesis)
  pro = pro_compelling_BF_two(t , n1=lo,r=r , model=model , location=location ,scale=scale,dff = dff, hypothesis=hypothesis)
  
  if (pro>target){
    N = 2
    return(N)
  }
  
  up = 100000
  

  if (de_an_prior == 1){
  Power_root_two <- function(n1){
    
    t = BF_bound_10_two(D=D,n1=n1,r=r,model=model,location=location ,scale=scale ,dff=dff, hypothesis =hypothesis)
    pro = pro_compelling_BF_two(t ,n1=n1,r=r , model=model , location=location ,scale=scale,dff=dff, hypothesis=hypothesis)
    return(target- pro)}
  }else {
    Power_root_two <- function(n1){
      t = BF_bound_10_two(D=D,n1=n1,r=r,model=model,location=location ,scale=scale ,dff=dff, hypothesis =hypothesis)
      pro = pro_compelling_BF_two(t ,n1=n1,r=r , model=model_d , location=location_d ,scale=scale_d,dff=dff_d, hypothesis=hypothesis_d)
      return(target- pro)}
  }
  if (target == .75) {
    up = 100000
  }
  N = uniroot(Power_root_two,lower = 2,upper =  up)$root
  return(N)
}

Table_two <- function(D,r,target,model,location,scale,dff, hypothesis, model_d,location_d,scale_d,dff_d, hypothesis_d,de_an_prior,N1,n2,mode2){
  hypothesis = hypothesis_d
  if (mode2 == 0){
    r = n2/N1
    n1 =N1
  }else {
    
    n1 = N1_finder_two(D,r,target,model,location ,scale  ,dff ,hypothesis, model_d,location_d,scale_d,dff_d, hypothesis_d,de_an_prior )
  }

  t = BF_bound_10_two(D,n1,r,model ,location ,scale ,dff ,hypothesis )
  
  if (de_an_prior == 1 ){
  p_BF10_D_H11 = pro_compelling_BF_two(t,n1,r ,model ,location ,scale ,dff ,hypothesis )
  p_BF10_D_H01 = false_positive_evidence_two(t,n1,r ,model ,location ,scale ,dff ,hypothesis )

  max_BF = 1/BF10_two_sample(0, n1,r ,model ,location ,scale ,dff ,hypothesis )
  BF_D = BF_bound_01_two(D , n1,r,model ,location ,scale,dff , hypothesis)
  
  if (any(hypothesis == "!=" & max_BF<D |BF_D == "bound cannot be found"  )) {
    p_BF10_D_H10 = 0
    p_BF10_D_H00 = 0
  }else{
    p_BF10_D_H10 = false_negative_BF_two(BF_D,n1,r ,model ,location ,scale ,dff ,hypothesis )
    p_BF10_D_H00 = pro_compelling_BF_null_two(BF_D,n1,r,model ,location ,scale ,dff ,hypothesis)
  }
  table <- data.frame(
    H11 = p_BF10_D_H11,
    H10 = p_BF10_D_H10,
    H00 = p_BF10_D_H00,
    H01 = p_BF10_D_H01,
    N1 =  ceiling(n1),
    N2 = ceiling(n1*r) 
  )
  colnames(table) <- c(sprintf("p(BF10> %0.f|H1)",D), sprintf("p(BF01> %0.f|H1)",D), sprintf("p(BF01> %0.f|H0)",D), sprintf("p(BF10> %0.f|H0)",D), " Reqiured N in group 1", "Reqiured N in group 2")
  }else{
    
    p_BF10_D_H11 = pro_compelling_BF_two(t,n1,r ,model_d ,location_d ,scale_d ,dff_d ,hypothesis_d )
    p_BF10_D_H01 = false_positive_evidence_two(t,n1,r ,model_d ,location_d ,scale_d ,dff_d ,hypothesis_d )

    max_BF = 1/BF10_two_sample(0, n1,r ,model ,location ,scale ,dff ,hypothesis )
    BF_D = BF_bound_01_two(D , n1,r,model ,location ,scale,dff , hypothesis)
    
    if (any(hypothesis == "!=" & max_BF<D |BF_D == "bound cannot be found"  )) {
      p_BF10_D_H10 = 0
      p_BF10_D_H00 = 0
    }else{
      p_BF10_D_H10 = false_negative_BF_two(BF_D,n1,r ,model_d ,location_d ,scale_d ,dff_d ,hypothesis_d )
      p_BF10_D_H00 = pro_compelling_BF_null_two(BF_D,n1,r,model_d ,location_d ,scale_d ,dff_d ,hypothesis_d)
    }
    table <- data.frame(
      H11 = p_BF10_D_H11,
      H10 = p_BF10_D_H10,
      H00 = p_BF10_D_H00,
      H01 = p_BF10_D_H01,
      N1 =  ceiling(n1),
      N2 = ceiling(n1*r) 
    )
    colnames(table) <- c(sprintf("p(BF10> %0.f|H1)",D), sprintf("p(BF01> %0.f|H1)",D), sprintf("p(BF01> %0.f|H0)",D), sprintf("p(BF10> %0.f|H0)",D), " Reqiured N in group 1", "Reqiured N in group 2")
    
  }
  return(table)
}

bf10_two <-function(D ,n1,r, target,model,location ,scale,dff = 1, hypothesis ){
  
  tt= seq(from = -5,to = 5,.1)
  #df = N_finder(D,target ,model,location,scale,dff,hypothesis)
  BF_D = BF_bound_10_two(D,n1,r,model ,location ,scale ,dff ,hypothesis )
  BF10 = BF10_two_sample(tt,n1,r,model ,location,scale,dff,hypothesis)
  
  
  if (length(BF_D) == 1){
    main =  bquote(bold("BF"[10]~"="~.(D) ~"when t = "~.(format(BF_D, digits = 4))))
    #sprintf("BF10 = %.0f when t = %.3f ",D,BF_D)
  } else {
    main =  bquote(bold("BF"[10]~"="~.(D) ~"when t = "~.(format(BF_D[1], digits = 4))~"or"~.(format(BF_D[2], digits = 4))))
    #sprintf("BF10 = %.0f when t = %.3f or %.3f ",D,BF_D[1],BF_D[2])
  }
  par(mfrow = c(1, 2))
  plot(tt,log(BF10),xlab= "t-value",type="l", ylab = expression("logarithm of BF"[10]),main =   main,frame.plot = FALSE,xaxt = "n")
    abline(v = BF_D)
  axis(1, c(-5,5))
  if (length(BF_D) != 0 ){
    axis(1, round(BF_D,2))}
  
  max_BF = 1/BF10_two_sample(0, n1,r ,model ,location ,scale ,dff ,hypothesis )
  BF_D = BF_bound_01_two(D , n1,r,model ,location ,scale,dff , hypothesis)
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


