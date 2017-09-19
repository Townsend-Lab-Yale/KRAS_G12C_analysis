
# Calculates lambda, the rate of substitution, using the golden section method

f <- function(lam,n.1,n.0){
  (1-exp(-lam))^n.1 * exp(-lam)^n.0
}

golden_section <- function(f,a=0,b=(n.1/n.0)*2,n.1,n.0,tol=1e-30){
  
  phi <- (sqrt(5)+1)/2
  c <- b-((b-a)/phi)
  d <- a+((b-a)/phi)
  
  while(abs(c-d)>tol){
    if(f(lam=c,n.1=n.1,n.0=n.0) > f(lam=d,n.1=n.1,n.0=n.0)){
      b <- d
    }else{
      a <- c
    }
    c <- b-((b-a)/phi)
    d <- a+((b-a)/phi)
  }
  return((b+a)/2)
}


lambda.calc <- function(n.one,n.zero){
  
  return(golden_section(f,a=0,b=(n.one/n.zero)*2,n.1=n.one,n.0=n.zero))
}