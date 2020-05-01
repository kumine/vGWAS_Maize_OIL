#####---useful function------------------------------------------------
#Levene's test for homogeneity of variance across groups
levene.test=function(y,group){
  group <- as.factor(group)
  valid <- complete.cases(y, group)
  y=y[valid]
  group=group[valid]
  meds <- tapply(y, group, median)
  resp <- abs(y - meds[group])
  res<- anova(lm(resp ~ group))
  res[1,5]
}

##======
#Estimate the inflation factor for a distribution of P-values.
#The major use of this procedure is the Genomic Control
# modified from estlambda {GenABEL}
vGWAS.gc <- function(object, proportion = 1) 
{
  if (proportion > 1 || proportion <= 0) 
    stop('proportion argument should be greater then zero and less than or equal to one.')
  ntp <- round(proportion * length(object$p.value))
  if (ntp <= 1) 
    stop('too few valid measurments.')
  if (ntp < 10) 
    warning(paste('number of points is fairly small:', ntp))
  if (min(object$p.value) < 0) 
    stop('data argument has values <0')
  if (max(object$p.value) <= 1) {
    data <- data0 <- qchisq(object$p.value, 1, lower.tail = FALSE)
  }
  data <- sort(data)
  ppoi <- ppoints(data)
  ppoi <- sort(qchisq(1 - ppoi, 1))
  data <- data[1:ntp]
  ppoi <- ppoi[1:ntp]
  s <- summary(lm(data ~ 0 + ppoi))$coeff
  out <- object
  
  out$p.value.gc <- pchisq(data0/s[1, 1], 1, lower.tail = FALSE)
  return(out)
}
