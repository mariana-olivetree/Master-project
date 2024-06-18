LB <- function(data, outcome, transform = F, alpha = 0.05){
  
  # data: a matrix records the count datA
  # outcome: a categorical vector of a specific kind of clinical condition
  # transform: square-root transform of the compositional matrix
  # alpha: customized threshold while calculating q values
  
data_filtered = which(colSums(data)>2*dim(data)[1])

final_data = data[, data_filtered]

sample_total = rowSums(final_data)

taxa = dim(final_data)[2]

beta = matrix(data = NA,taxa,2)

for (i in 1:taxa) {
  
  x.prop = final_data[,i]/sample_total
  
  if (transform == T){
    
    x.prop=sqrt(x.prop)
  }
  
  bereg=gamlss(x.prop ~ outcome, family=BEZI(sigma.link="identity"),trace=FALSE,control = gamlss.control(n.cyc = 100))
  
  out=summary(bereg)
  
  beta[i,]=out[2,c(1,4)]}  ## get all the coefficients and P values in beta regression

pvalues=beta[,2]

qvalues=calc_qvalues(pvalues)

sig=which(qvalues<alpha)

sigFeature=colnames(X)[sig]

  list(sigFeature=sigFeature,useFeature=P,qvalues=qvalues, pvalues=pvalues) }




