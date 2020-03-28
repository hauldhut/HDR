#F = w*t(FDr)+ (1 − w)*FD
#FDr= WDr(WDr + nhetaDr*IDr)*t(A)
#FD= WD(WD + nhetaD*ID)*A

#setwd('/Users/admin/Data/DR/Fdatasets/')
#setwd('/Users/admin/Data/DR/Cdatasets/')
#setwd('/Users/admin/Data/DR/DNdatasets/')

RLS <- function(A, WDr, WD, nhetaDr = 1, nhetaD = 1, w=0.5) {
  IDr<-diag(nrow((WDr)))
  FDr<-WDr%*%(WDr + nhetaDr*IDr)%*%A
  ID<-diag(nrow(WD))
  FD<-WD%*%(WD + nhetaD*ID)%*%t(A)
  F=(1-w)*FDr + w*t(FD)
  return(F)
}