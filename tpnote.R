
# Question 1
# ==========

# Chargement des données
velocite<-read.csv(file = "whales.csv",header =FALSE,dec=".")$V1


# Représentation de l'histogramme
# Ajuster la largeur d'aggrégation peut changer fortement
# l'apparence de la distribution
hist(velocite)
hist(velocite, breaks=18)

nlogvrais=function(theta,x=velocite){
a<-theta[1]
s<-theta[2]
sum(a*log(s)+lgamma(a) -(a-1)*log(x)+x/s)
}

freq<-function(a,s,x){
  1/(s**a*gamma(a))*x**(a-1)*exp(-x/s)
  
  
}

mle<-optim(c(1,2),nlogvrais)
hist(velocite,breaks=18,freq = FALSE)
curve(freq(mle$par[1],mle$par[2],x),from =0.01,to=5,add=T,col="red")


qqplot(qgamma(ppoints(500),shape=mle$par[1],scale=mle$par[2]),y=velocite)
qqline(y=velocite, distribution=function(p) qgamma(p,shape=mle$par[1],scale=mle$par[2]))


rmle<-function(){
  rvelocite<-sample(x=velocite,size = length(velocite),replace =T)

  mle<-optim(c(1,2),function(theta) nlogvrais(theta,rvelocite))
  
  return (mle$par)
  
}

bootstrap<-replicate(1000,rmle())


hist(bootstrap[1,], main='Distribution bootstrap (shape)', xlab='shape (MLE)')
hist(bootstrap[2,], main='Distribution bootstrap (scale)', xlab='scale (MLE)')


et1<-sd(bootstrap[1,])
et1
et2<<-sd(bootstrap[2,])
et2

#borne inférieur
(mle$par[1]-qnorm(0.975)*et1)

#borne supérieur
(mle$par[1]+qnorm(0.975)*et1)

#borne inférieur
(mle$par[2]-qnorm(0.975)*et2)
#borne supérieur
(mle$par[2]+qnorm(0.975)*et2)



library(pracma)
Ifisher=function(x=velocite){
  
  nll<-function(theta,x){
  
     a<-theta[1]
     s<-theta[2]
     (a*log(s)+lgamma(a) -(a-1)*log(x)+x/s)
  
  }

c(shape=mean(sapply( x,function (x ) hessian(function(theta) nll(theta,x),mle$par)[1,1])),scale=mean(sapply( x,function (x ) hessian(function(theta) nll(theta,x),mle$par)[2,2])))
    
}
Ifisher()
sqrt(1/Ifisher())

