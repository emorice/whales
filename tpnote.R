#
# Question 1
# ==========

# Chargement des données
velocite<-read.csv(file = "whales.csv",header =FALSE,dec=".")$V1


# Représentation de l'histogramme
# Ajuster la largeur d'aggrégation peut changer fortement
# l'apparence de la distribution
hist(velocite)
hist(velocite, breaks=18)

# On définit la fonction à optimiser :
# L'opposé de la log-vraisemblance
# Le parametre est recodé dans un vecteur c(shape, scale)
nlogvrais=function(theta,x=velocite){
	a<-theta[1]
	s<-theta[2]
	# Stabilisation: il est assez peu pratique d'imposer aux méthodes
	# d'optimisation de ne donner que des valeurs positives aux paramètres,
	# on "clippe" donc les valeurs interdites en entrée. Bien sûr, si une
	# méthode d'optimisation ressort des valeurs négatives, il faudra
	# comprendre que l'algorithme a divergé.
	# Les autres solutions consisteraient à changer le paramétrage pour
	# travailler avec le log des paramètres, ou à choisir des algorithmes
	# d'optimisation qui acceptent des contraintes en entrée.
	if(a <= 0)
		a = 1e-7
	if(s <= 0)
		s = 1e-7
	sum(a*log(s)+lgamma(a) -(a-1)*log(x)+x/s)
}

# On definit au passage la densité associée
freq<-function(a,s,x){
  1/(s**a*gamma(a))*x**(a-1)*exp(-x/s)
}

# Remarque: on ici codé littéralement les [log-]densités, mais on aurait aussi
# pu utiliser de façon équivalente la fonction dgamma de R

# Estimation : on cherche numériquement le minimum de nlogvrais
# Les parametres initiaux sont choisi arbitrairement pour donner une
# approximation initiale de la bonne forme
# La méthode d'optimisation par défaut utilisée est une méthode boîte noire de
# Nelder-Mead qui ne recourt pas aux dérivées de la fonction, elle est assez peu
# efficace, mais robuste et suffit largement pour la taille des données
# considérées ici.

mle<-optim(c(1,2),nlogvrais)

cat("Valeurs obtenues par MLE\n")
c(shape=mle$par[1], scale=mle$par[2])
cat('\n')

# Variante plus compacte avec les built-in de R:
#mle_R<-optim(c(1,2),function(t) -sum(log(dgamma(velocite, shape=t[1], scale=t[2]))))
#mle_R

# Diagnostiques
## Histogramme + fit
hist(velocite,breaks=18,freq = FALSE)
curve(freq(mle$par[1],mle$par[2],x),from =0.01,to=5,add=T,col="red")
#curve(dgamma(shape=mle$par[1],scale=mle$par[2],x=x),from =0.01,to=5,add=T,col="blue")

## QQ-plot
# Ici, recoder les quantiles de la loi gamma est peu pratique, on a donc utilisé
# qgamma
qqplot(qgamma(ppoints(500),shape=mle$par[1],scale=mle$par[2]),y=velocite)
qqline(y=velocite, distribution=function(p) qgamma(p,shape=mle$par[1],scale=mle$par[2]))

# Question 2
# ==========

# On définit une étape du bootstrap : rééchantillonage et MLE.
rmle<-function(){
  rvelocite<-sample(x=velocite,size = length(velocite),replace =T)

  mle<-optim(c(1,2),function(theta) nlogvrais(theta,rvelocite))
  
  return (mle$par)
  
}

# Puis on l'exécute un bon nombre de fois
bootstrap<-replicate(1000,rmle())

# Histogrammes des marginales
hist(bootstrap[1,], main='Distribution bootstrap (shape)', xlab='shape (MLE)')
hist(bootstrap[2,], main='Distribution bootstrap (scale)', xlab='scale (MLE)')

# Evaluation des écart-types bootstrap
et1<-sd(bootstrap[1,])
et2<<-sd(bootstrap[2,])
cat('Écarts-types des deux distributions des estimateurs par bootstrap\n')
c(shape=et1, scale=et2)
cat('\n')


# Question 3
# ==========
# On construit des intervalles de confiance gaussiens symétriques à 5% de risque
# à partir des écarts-types estimés :
cat('Intervalles de confiance bootstrap à 95% (via estimation gaussienne)\n')
c(shape=c(
	#borne inférieure
	lower=(mle$par[1]-qnorm(0.975)*et1),
	#borne supérieure
	upper=(mle$par[1]+qnorm(0.975)*et1)
	),
  scale=c(
	#borne inférieure
	lower=(mle$par[2]-qnorm(0.975)*et2),
	#borne supérieure
	upper=(mle$par[2]+qnorm(0.975)*et2)
	)
  )
cat('\n')

# Intervalles empiriques
cat('Intervalles de confiance bootstrap à 95% (complètement non-paramétriques)\n')
cat('Shape\n')
quantile(bootstrap[1,], c(0.025, 0.975)) - mean(bootstrap[1,]) + mle$par[1]
cat('Scale\n')
quantile(bootstrap[2,], c(0.025, 0.975)) - mean(bootstrap[2,]) + mle$par[2]
cat('\n')

# Question 4
# ==========

# On estime l'information de Fischer en moyennant la hessienne de la
# vraisemblance au point du MLE sur le jeu de données.
# Ici, on calcule numériquement la hessienne

library(pracma)
Ifisher=function(x=velocite){
  
	# negative log-likelihood en un seul point
	nll<-function(theta,x){
		a<-theta[1]
		s<-theta[2]
		(a*log(s)+lgamma(a) -(a-1)*log(x)+x/s)
	}

	# Hessienne
	Hflat = rowMeans(
		sapply(
			x,
			function(x) hessian(function(theta) nll(theta,x), mle$par)
			)
		)
	matrix(Hflat, nrow=2)
}
cat('Information de Fisher estimée\n')
If = Ifisher()
If
cat('\n')

cat('Covariance estimée de l\'EMV\n')
C = solve(If)/length(velocite)
C
cat('\n')

cat('Écarts-type marginaux estimés pour les EMV de façon classique\n')
c(shape=sqrt(C[1,1]), scale=sqrt(C[2,2]))
cat('\n')

# Il est bien sûr possible d'obtenir ces résultats plus directement avec les
# paquets R adaptés:
#library(stats4)
#fit = mle(function(shape, scale) -sum(log(dgamma(velocite, shape=shape, scale=scale))), list(shape=1, scale=1))
#summary(fit)

# Autres
# ======
# Normalité des bootsraps
library(moments)
qqnorm(bootstrap[1,])
qqline(bootstrap[1,])
#skewness(bootstrap[1,])
qqnorm(bootstrap[2,])
qqline(bootstrap[2,])
#skewness(bootstrap[2,])
