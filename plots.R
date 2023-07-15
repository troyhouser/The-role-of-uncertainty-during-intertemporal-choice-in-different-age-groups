m = 0.2
m2 = -0.2
std = 10
N=100
dv1=dv2=c(0)


for(i in 2:N){
  evidence = rnorm(1,m,std)
  dv1 = c(dv1,dv1[i-1]+evidence)
  evidence = rnorm(1,m2,std)
  dv2 = c(dv2,dv2[i-1]+evidence)
}
plot(dv1,type="l",ylim=c(-100,110),lwd=5,col="blue")
lines(dv2,lwd=5,col="orange")
abline(h=110,col="red",lwd=4,lty=2)
abline(h=-100,col="red",lwd=4,lty=2)
abline(a=0,b=1.2,col=adjustcolor("orange",alpha=0.5),lwd=3)
abline(a=0,b=-.75,col=adjustcolor("blue",alpha=0.5),lwd=3)


beta = seq(1,5,length.out=100)
kfunc = function(beta,reward){return(1/(beta*reward))}
plot(beta,kfunc(beta,1),type="l",ylim=c(0,1),col="purple",lwd=3,ylab="temporal discount factor (k)",
     xlab = expression(paste("reward sensitivity (",beta,")")))
lines(beta,kfunc(beta,2),col="darkgreen",lwd=3)
lines(beta,kfunc(beta,3),col="red",lwd=3)
legend(3,.9,legend=c("small reward","medium reward","large reward"),
       col=c("purple","darkgreen","red"),lty=1,cex=0.75)


n=rep(1:100,2)
a=0
b = 1
sigma2 = n^1.3
eps = rnorm(n,mean=0,sd=sqrt(sigma2))
y=a+b*50 + eps
plot(n,y,ylab=c("Reward"),xlab=c("Delay"))
abline(a=50,b=0.4,lwd=3,col="red",lty=2)
abline(a=50,b=-0.5,lwd=3,col="red",lty=2)
abline(v=40,lwd=2,lty=3,col="darkgreen")
