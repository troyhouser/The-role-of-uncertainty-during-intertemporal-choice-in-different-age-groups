itc = read.csv("~/Dropbox (University of Oregon)/ITC_dev/preprocessed_choicedata_devITCstudy.csv")
itc$LL = itc$Choice_SS0LL1
S = unique(itc$ppnr)
itc$X1 = itc$SSamount
itc$X2 = itc$LLamount
itc$T1 = itc$SStime/365
itc$T2 = itc$LLtime/365
itc$s2_u=0

safelog = function(x){
  x[x==0]=1e-100
  y=log(x)
  return(y)
}
R1 = function(pars,dat){
  beta = pars[1]; alpha = pars[2]; lapse = 0
  r = cbind(dat$X2,dat$X1)
  r[r==0]=1e-6
  t = cbind(dat$T2,dat$T1)
  s2_e = dat$s2_u/(abs(r)*beta)
  k = s2_e/dat$s2_u
  V = r/(1+(k*t))
  P = (1-lapse)*pnorm(alpha*(V[,1]-V[,2]))+lapse*0.5
  lik = safelog(t(P))*dat$LL+safelog(1-t(P))*(1-dat$LL)
  fit = sum(-lik,na.rm=T)
  return(fit)
}
R2 = function(pars,dat){
  beta = pars[1]; lapse = 0
  r = cbind(dat$X2,dat$X1);r[r==0]=1e-6
  t = cbind(dat$T2,dat$T1)
  s2_e = dat$s2_u/(abs(r)*beta)
  k = s2_e/dat$s2_u
  D = 1 / (1+(k*t))
  V = r*D
  alpha = 1/sqrt(rowSums((D^2)*s2_e))
  P = (1-lapse)*pnorm(alpha*(V[,1]-V[,2]))+lapse*0.5
  lik = safelog(t(P))*dat$LL+safelog(1-t(P))*(1-dat$LL)
  fit = sum(-lik,na.rm=T)
  return(fit)
}
R3 = function(pars,dat){
  beta = pars[1]; alpha = pars[2]; lapse = 0; k0 = pars[3]
  r = cbind(dat$X2,dat$X1);r[r==0]=1e-6
  t = cbind(dat$T2,dat$T1)
  s2_e = dat$s2_u/(abs(r)*beta)
  k = k0+s2_e/dat$s2_u
  V = r/(1+(k*t))
  P = (1-lapse)*pnorm(alpha*(V[,1]-V[,2]))+lapse*0.5
  lik = safelog(t(P))*dat$LL+safelog(1-t(P))*(1-dat$LL)
  fit = sum(-lik,na.rm=T)
  return(fit)
}
QH = function(pars,dat){
  d = pars[1]; b = pars[2]; alpha = pars[3]
  lapse = 0
  r = cbind(dat$X2,dat$X1);r[r==0]=1e-6
  t = cbind(dat$T2,dat$T1)
  V = r*b*(d^t)
  P = (1-lapse)*pnorm(alpha*(V[,1]-V[,2]))+lapse*0.5
  lik = safelog(t(P))*dat$LL+safelog(1-t(P))*(1-dat$LL)
  fit = sum(-lik,na.rm=T)
  return(fit)
}
H0 = function(pars,dat){
  alpha = pars[1]; k = pars[2]
  lapse = 0
  r = cbind(dat$X2,dat$X1);r[r==0]=1e-6
  t = cbind(dat$T2,dat$T1)
  V = r/(1+(k*t))
  P = (1-lapse)*pnorm(alpha*(V[,1]-V[,2]))+lapse*0.5
  lik = safelog(t(P))*dat$LL+safelog(1-t(P))*(1-dat$LL)
  fit = sum(-lik,na.rm=T)
  return(fit)
}

rtnorm <- function(n, mean, sd, a = -Inf, b = Inf){
  qnorm(runif(n, pnorm(a, mean, sd), pnorm(b, mean, sd)), mean, sd)
}

r1=r2=r3=h0=qh=list()
ntrials = c()
for(s in 1:length(S)){
  ix = itc[itc$ppnr==S[s],]
  ntrials[s] = nrow(ix)
  for(n in 1:nrow(ix)){
    ix$s2_u[n] = max(1e-6,var(c(ix$X1[1:n],ix$X2[1:n]),na.rm=T))
  }
  r1[[s]] = nlminb(start = rtnorm(2,0,1,0,1), objective = R1, lower = c(1e-6,1e-6), upper = c(100,50), dat = ix)
  r2[[s]] = nlminb(start = rtnorm(1,0,1,0,1), objective = R2, lower = c(1e-6), upper = c(100), dat = ix)
  r3[[s]] = nlminb(start = rtnorm(3,0,1,0,1), objective = R3, lower = rep(1e-6,3), upper = c(50,100,50), dat = ix)
  h0[[s]] = nlminb(start = rtnorm(2,0,1,0,1), objective = H0, lower = rep(1e-6,2), upper = c(50,100), dat = ix)
  qh[[s]] = nlminb(start = rtnorm(3,0,1,0,1), objective = QH, lower = rep(1e-6,3), upper = c(1,1,10), dat = ix)
}
bic = matrix(0,length(S),5)
par_est = matrix(0,length(S),2)
ages_cont=ages_fact = ll = c()
for(i in 1:length(S)){
  bic[i,1] = 2*r1[[i]]$objective+log(ntrials[i])*2
  bic[i,2] = 2*r2[[i]]$objective+log(ntrials[i])
  bic[i,3] = 2*r3[[i]]$objective+log(ntrials[i])*3
  bic[i,4] = 2*h0[[i]]$objective+log(ntrials[i])*2
  bic[i,5] = 2*qh[[i]]$objective+log(ntrials[i])*3
  
  par_est[i,1] = r1[[i]]$par[1]
  par_est[i,2] = r1[[i]]$par[2]
  
  ages_cont[i] = itc$age[itc$ppnr==S[i]]
  ages_fact[i] = itc$agegroup[itc$ppnr==S[i]]
  ll[i] = mean(itc$LL[itc$ppnr==S[i]],na.rm=T)
}
colMeans(bic);apply(bic,2,sd)
t.test(bic[,1],bic[,4],paired=T)
lsr::cohensD(bic[,1],bic[,4])
t.test(bic[-21,1],bic[-21,4],paired=T)
lsr::cohensD(bic[-21,1],bic[-21,4])
bic_df = data.frame(BIC = c(bic[,1],bic[,2],bic[,3],bic[,4],bic[,5]),
                    sub = rep(1:nrow(bic),5),
                    model = rep(c("RI1","RI2","RI3","HB","QH"),each=nrow(bic)))

p1a = ggplot(bic_df[bic_df$BIC<100,],aes(x=model,y=BIC,fill=model,col=model))+
  geom_line(aes(group=sub),alpha=0.2,col="gray")+
  geom_bar(stat="summary",col="black")+
  geom_jitter(shape=21,col="black",width=0.2)+
  geom_errorbar(stat="summary",width=0.2,col="black",linewidth=1)+
  scale_y_continuous(limits = c(),expand= c(0,0))+
  scale_fill_manual(values=c("#EBD913","#B8AA0F","#3387B8","#42ADEB","#3B9AD1"))+
  ylab("Bayesian Information Criterion")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),
        axis.title=element_text(size=15))+
  theme(legend.position = "None")
p1a
ggsave(plot=p1a,filename="p1a.png",units="px",width=2400,height=2400,
       path="~/Dropbox (University of Oregon)/ITC_dev/plots/")

df = data.frame(sub = S,
                agegroup = ages_fact,
                age_num = ages_cont,
                beta = par_est[,1],
                alpha = par_est[,2],
                LL = ll)
m_beta=mean(df$beta);std_beta=sd(df$beta)
m_alpha=mean(df$alpha);std_alpha=sd(df$alpha)
crit_beta = m_beta+2*std_beta;crit_alpha = m_alpha+2*std_alpha
df = df[!df$beta>crit_beta,]
df = df[!df$alpha>crit_alpha,]
nrow(df)
df$agegroup=factor(df$agegroup,levels=c("children","adolescents","young adults"))
p1=ggplot(df,aes(x=agegroup,y=beta,fill=agegroup))+
  geom_bar(stat="summary",col="black")+
  geom_errorbar(stat="summary",width=0.2,col="black",linewidth=1)+
  scale_y_continuous(limits = c(),expand= c(0,0))+
  scale_fill_manual(values=c("#07E34F","#FF9800","#8C19FF"))+
  ylab(expression(paste("Reward Sensitivity (",beta,")")))+
  xlab("age group")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),
        axis.title=element_text(size=15))+
  theme(legend.position = "None")

p2=ggplot(df,aes(x=agegroup,y=alpha,fill=agegroup))+
  geom_bar(stat="summary",col="black")+
  geom_errorbar(stat="summary",width=0.2,col="black",linewidth=1)+
  scale_y_continuous(limits = c(),expand= c(0,0))+
  scale_fill_manual(values=c("#07E34F","#FF9800","#8C19FF"))+
  ylab(expression(paste("Temperature (",alpha,")")))+
  xlab("")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),
        axis.title=element_text(size=15))+
  theme(legend.position = "None")
p3=ggplot(df,aes(x=agegroup,y=LL,fill=agegroup))+
  geom_bar(stat="summary",col="black")+
  geom_errorbar(stat="summary",width=0.2,col="black",linewidth=1)+
  scale_y_continuous(limits = c(),expand= c(0,0))+
  scale_fill_manual(values=c("#07E34F","#FF9800","#8C19FF"))+
  ylab("Proportion of LL Choices")+
  xlab("")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),
        axis.title=element_text(size=15))+
  theme(legend.position = "None")
p123 = ggpubr::ggarrange(p3,p1,p2,nrow=1,labels=c("A","B","C"))
p123

ggsave(plot=p123,filename="p123.png",units="px",width=5400,height=2400,
       path="~/Dropbox (University of Oregon)/ITC_dev/plots/")


#### regression predicting age from parameter estimates
sd(df$beta[df$agegroup=="children"])
m1 = lm(age_num~beta+alpha,df)
summary(m1)
cor.test(df$beta[df$agegroup=="children"],
         df$age_num[df$agegroup=="children"])
cor.test(df$beta[df$agegroup=="adolescents"],
         df$age_num[df$agegroup=="adolescents"])
cor.test(df$beta[df$agegroup=="young adults"],
         df$age_num[df$agegroup=="young adults"])

#############################################################
###################ddm
#############################################################
###################ddm
#############################################################
###################ddm
#############################################################
###################ddm ddm ddm ddm ddm ddm ddm ddm ddm ddm ddm

library(RWiener)
library(tidyverse)
library(tidyr)
library(reshape2)
library(dplyr)

##### function for optimizing ddm parameters via first wiener diffusion passage
many_drifts <- function (x, datlist) {
  loss=0
  for (myVar in 1:length(datlist))
    loss = loss - logLik(wdm(datlist[[myVar]], 
                             alpha=x[myVar], tau=x[myVar+1], beta=x[myVar+2], delta=x[myVar+3]))
  return(loss)
}
p_ddm2 = rep(c(1,0.1,0.5,0),each=1)##parameter initializations

data1.ddm2 = list()
subs = unique(itc$ppnr)
x=x2=c() # empty vectors to later compared observed (x2) & ddm.model-predicted (x) LL responses
for(sub in 1:length(subs)){
  df = itc[itc$ppnr==subs[sub],]
  if(length(unique(df$LL))==1) next # optimization needs multiple response types
  df = df[!df$RT<200,]##get rid of RTs below 200 ms
  q = df$RT/1000##turn RTs into seconds
  resp = ifelse(df$LL==1,"upper","lower")##LL==upper boundary; SS==lower boundary
  
  #parameter estimation
  data1.w = data.frame(q,resp)
  data1.w$resp = factor(data1.w$resp)
  data1.w = as.wiener(data1.w)
  data1.both = list(data1.w)
  data1.ddm2[[sub]] = optim(par=p_ddm2,f=many_drifts,dat=data1.both)
  
  #parameter recovery
  parm_recov_s1 <- cbind(rwiener(n=nrow(df), alpha=data1.ddm2[[sub]]$par[1],
                                 tau=data1.ddm2[[sub]]$par[2],
                                 beta=data1.ddm2[[sub]]$par[[3]],
                                 delta=data1.ddm2[[sub]]$par[4]))

  parm_recov.w <- as.wiener(parm_recov_s1)
  x[sub] <- mean(parm_recov.w$resp=="upper")
  x2[sub] = mean(df$LL,na.rm=T)
  if(length(unique(x))==1 | length(unique(x2))==1) next 
  print(sub)
}

####################### parameter recovery plot
X = data.frame(observed = x2,predicted = x)
sqrt(mean((X$observed-X$predicted)^2,na.rm=T))
p00 = ggplot(X,aes(x=predicted,y=observed))+
  geom_smooth(method="lm",col="black")+
  geom_point(shape=21)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),
        axis.title=element_text(size=15))+
  theme(strip.text.x = element_text(size = 15))+
  theme(legend.title=element_text(size=15),
        legend.text = element_text(size=15))+
  ylab("observed proportion of LL responses")+
  xlab("predicted proportion of LL responses")
p00
ggsave(plot=p00,filename="p00.png",units="px",width=2400,height=2400,
       path="~/Dropbox (University of Oregon)/ITC_dev/plots/")

#################### parameter evaluation plots
pars = matrix(0,length(data1.ddm2),4)
for(i in 1:length(data1.ddm2)){
  if(is.null(data1.ddm2[[i]])){
    pars[i,] = NA
    next
  }
  pars[i,1] = data1.ddm2[[i]]$par[1]
  pars[i,2] = data1.ddm2[[i]]$par[2]
  pars[i,3] = data1.ddm2[[i]]$par[3]
  pars[i,4] = data1.ddm2[[i]]$par[4]
}
pars = as.data.frame(pars)
df = data.frame(sub = S,
                agegroup = ages_fact,
                age_num = ages_cont,
                beta = par_est[,1],
                alpha = par_est[,2],
                LL = ll)

##alpha = threshold; tau = nondecision time; delta = drift rate
pars$agegroup = df$agegroup

ddm_df = data.frame(drift_rate = pars$V4,
                    threshold = pars$V1,
                    bias = pars$V3,
                    nondecision = pars$V2,
                    beta = df$beta,
                    alpha = df$alpha,
                    agegroup = df$agegroup,
                    agenum = df$age_num,
                    LL = df$LL,
                    sub = df$sub)
mean(ddm_df$beta);sd(ddm_df$beta)
idx = ddm_df$beta<mean(ddm_df$beta)+2*sd(ddm_df$beta)
idx2 = ddm_df$drift_rate<mean(ddm_df$drift_rate,na.rm=T)+2*sd(ddm_df$drift_rate,na.rm=T)
idx3 = ddm_df$alpha<mean(ddm_df$alpha)+2*sd(ddm_df$alpha)
ddm_df$agegroup=factor(ddm_df$agegroup,levels=c("children","adolescents","young adults"))


m1b = lm(agenum~drift_rate+threshold+nondecision+bias,ddm_df)
summary(m1b)

ddm_df = ddm_df[idx==T&idx==T&idx3==T,]
m3 = lm(beta~drift_rate*agenum,ddm_df)
summary(m3)

cor.test(ddm_df$drift_rate,ddm_df$beta)
cor.test(ddm_df$drift_rate[ddm_df$agegroup=="children"],
         ddm_df$beta[ddm_df$agegroup=="children"])
cor.test(ddm_df$drift_rate[ddm_df$agegroup=="adolescents"],
         ddm_df$beta[ddm_df$agegroup=="adolescents"])
cor.test(ddm_df$drift_rate[ddm_df$agegroup=="young adults"],
         ddm_df$beta[ddm_df$agegroup=="young adults"])

cor.test(ddm_df$drift_rate,ddm_df$alpha)
cor.test(ddm_df$drift_rate[ddm_df$agegroup=="children"],
         ddm_df$alpha[ddm_df$agegroup=="children"])
cor.test(ddm_df$drift_rate[ddm_df$agegroup=="adolescents"],
         ddm_df$alpha[ddm_df$agegroup=="adolescents"])
cor.test(ddm_df$drift_rate[ddm_df$agegroup=="young adults"],
         ddm_df$alpha[ddm_df$agegroup=="young adults"])

###correlate drift rates and reward sensitivities
p4 = ggplot(ddm_df[idx==T&idx2==T,],aes(x=drift_rate,y=beta,col=agegroup,fill=agegroup))+
  geom_smooth(method="lm")+
  scale_y_continuous(limits = c(),expand= c(0,0))+
  ylab(expression(paste("reward sensitivity (",beta,")")))+
  xlab("drift rate")+
  scale_color_manual(values=c("#07E34F","#FF9800","#8C19FF"))+
  scale_fill_manual(values=c("#07E34F","#FF9800","#8C19FF"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),
        axis.title=element_text(size=15))+
  theme(strip.text.x = element_text(size = 15))+
  theme(legend.title=element_text(size=15),
        legend.text = element_text(size=15))+
  theme(legend.position = c(0.75, 0.8))
p4
### correlate drift rates and choice stochasticities
p5 = ggplot(ddm_df[idx3==T&idx2==T,],aes(x=drift_rate,y=alpha,col=agegroup,fill=agegroup))+
  geom_smooth(method="lm")+
  scale_y_continuous(limits = c(),expand= c(0,0))+
  ylab(expression(paste("temperature (",alpha,")")))+
  xlab("drift rate")+
  scale_color_manual(values=c("#07E34F","#FF9800","#8C19FF"))+
  scale_fill_manual(values=c("#07E34F","#FF9800","#8C19FF"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),
        axis.title=element_text(size=15))+
  theme(strip.text.x = element_text(size = 15))+
  theme(legend.title=element_text(size=15),
        legend.text = element_text(size=15))+
  theme(legend.position = c(0.25, 0.85))
p5
p45 = ggpubr::ggarrange(p4,p5,nrow=1,labels=c("A","B"))
ggsave(plot=p45,filename="p45.png",units="px",width=5400,height=2400,
       path="~/Dropbox (University of Oregon)/ITC_dev/plots/")


ddm2 = data.frame(pe = c(pars[,1],pars[,2],pars[,3],pars[,4]),
                  par = rep(c("thr","non","bias","drift"),each=nrow(pars)),
                  agegroup = rep(df$agegroup,4),
                  sub = rep(df$sub,4))
ddm2$agegroup=factor(ddm2$agegroup,levels=c("children","adolescents","young adults"))

ddm_pe1=ggplot(ddm2[ddm2$par=="thr",],aes(x=agegroup,y=pe,fill=agegroup))+
  geom_jitter(width=0.2,shape=21)+
  geom_boxplot(notch=T,outlier.shape=NA)+
  scale_color_manual(values=c("#07E34F","#FF9800","#8C19FF"))+
  scale_fill_manual(values=c("#07E34F","#FF9800","#8C19FF"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),
        axis.title=element_text(size=15))+
  theme(strip.text.x = element_text(size = 15))+
  theme(legend.title=element_text(size=15),
        legend.text = element_text(size=15))+
  theme(legend.position = "None")+
  ylab("threshold")+xlab("")

ddm_pe2=ggplot(ddm2[ddm2$par=="non",],aes(x=agegroup,y=pe,fill=agegroup))+
  geom_jitter(width=0.2,shape=21)+
  geom_boxplot(notch=T,outlier.shape=NA)+
  scale_color_manual(values=c("#07E34F","#FF9800","#8C19FF"))+
  scale_fill_manual(values=c("#07E34F","#FF9800","#8C19FF"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),
        axis.title=element_text(size=15))+
  theme(strip.text.x = element_text(size = 15))+
  theme(legend.title=element_text(size=15),
        legend.text = element_text(size=15))+
  theme(legend.position = "None")+
  ylab("nondecision time")+xlab("")

ddm_pe3=ggplot(ddm2[ddm2$par=="bias",],aes(x=agegroup,y=pe,fill=agegroup))+
  geom_jitter(width=0.2,shape=21)+
  geom_boxplot(notch=T,outlier.shape=NA)+
  scale_color_manual(values=c("#07E34F","#FF9800","#8C19FF"))+
  scale_fill_manual(values=c("#07E34F","#FF9800","#8C19FF"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),
        axis.title=element_text(size=15))+
  theme(strip.text.x = element_text(size = 15))+
  theme(legend.title=element_text(size=15),
        legend.text = element_text(size=15))+
  theme(legend.position = "None")+
  ylab("bias")+xlab("")

ddm_pe4=ggplot(ddm2[ddm2$par=="drift",],aes(x=agegroup,y=pe,fill=agegroup))+
  geom_jitter(width=0.2,shape=21)+
  geom_boxplot(notch=T,outlier.shape = NA)+
  scale_color_manual(values=c("#07E34F","#FF9800","#8C19FF"))+
  scale_fill_manual(values=c("#07E34F","#FF9800","#8C19FF"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),
        axis.title=element_text(size=15))+
  theme(strip.text.x = element_text(size = 15))+
  theme(legend.title=element_blank(),
        legend.text = element_text(size=15))+
  theme(legend.position = c(0.2, 0.95))+
  ylab("drift rate")+xlab("")
ddm1234=ggpubr::ggarrange(ddm_pe2,
                  ddm_pe1,
                  ddm_pe3,
                  ddm_pe4,
                  nrow=2,ncol=2,
                  labels = c("A","B","C","D"))
ggsave(plot=ddm1234,filename="ddm1234.png",units="px",width=2700,height=2400,
       path="~/Dropbox (University of Oregon)/ITC_dev/plots/")

