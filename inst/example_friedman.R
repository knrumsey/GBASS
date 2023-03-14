f<-function(x){
  10*sin(pi*x[,1]*x[,2])+20*(x[,3]-.5)^2+10*x[,4]+5*x[,5]
}
n<- 500 # number of observations
x<-matrix(runif(n*10),n,10) #10 variables, only first 5 matter
xx<-matrix(runif(1000*10),1000,10)
y<-f(x)

#profvis({gmod<-gbass(x,y, nmcmc=1000, nburn=901, v_prior = build_prior('GIG',c(-15,0,30),.5),w_prior = build_prior('GIG',c(0,0,0),.5))})
gmod<-gbass(x,y, nmcmc=1000, nburn=901, v_prior = build_prior('GIG',c(-15,0,30),.5),w_prior = build_prior('GIG',c(0,0,0),.5))


#gmod2<-gbass(x,y,v_prior = build_prior('GIG',c(-5,0,1),.5),w_prior = build_prior('GIG',c(0,0,0),.5))
gmod3<-gbass(x,y,v_prior = build_prior('GIG',c(5,10,0),.5),w_prior = build_prior('GIG',c(0,0,0),.5))
#gmod4<-gbass(x,y,v_prior = build_prior('GIG',c(3,6,0),.5),w_prior = build_prior('GIG',c(0,0,0),.5))
#gmod5<-gbass(x,y,v_prior = build_prior('GIG',c(300,600,0),.5),w_prior = build_prior('GIG',c(0,0,0),.5))
#gmod6<-gbass(x,y,v_prior = build_prior('GIG',c(.3,.6,0),.5),w_prior = build_prior('GIG',c(0,0,0),.5))
#gmod7<-gbass(x,y,v_prior = build_prior('GIG',c(2.5,15/pi^2,0),.5),w_prior = build_prior('GIG',c(0,0,0),.5))
#gmod8<-gbass(x,y,v_prior = build_prior('GIG',c(1,2*5,0),.5),w_prior = build_prior('GIG',c(0,0,0),.5))
#gmod9<-gbass(x,y,v_prior = build_prior('GIG',c(-.5,2/4,2,0),.5),w_prior = build_prior('GIG',c(0,0,0),.5))
bmod<-BASS::bass(x,y,h1=1,h2=.001)

sqrt(mean((f(xx)-rowMeans(predict(gmod,xx)))^2))
#sqrt(mean((f(xx)-rowMeans(predict(gmod2,xx)))^2))
sqrt(mean((f(xx)-rowMeans(predict(gmod3,xx)))^2))
#sqrt(mean((f(xx)-rowMeans(predict(gmod4,xx)))^2))
#sqrt(mean((f(xx)-rowMeans(predict(gmod5,xx)))^2))
#sqrt(mean((f(xx)-rowMeans(predict(gmod6,xx)))^2))
#sqrt(mean((f(xx)-rowMeans(predict(gmod7,xx)))^2))
#sqrt(mean((f(xx)-rowMeans(predict(gmod8,xx)))^2))
#sqrt(mean((f(xx)-rowMeans(predict(gmod9,xx)))^2))
sqrt(mean((f(xx)-colMeans(predict(bmod,xx)))^2))

plot(gmod3$v[,2],type='l')

matplot(gmod$v,type='l')
ts.plot(gmod$w)

plot(bmod)
