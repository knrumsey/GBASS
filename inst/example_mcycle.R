library(MASS)
library(BASS)
library(parallel)
library(RColorBrewer)
mcycle$x <- BASS:::scale.range(mcycle$times)
xx <- matrix(seq(0,1,by=.01), ncol=1)

mods <- mclapply(c(.1,.25,.5,.75,.9),
                 function(qq) qbass(matrix(mcycle$x),
                                    mcycle$accel,q=qq,
                                    maxInt=1,
                                    w_prior=list(type="GIG", p=0, a=0, b=0, prop_sigma=0.2),
                                    a_lambda=.03, b_lambda=.03,
                                    nmcmc=10000, nburn=8001, thin=4),
                 mc.cores = 3, mc.preschedule = F)

mods2 <- mclapply(c(.1,.25,.5,.75,.9),
                 function(qq) qbass(matrix(mcycle$x),
                                    mcycle$accel,q=qq,
                                    maxInt=1,
                                    w_prior=list(type="GIG", p=0, a=0, b=0, prop_sigma=0.2),
                                    a_lambda=1, b_lambda=1,
                                    nmcmc=10000, nburn=8001, thin=4),
                 mc.cores = 3, mc.preschedule = F)



bob.ross <- brewer.pal(8, "Dark2")
plot(mcycle$x, mcycle$accel, pch=16, xlab="Time", ylab='Acceleration', cex.lab=1.5, ylim=c(-150, 60))
grid()
for(i in 1:5){
  yhat <- rowMeans(predict(mods[[i]], xx))
  yhat2 <- rowMeans(predict(mods2[[i]], xx))
  lines(xx, yhat, lwd=2, col=bob.ross[i])
  lines(xx, yhat2, lty=2, lwd=2, col=adjustcolor(bob.ross[i], alpha.f=0.5))
}
