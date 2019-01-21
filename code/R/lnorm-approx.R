random_com <- function(N,p,Tdiff){
  U <- rlnorm(N,log(p$U0) - p$EU*Tdiff, abs(p$sU*Tdiff))
  R <- rlnorm(N,log(p$R0) - p$ER*Tdiff, abs(p$sR*Tdiff))
  
  a_sample <- rlnorm(N,log(0.01)-0*Tdiff, abs(0.1*Tdiff))
  P <- replicate(N,sum(sample(a_sample,100)))
  
  C <- U - R - P*(mean(U)-mean(R))/(mean(P)+1)
  R_flux <- C*R
  
  #approx the C distribution
  u_app <- log(mean(C) / sqrt( (var(C)/mean(C)^2)+1 ) )
  v_app <- sqrt(log( (var(C)/mean(C)^2)+1 ))
  
  #get the product
  u_RC <- log(p$R0) - p$ER*Tdiff + u_app
  sR <- abs(p$sR*Tdiff)
  sC <- v_app
  s_RC <- sR^2 + sC^2 + 2*cor(C,R)*sR*sC
  
  return(data.frame(RC_approx = rlnorm(N,u_RC,sqrt(s_RC)), R_real = R_flux))
  
}

p <- list(U0 = 10,EU = 0.32, sU = 0.1,
          R0 = 1.,ER = 0.64, sR = 0.1,
          a0 = 0.01,Ea=0.0,sa=0.1)

images <- list(length=100)
Tvec <- seq(2,-2,length.out = 100)

for(i in 1:100){
images[[i]] <- random_com(1e4,p,Tvec[i]) %>%
  gather() %>%
  ggplot(aes(x=value,fill=key))+
  geom_histogram(position='identity',alpha=0.5,stat='density')+
  xlim(0,150)+
  ggtitle(paste(Tvec[i]))
}

setwd("Desktop/")
pdf("plots.pdf")
for (i in 1:100) {
  print(images[[i]])
}
dev.off()

