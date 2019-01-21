#load libraries
library(tidyverse)

#set up expectations
nSp <- 100

p <- read.csv("results/meanfield_simulations/random_matrix/interaction_random_sims/params.csv",header = F) %>%
  spread("V1","V2")

k <- 8.6173303e-5
Tref <- 1 / (k*285.0)
Tdiff <- 1 / (k*seq(275.0,295.0,length.out = 10)) - Tref
a_seq <- seq(0.1,10.0,length.out = 10) / 150

expected_results <- data.frame(Temp = 1:(length(Tdiff)*length(a_seq)),
                               a    = 1:(length(Tdiff)*length(a_seq)),
                               mean = 1:(length(Tdiff)*length(a_seq)),
                               var  = 1:(length(Tdiff)*length(a_seq)),
                               tot = 1:(length(Tdiff)*length(a_seq)),
                               flux = 1:(length(Tdiff)*length(a_seq)))
for(a in 1:length(a_seq)){
  for(i in 1:length(Tdiff)){

    uU<-log(p$U0)-p$uEU*Tdiff[i];   sU<-abs(p$sEU * Tdiff[i])
    uR<-log(p$R0)-p$uER*Tdiff[i];   sR<-abs(p$sER * Tdiff[i])
    ua<-log(a_seq[a])-p$uEa*Tdiff[i];   sa<-abs(p$sEa * Tdiff[i])

    uU_real <- exp(uU + (sU^2/2))
    sU_real <- exp(2*uU + (sU^2))*(exp((sU^2)) - 1)

    uR_real <- exp(uR + (sR^2/2))
    sR_real <- exp(2*uR + (sR^2))*(exp((sR^2)) - 1)

    ua_real <- exp(ua + (sa^2/2))
    sa_real <- exp(2*ua + (sa^2))*(exp((sa^2)) - 1)

    uP <- ua_real * 150
    sP <- sa_real * 150

    expected_results$Temp[i + (a-1)*length(Tdiff)] <- Tdiff[i]
    expected_results$a[i + (a-1)*length(Tdiff)] <- a_seq[a]
    expected_results$mean[i + (a-1)*length(Tdiff)] <- (uU_real - uR_real)/(uP + 1)
    expected_results$var[i + (a-1)*length(Tdiff)] <- sU_real + sR_real + (sP*(uU_real - uR_real)/(uP + 1))

    expected_results$tot[i + (a-1)*length(Tdiff)] <- (uU_real - uR_real - expected_results$mean[i + (a-1)*length(Tdiff)]) / ua_real
    expected_results$flux[i + (a-1)*length(Tdiff)] <- expected_results$tot[i + (a-1)*length(Tdiff)] * uR_real

  }
}

expected_results %>%
  gather('param','value',c(-Temp,-a)) %>%
  ggplot(aes(x=Temp,y=value,colour=a,group=a))+
  geom_line()+
  facet_grid(param~a,scales = 'free_y')



#get simulation results
filenames <- list.files("results/meanfield_simulations/random_matrix/interaction_random_sims/")
filenames <- filenames[grep("^Int",filenames)]


results <- data.frame(Temp = 1:length(filenames),
                      a = 1:length(filenames),
                      sim = 1:length(filenames),
                      mean = 1:length(filenames),
                      var = 1:length(filenames),
                      tot = 1:length(filenames),
                      flux = 1:length(filenames))

end_bio <- data.frame(matrix(ncol = length(filenames),nrow = 150))
R_df <- data.frame(matrix(ncol = length(filenames),nrow = 150))


for(i in 1:length(filenames)){
  files <- str_split(filenames[i],pattern = "_",simplify = T)
  df <- read.csv(paste("results/meanfield_simulations/random_matrix/interaction_random_sims/",filenames[i],sep = '/'),header = F)
  R <-  read.csv(paste("results/meanfield_simulations/random_matrix/interaction_random_sims",paste0("Resp_",filenames[i]),sep = '/'),header = F)

  endbio <- as.numeric(df[nrow(df),])
  end_bio[,i] <- endbio

  #get values
  avg <- mean(endbio)
  var <- var(endbio)

  Tdiff <- (1 / (k*as.numeric(files[4]))) - Tref

  results$Temp[i] <- Tdiff
  results$a[i] <- as.numeric(files[2])
  results$sim[i] <- as.numeric(str_split(files[6],pattern = ".c")[[1]][1])
  results$mean[i] <- avg
  results$var[i] <- var
  results$tot[i] <- sum(endbio)
  results$flux[i] <- sum(R*endbio)

  print(i)

  end_bio[,i] <- as.numeric(df[nrow(df),])
  R_df[,i] <- R
  }

plot_sim <- function(i){
 read.csv(paste("results/meanfield_simulations/random_matrix/interaction_meanfield_sims",filenames[i],sep = '/'),header = F) %>%
    mutate(t = 1:n()) %>%
    gather("species","biomass",-t)%>%
    ggplot(aes(x=t,y=biomass,group=species))+
    geom_line()
}

merge(results,expected_results,by = c('Temp','a')) %>%
  rename(sim_mean=mean.x,sim_var=var.x,sim_tot=tot.x,sim_flux=flux.x,
         exp_mean=mean.y,exp_var=var.y,exp_tot=tot.y,exp_flux=flux.y) %>%
  gather("param","value",-Temp,-sim,-a) %>%
  separate(param,c("type",'param')) %>%
  spread(type,value) %>%
  ggplot(aes(x=a))+
  facet_grid(param~Temp,scales = 'free')+
  geom_point(aes(y=sim,colour=Temp))+
  geom_line(aes(y=exp,group=Temp),color='red')+
  scale_y_log10()

model <- lm(log(flux)~Temp*a,data = results)

a <- 0.0
Temp <- c(-2,2)

b <- coef(model)[[1]] + coef(model)[[2]]*(Temp) + coef(model)[[3]]*a + coef(model)[[4]]*Temp*a

diff(b)/4


