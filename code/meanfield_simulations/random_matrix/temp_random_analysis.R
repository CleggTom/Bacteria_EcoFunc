#load libraries
library(tidyverse)

#set up expectations
nSp <- 100

p <- read.csv("results/meanfield_simulations/random_matrix/temp_meanfield_sims/params.csv",header = F) %>%
  spread("V1","V2")

k <- 8.6173303e-5
Tref <- 1 / (k*285.0)
Tdiff <- 1 / (k*seq(275.0,295.0)) - Tref

expected_results <- data.frame(Temp = 1:length(Tdiff),
                               mean = 1:length(Tdiff),
                               var = 1:length(Tdiff),
                               tot = 1:length(Tdiff),
                               flux = 1:length(Tdiff))

for(i in 1:length(Tdiff)){
  uU<-log(p$U0)-p$uEU*Tdiff[i];   sU<-abs(p$sEU * Tdiff[i])
  uR<-log(p$R0)-p$uER*Tdiff[i];   sR<-abs(p$sER * Tdiff[i])
  ua<-log(p$a0)-p$uEa*Tdiff[i];   sa<-abs(p$sEa * Tdiff[i])

  uU_real <- exp(uU + (sU^2/2))
  sU_real <- exp(2*uU + (sU^2))*(exp((sU^2)) - 1)

  uR_real <- exp(uR + (sR^2/2))
  sR_real <- exp(2*uR + (sR^2))*(exp((sR^2)) - 1)

  ua_real <- exp(ua + (sa^2/2))
  sa_real <- exp(2*ua + (sa^2))*(exp((sa^2)) - 1)

  uP <- ua_real * 100
  sP <- sa_real * 100

  expected_results$Temp[i] <- Tdiff[i]
  expected_results$mean[i] <- (uU_real - uR_real)/(uP + 1)
  expected_results$var[i] <- sU_real + sR_real + expected_results$mean[i]*sP
  expected_results$tot[i] <- (uU_real - uR_real - expected_results$mean[i]) / ua_real
  expected_results$flux[i] <- expected_results$tot[i] * uR_real

}

expected_results %>%
  gather("param","value",-Temp) %>%
  ggplot(aes(x=Temp,y=value))+
  geom_point()+
  facet_wrap(~param,scales = 'free_y')

#get simulation results
filenames <- list.files("results/meanfield_simulations/random_matrix/temp_meanfield_sims")
filenames <- filenames[grep("^Temp",filenames,)]


results <- data.frame(Temp = 1:length(filenames),
                      sim = 1:length(filenames),
                      mean = 1:length(filenames),
                      var = 1:length(filenames),
                      tot = 1:length(filenames),
                      flux = 1:length(filenames))

end_bio <- data.frame(matrix(ncol = length(filenames),nrow = 100))

for(i in 1:length(filenames)){
  files <- str_split(filenames[i],pattern = "_",simplify = T)
  df <- read.csv(paste("results/meanfield_simulations/random_matrix/temp_meanfield_sims/",filenames[i],sep = '/'),header = F)
  R <-  read.csv(paste("results/meanfield_simulations/random_matrix/temp_meanfield_sims",paste0("Resp_",filenames[i]),sep = '/'),header = F)

  endbio <- as.numeric(df[nrow(df),])
  end_bio[,i] <- endbio

  #get values
  avg <- mean(endbio)
  var <- var(endbio)


  Tdiff <- (1 / (k*as.numeric(files[2]))) - Tref


  results$Temp[i] <- Tdiff
  results$sim[i] <- as.numeric(str_split(files[4],pattern = ".c")[[1]][1])
  results$mean[i] <- avg
  results$var[i] <- var
  results$tot[i] <- sum(endbio)
  results$flux[i] <- sum(R*endbio)

  print(i)

  end_bio[,i] <- as.numeric(df[nrow(df),])
  }

plot_sim <- function(i){
 read.csv(paste("results/meanfield_simulations/random_matrix/temp_meanfield_sims/",filenames[i],sep = '/'),header = F) %>%
    mutate(t = 1:n()) %>%
    gather("species","biomass",-t)%>%
    ggplot(aes(x=t,y=biomass,group=species))+
    geom_line()
}

#looks pretty good!
results %>%
  gather("param",'value',c(-sim,-Temp)) %>%
  ggplot(aes(x=Temp,y=value))+
  geom_point()+
  facet_wrap(~param,scales = 'free_y')+
  geom_line(data = expected_results %>% gather("param","value",-Temp) ,colour='red')
