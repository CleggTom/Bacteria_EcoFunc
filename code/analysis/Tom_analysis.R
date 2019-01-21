library(tidyverse)
library(fitdistrplus)

k <- 8.6173303e-5
Tvec <- seq(1/(273.15*k),1/(295.15*k),length.out = 100)
Tnorm <- 1/(k*273.15)
df <- read_csv("data/Tom/summary_R_biomass.csv")

#Plotting to check fits
plot_tpc <- function(i,Tvec){
  assertthat::assert_that(i<= length(unique(df$strain)),msg = "strain index wrong")
  strain_name <- unique(df$strain)[i]
  
  params <- df %>%
    filter(strain == strain_name) %>%
    select(strain,trait,E_sch,B0_sch)
    
  g_p <- params[params$trait == "Growth Rate",]
  r_p <- params[params$trait == "Respiration Rate",]
  
  g <- g_p$B0_sch * exp(-g_p$E_sch * (Tvec-Tnorm))
  r <- r_p$B0_sch * exp(-r_p$E_sch * (Tvec-Tnorm))
  
  B0_lines <- data.frame(key = c('g','r'), B0 = c(g_p$B0_sch , r_p$B0_sch))
  
  data.frame(g,r,Tvec) %>%
      gather(key,value,-Tvec) %>%
      ggplot(aes(x=Tvec,y=value))+
        geom_line()+
        facet_wrap(~key,scales = 'free_y')+
       # scale_y_log10()+
        geom_hline(data=B0_lines,aes(yintercept=B0),color='red')
  
}

plot_tpc(4,Tvec)

#covarience
df %>%
  filter(temps_before_peak > 3) %>%
  ggplot(aes(x=(E_sch),y=log(B0_sch)))+
  geom_point()+
  facet_wrap(~trait)

df %>%
  filter(temps_before_peak > 3) %>%
  group_by(trait) %>%
  summarise(cor = cor((.$E_sch),log(.$B0_sch)))

#df %>%
#  filter(strain != '50_RT_03') %>%
#  mutate(B0_sch = log(B0_sch),E_sch = log(E_sch)) %>%
#  gather('param','value',E_sch,B0_sch) %>%
#  group_by(param,trait) %>%
#  summarise(mean = median(value))

df %>%
  filter(strain != '50_RT_03') %>%
  filter(log(B0_sch) > -10) %>%
  mutate(B0_sch_log = log(B0_sch),E_sch = (E_sch)) %>%
  gather('param','value',E_sch,B0_sch_log) %>%
  ggplot(aes(x=value,group=trait,fill=trait))+
    geom_histogram(position='identity',binwidth = 0.5)+
    facet_grid(trait~param,scales = 'free_x')

cov_E_B0 <- df %>%
  filter(temps_before_peak > 3) %>%
  filter(trait == 'Respiration Rate') %>%
  summarise(cov = cov(E_sch,log(B0_sch)))


a <- df %>%
  filter(temps_before_peak > 3) %>%
  filter(trait == 'Respiration Rate') %>%
  select(E_sch,B0_sch) %>%
  mutate(B0_sch = log(B0_sch), l_B = B0_sch - E_sch*1) %>%
  gather() %>%
  group_by(key) %>%
    summarise(mean = mean(value),
              var = var(value))
#pred mean
a$mean[1] - a$mean[2]

#pred var
a$var[1] + a$var[2] - 2*cov_E_B0$cov

df %>%
  filter(temps_before_peak > 3) %>%
  filter(trait == 'Respiration Rate') %>%
  select(E_sch,B0_sch) %>%
  mutate(B0_sch = log(B0_sch), l_B = B0_sch - E_sch*1) %>%
  ggplot(aes(x=exp(l_B)))+
  geom_histogram(bins=10)
