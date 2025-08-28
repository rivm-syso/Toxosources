# Load required libraries
library( tidyverse )
library( here )
library( readxl )
library( rstan )
library( tidybayes )
library( patchwork )
library( knitr )
library( furrr )
library( writexl )
library( scales )
library( glue )

setwd(here())

options(mc.cores = 1)
rstan_options(auto_write = TRUE)
plan(sequential)
options(warn=0)

source( "functions.R")

modeltype <- "SIS"
filepath <- here("output", str_c("_", modeltype ))
if( !dir.exists(filepath)) dir.create(filepath, recursive = TRUE)


df_animals <- read_animals() %>% 
  mutate( region = case_when(country == "the United Kingdom" ~ "United Kingdom",
                             .default = region)) %>% # Separate UK from West region
  mutate( age_high= ifelse( age_low==age_high, age_high*1.1, age_high )) %>%
  mutate( age_high= ifelse( age_low==age_high & age==age_high, age_high*1.1, age_high )) %>%
  mutate( age_low = ifelse( age_low==age_high & age==age_low,  age_low*0.9, age_low )) %>%
  mutate( age = ifelse( age_low==age, (age_high+age_low)/2, age )) %>%
  mutate( age = ifelse( age_high==age, (age_high+age_low)/2, age )) %>%
  ungroup() %>% 
  droplevels()

lm(prev ~ 1 + pop_group, data=df_animals) %>% summary()

#
# Stan model for prevalence in humans and regions, censored age
#

m_prevalence <- stan( file="prevalence.stan",
                      iter=3000,
                      seed=3141,
                      data=compose_data( df_animals,
                                         use_sis_model=(modeltype=="SIS")))

my_pars <- c( "lambda_region[1]", "lambda_region[2]",
              "lambda_baseline",
              "sigma_lambda_region",
              "lp__")

if( modeltype!="SI") my_pars <- c(my_pars, "gamma_baseline")

png( file.path(filepath,"traceplot.png" ) )
print(traceplot(m_prevalence, pars = my_pars, inc_warmup=FALSE ))
dev.off()


df_prevalence <- m_prevalence %>% 
  recover_types( df_animals ) %>% 
  spread_draws( lambda_baseline,
                lambda_region[region],
                sigma_lambda_region,
                gamma_baseline[0] ) %>% 
  mutate( across( contains("gamma"), as.numeric )) %>%
  mutate( lambda_region = sigma_lambda_region * lambda_region,
          lambda = exp(lambda_baseline + lambda_region),
          gamma = exp( gamma_baseline  ) ) %>% 
  #  Rescale lambda_ and gamma_
  # variables because of how they are in the formula.
  mutate( across( starts_with(c("lambda_", "gamma_" )), exp ) ) %>% 
  ungroup()

# In the stan code lambda = exp(lambda_baseline + lambda_region );
# Above, lambda_baseline and the rest is exponentiated!

saveRDS( df_prevalence, file= file.path(filepath, "df_prevalence.RData"))
saveRDS( m_prevalence, file= file.path(filepath, "m_prevalence.RData"))

m_prevalence <- readRDS( file= file.path(filepath, "m_prevalence.RData"))
df_prevalence <- readRDS( file.path(filepath, "df_prevalence.RData"))

df_prevalence %>%
  group_by( region ) %>% 
  mean_qi( lambda_region )  %>% 
  arrange( lambda_region ) 

df_prevalence %>%
  mean_qi( lambda_baseline )  

df_prevalence %>%
  mean_qi( gamma_baseline ) 

df_age <- m_prevalence %>% 
  recover_types( df_animals ) %>% 
  spread_draws( age_raw[i] ) %>% 
  group_by(i) %>% 
  summarize( age_raw=mean(age_raw), .groups="drop") %>% 
  cbind( ind=df_animals$ind ) %>% 
  left_join( df_animals ) %>% 
  mutate( agedist = (age_raw*(age_high-age_low) + age_low )) %>% 
  select( ind, agedist )


df_age %>% 
  left_join( df_animals ) %>% 
  mutate( ind=as.factor(ind),
          ind=fct_reorder(ind, age)) %>% 
  ggplot() +
    geom_errorbar( aes(xmin=age_low,xmax=age_high, y=ind) ) +
    geom_point( aes( x=age, y=ind), color="black", size=2 ) +
    geom_point( aes( x=agedist, y=ind), color="blue", size=2 ) 
ggsave(file.path(filepath, "age_shift.png"), 
       width=30, height=30, units="cm" )


df_prevalence %>% 
  ggplot() +
  stat_halfeye( aes(sigma_lambda_region), fill="green" ) +
  ggtitle( "Hyperparameters lambda region")  

ggsave( file.path(filepath, "hyperparameters.png"), 
              width=15, height=10, units="cm" )

plot_posteriors( df_prevalence )
ggsave( file.path(filepath, "posteriors.png"), 
        width=25, height=20, units="cm", dpi=600 )   
ggsave( file.path(filepath, "posteriors.svg"), 
        width=25, height=20, units="cm", dpi=600 )   


l_plots <- df_prevalence %>% 
  group_by( region ) %>% 
  expand_grid( x=seq(0, 100, length.out=20)) %>% 
  mutate(p=sis(lambda, gamma, x),
         region = fct_relevel(region, "West", "East", "North", "Southwest", "Southeast")) %>% 
  group_by( region ) %>% 
  group_map( \(y, gr){
    df_animals %>%
      filter( region == gr[[1,1]] ) %>% 
      left_join( df_age, by="ind" ) %>% 
        ggplot( ) + 
          scale_x_continuous( "Age", limits=c(  0, 100 )) +
          scale_y_continuous( "Prevalence",  breaks=seq(0,1,by=0.2), labels=scales::percent ) +
          coord_cartesian( ylim=c(0,1) ) +
          # Data
          geom_segment( aes( x=age_low, xend=age_high, y=prev, yend=prev), color="gray", alpha=0.5) +
          geom_point( aes( x=age, y=prev, size=log10(n_tot) ), color="gray" )+
          scale_size(guide = 'none') +
          scale_fill_discrete(guide = "none" ) +
          # Model fit
          stat_lineribbon( data=y, 
                           mapping=aes(x, p ), 
                           .width=0.95, point_interval=mean_hdci, alpha=0.3 ) +
          geom_point( aes(x=agedist, y=prev, color=country ) ) +
          theme_minimal() +
          theme( legend.position = "bottom", legend.title = element_blank(),
                 legend.key.spacing.y = unit(0.01,"cm"),
                 legend.box.spacing = unit(0,"cm")) +
          guides(color=guide_legend(nrow=2,byrow=TRUE)) +
          ggtitle( glue("{gr[[1]]}"))})

wrap_plots( l_plots, nrow=2 )
ggsave(file.path(filepath, "prevalence_curves.png"), 
       width=30, height=30, units="cm", dpi=600)
ggsave(file.path(filepath, "prevalence_curves.svg"), 
       width=30, height=30, units="cm", dpi=600)

# With or without gamma
df_prevalence %>% 
  group_by( region ) %>% 
  expand_grid( x=seq(0, 100, length.out=20)) %>% 
  mutate(p=sis(lambda, gamma, x),
         p2=sis(lambda, 0, x)) %>%
  ggplot() +
  stat_lineribbon( aes(x, p ), 
                   .width=0.95, point_interval=mean_hdci, alpha=0.3 ) +
  stat_lineribbon( aes(x, p2 ), 
                   .width=0.95, point_interval=mean_hdci, alpha=0.3, fill="blue" )
  
  
  
rbind(
  df_prevalence %>% 
    group_by( region ) %>% 
    expand_grid( age=seq(0, 1.1* 99, length.out=20)) %>%
    mutate( age_group = case_when(age<25 ~ "0-25",
                                   age<50 ~ "25-50",
                                   .default = "50+")),
  df_prevalence %>% 
    group_by( region ) %>% 
    expand_grid( age=seq(0, 1.1* 99, length.out=20)) %>%
    mutate( age_group = "Overall" )
  ) %>% 
  mutate(p=sis(lambda, gamma, age) ) %>%
  group_by(age_group, region) %>%
  summarise(mean_val = mean(p), 
            q025 = quantile(p, probs = .025),
            q975 = quantile(p, probs = .975)) %>% 
  mutate( mean = glue( "{round(mean_val*100)}%" ),
          ci = glue( "{round(q025*100)}-{round(q975*100)}%") ) %>% 
  select( region, age_group, mean, ci ) %>% 
  pivot_wider( names_from=age_group, values_from = c(mean, ci ) ) %>%
  knitr::kable( format="html",digits=2) %>% 
  cat( file=here("output", "table2.html") )

# Linear model

df_lm <- df_prevalence %>%
  filter( region=="West" ) %>% 
  expand_grid( age=seq(0, 90, length.out=200)) %>% 
  mutate(p=sis(lambda, gamma, age) ) %>% 
  select( age, p )

lm( p ~ age-1, data=df_lm %>% filter( age <= 17 ))
lm( p ~ age, data=df_lm %>% filter( age > 17 ))


df_quantiles <- df_prevalence %>% 
  mutate( avg = 1/lambda,
          q10 = qexp( 0.1, lambda),
          q50 = qexp( 0.5, lambda),
          q90 = qexp( 0.9, lambda)) 

rbind( 
  df_quantiles %>% 
    group_by( region ) %>% 
    summarize( across( c(avg,starts_with("q")), mean ), .groups="drop" ) %>%
    select( region, avg, starts_with("q") ),
  df_quantiles %>% 
    summarize( across( c(avg,starts_with("q")), mean ), .groups="drop" ) %>%
    mutate( region="Overall") %>% 
    select( region, avg, starts_with("q") )) %>% 
    knitr::kable( format="html",digits=2) %>% 
    cat( file=here("output", "table_quantiles.html") )
