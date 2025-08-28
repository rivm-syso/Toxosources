read_animals <- function(){
   read_excel(here("data", "humanprevalence_extracteddata_age15-50_only_160524.xlsx"), 
                   trim_ws = TRUE) %>% 
          select( 
            n_pos = "Total number of seropositive participants (n)",
            n_tot =  "Total number of tested participants (N)",
            test = Test,
            country=Country,
            year = `End of sample collection (Year)`,
            pop_group = `population group`,
            region = `Region`,
            age_low = `Age - lower [years]`,
            age_high = `Age - upper [years]`,
            age_best = `Age - most probable [years]`,
            pop_id = population,
            study = `Number of the article`
             ) %>% 
          mutate( pop_id = as.character(pop_id),
                  prev = n_pos/n_tot,
                  age_best = as.numeric( age_best),
                  across( c(test,study, region), as.factor)) %>% 
        select( region, test,year,pop_group,  country,
              age_low, age_high, age=age_best,pop_id,
              prev, study,
              n_pos, n_tot) %>%
        filter( n_tot > 0 ) %>%
        filter( !is.na(n_tot), !is.na(n_pos)) %>%
    group_by( pop_id, region, test, year, pop_group, country ) %>% 
    mutate( population_id = as.factor(as.character(cur_group_id())),
            weight = 1 ) %>% 
    ungroup() %>% 
    mutate( ind =1:n()) %>% 
    select(-pop_id ) %>%
    droplevels()
}

sir <- function(lambda, gamma=0, x){
  lambda*(exp(-lambda*x)-exp( -gamma*x) )/(gamma-lambda)
}

sis <- function(lambda, gamma, x){
  lambda*(1-exp( -(lambda+gamma) * x) )/(lambda+gamma)
}

plot_hyperparameters <- function( df_prevalence ){
 p1 <- df_prevalence %>% 
      ggplot() +
      stat_halfeye( aes(sigma_lambda_region), fill="green" ) +
      ggtitle( "Hyperparameters lambda region")
 p2 <- df_prevalence %>% 
      ggplot() +
      stat_halfeye( aes(sigma_lambda_species), fill="green" ) +
      ggtitle( "Hyperparameters lambda species")
 p3 <- df_prevalence %>% 
   ggplot() +
   stat_halfeye( aes(sigma_lambda_sample_type), fill="green" ) +
   ggtitle( "Hyperparameters lambda sample type")
 p4 <- df_prevalence %>% 
   ggplot() +
   stat_halfeye( aes(sigma_lambda_outdoor), fill="green" ) +
   ggtitle( "Hyperparameters lambda outdoor")
 p5 <- df_prevalence %>% 
   ggplot() +
   stat_halfeye( aes(sigma_gamma_species), fill="green" ) +
   ggtitle( "Hyperparameter sigma_gamma_species")

  if( modeltype=="SI"){  
    reduce( list(p1,p2,p3,p4), `+` ) + plot_layout( ncol=2 )
  }else{
    reduce( list(p1,p2,p3,p4,p5), `+` ) + plot_layout( ncol=2 )
  }
    
}

plot_posteriors <- function( df_prevalence ){

  p1 <- df_prevalence %>% 
    ggplot() +
      stat_halfeye( aes(lambda_baseline), fill="gray" ) +
      scale_x_continuous( "Baseline Force of Infection") +
      ggtitle(expression(lambda["baseline"])) +
      theme_minimal()
  p2 <- df_prevalence %>% 
    ggplot() +
      stat_gradientinterval( aes(lambda_region, y=region ),fill_type = "gradient" ) +
      scale_x_continuous( "Force of infection", limits=c(0,2.5) ) +
      ggtitle(expression(lambda["region"])) +
      theme_minimal()
  p6 <- df_prevalence %>% 
    ggplot() +
      stat_halfeye( aes(gamma_baseline), fill="gray" ) +
      scale_x_continuous( "Reversion", limits=c(0,0.005) ) +
      ggtitle(expression(gamma)) +
      theme_minimal()

  if(modeltype == "SI" ){
    reduce( list(p1,p2), `+` ) + plot_layout( ncol=2  )
  }else{
    reduce( list(p1,p2,p6), `+` ) + plot_layout( ncol=2  )
  }
}

plot_agedist <- function( df_age ){
  df_age %>%
    filter( ind %in% sample(1:max(.$ind), 10 )) %>% 
    ggplot() +
    stat_halfeye( aes(agedist, fill=as.factor(ind)) )
}


# Using estimated age
plot_dynamics_est_age <- function( df_animals, df_age, fun ){
  if( modeltype !="SI") {
    df_prevalence <- df_prevalence %>%
      group_by( region, species, sample_type) %>%  # outdoor ) %>%
      summarize( lambda = mean(lambda_baseline * 
                                      lambda_species * 
                                      lambda_region *
                                      lambda_sample_type #*
                                      #lambda_outdoor 
                               ), 
                  gamma=mean(gamma), .groups="drop" ) %>%
      expand_grid( x=0:20 ) %>%
      mutate(p=fun(lambda, gamma, x) )
  }else{
    df_prevalence <- df_prevalence %>%
      group_by( region, species, sample_type ) %>% #, outdoor ) %>%
      summarize( lambda = mean( exp(lambda_baseline + 
                                      lambda_species + 
                                      lambda_region + 
                                      #lambda_outdoor +
                                      #lambda_population +
                                      lambda_sample_type)), 
                 .groups="drop" ) %>%
      expand_grid( x=0:20 ) %>%
      mutate(p=fun(lambda, 0, x) )
  }
  
  df_animals %>%
    left_join( df_age) %>% 
    mutate( 
      ymax=qbinom( 0.975, size=n_tot, prob=prev )/n_tot,
      ymin=qbinom( 0.025, size=n_tot, prob=prev )/n_tot) %>%
    ggplot( ) + 
    scale_x_continuous( "Age", limits=c(  0,20)) +
    scale_y_continuous( "Prevalence") +
    ggtitle( "The data and fit, estimated age" ) +
    geom_rect( aes( xmin=age_low, xmax=age_high, ymin=ymin, ymax=ymax, color=sample_type ),
               fill=NA ) +
    geom_point( aes(x=agedist, y=prev, color=sample_type ), size=1 )+
    geom_line( data=df_prevalence,
               mapping=aes(x, p, color=sample_type), size=0.5) +
    facet_grid( vars(region), vars(species), scale="free") +
    scale_size(range = c(0, 5))
}


plot_dynamics_by_species <- function( df_prevalence, df_animals, df_age, fun ){

  if( modeltype =="SI") {
    df_prevalence$gamma <- 0
  }
  
  df_prevalence <- df_prevalence %>% 
    group_by( region ) %>% 
    expand_grid( x=seq(0, 1.1* 99, length.out=20)) %>% 
    mutate(p=fun(lambda, gamma, x) )
  
  p <- df_animals %>%
    left_join( df_age, by="ind" ) %>% 
    mutate( across( c(age, age_low, age_high), ~.x*age_max ) ) %>%
    ggplot( ) + 
        scale_x_continuous( "Age", limits=c(  0, 100 )) +
        scale_y_continuous( "Prevalence",  breaks=seq(0,1,by=0.2), labels=scales::percent ) +
        coord_cartesian( ylim=c(0,1) ) +
        # Data
        geom_segment( aes( x=age_low, xend=age_high, y=prev, yend=prev), color="gray", alpha=0.5)+
        geom_point( aes( x=age, y=prev, size=log10(n_tot) ), color="gray" )+
        # Model fit
        stat_lineribbon( data=df_prevalence, 
                   mapping=aes(x, p ), 
                   .width=0.95, point_interval=mean_hdci, alpha=0.3 ) +
        geom_point( aes(x=agedist, y=prev, color=country ) )+
        #geom_label( aes(x=agedist, y=prev, label=study) )+
        facet_wrap( vars(region) ) +
        ggtitle( "The data and fit, estimated age" )
  ggsave(file.path(filepath, "byspecies.png"), plot=p)
  return(p)
}

inits = function() {
  list(
    # sigma_lambda_region = 0.5,
    # sigma_lambda_species = 1,
    # sigma_lambda_sample_type = 0.1,
    # sigma_lambda_outdoor = 1,
    # sigma_gamma_species = array( 1, dim=(modeltype!="SI")),
    # mu_gamma_species = array( -1, dim=(modeltype!="SI")),
    lambda_baseline = 1,
    gamma_baseline = array( 0.5, dim=(modeltype!="SI"))
    # lambda_region = rep(log(1), length(levels(df_animals$region))),
    # lambda_species = rep(log(1), length(levels(df_animals$species))),
    # gamma_species = matrix(-1, nrow=(modeltype!="SI"), ncol=length(levels(df_animals$species))),
    # lambda_sample_type = rep(log(1), length(levels(df_animals$sample_type))),
    # lambda_outdoor = array( rep(0, length(levels(df_animals$outdoor)))),
    # lambda_population = array( rep(0, length(levels(df_animals$population_id)))),
    # lambda_population_raw = array( rep(1/length(levels(df_animals$population_id)), length(levels(df_animals$population_id))-1))
    
    # lambda_region_raw = rep(0, length(levels(df_animals$region))-1),
    # lambda_species_raw = rep(0, length(levels(df_animals$species))-1),
    # gamma_species_raw = matrix(0, nrow=(modeltype!="SI"), ncol=length(levels(df_animals$species))-1),
    # lambda_sample_type_raw = array(rep(0, length(levels(df_animals$sample_type))-1)),
    # lambda_outdoor_raw = array( rep(0, length(levels(df_animals$outdoor))-1)),
    
    # age_raw = runif( n=nrow(df_animals), min=0.1, max=0.9),
    # phi = runif( 1, 0.5, 1.5 ) 
  )
}

inits_sis = function() {
  return(list(
    sigma_lambda_region = 0.25,
    sigma_lambda_species = 0.5,
    
    sigma_lambda_sample_type = 1,
    
    lambda_baseline = -2,
    mu_gamma_species = array( -2, dim=(modeltype!="SI") ),
    sigma_gamma_species = array( 4, dim=(modeltype!="SI")),
    
    # lambda_region_raw = rep(-1, length(levels(df_animals$region))-1),
    # lambda_species_raw = rep(1, length(levels(df_animals$species))-1),
    # gamma_species_raw = matrix(1, nrow=(modeltype!="SI"), ncol=length(levels(df_animals$species))-1),
    # lambda_sample_type_raw = array(rep(-1, length(levels(df_animals$sample_type))-1)),
    # lambda_outdoor_raw = array( rep(0, length(levels(df_animals$outdoor))-1)),
    # 
    age_raw = rep(0.5, nrow(df_animals)),
    phi = 0.8 ))
  
}
