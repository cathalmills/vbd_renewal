#Functions Script for manuscript:
# "Renewal equations for vector-borne diseases"

#Palettes for plotting ----
blue_palette_setup_function <- function(){
  my_palette <- RColorBrewer::brewer.pal(n = 8, name = "BuPu")
  blue_pal <- my_palette[c(4, 6, 8)]
  return(blue_pal)
}
green_palette_setup_function <- function(){
  my_palette <- RColorBrewer::brewer.pal(n = 8, name = "Greens")
  green_pal <- my_palette[c(4, 6, 8)]
  return(green_pal)
}


#Analytical investigation----
#Relationship between R0 and r under assumption of exponentially distributed
# birth processes b_EM(t,a_E) and b_ME(t, a_M)

#Figure 2 (A)-(B)
#Function that allows derivation of r (numerically) and R (phi(0))----
phi_of_gamma <- function(gamma,
                         lambda_e = lambda_e,
                         lambda_m = lambda_m,
                         mu = mu){
  #Let t->infty
  res <- 1/(lambda_e + gamma) *(1/(lambda_m + mu + gamma))
  return(res)
}


R_r_relationship_exponential_diff_lambdas <- function(){
  #Function for deriving r under different assumed regimes
  
  #Can generalise this function in future works
  
  
  tmp_lambdas_dt <- data.table(lambda_e = c(0.5, 0.75, 1.0),
                               lambda_m = c(0.5, 0.75, 1.0),
                               mu = c(0.025, 0.125))
  lambdas_dt <- data.table(expand.grid(tmp_lambdas_dt))
  #Just going to look at three cases for simplicity 
  lambdas_dt <- lambdas_dt[which(lambda_e == lambda_m),]
  lambdas_dt[, IND:= seq(1, nrow(lambdas_dt))]
  number_lambdas <- nrow(lambdas_dt)
  
  ovr_gamma_dt <- NULL
  for(i in 1:number_lambdas){
    tmp_lambdas_dt <- subset(lambdas_dt, IND == i)
    tmp_lambda_e <- tmp_lambdas_dt$lambda_e
    tmp_lambda_m <- tmp_lambdas_dt$lambda_m
    tmp_mu <- tmp_lambdas_dt$mu

    #Numerically solve for r as unique root of phi(gamma) = 1
    tmp_r_res <- uniroot(function(gamma) phi_of_gamma(gamma,
                                                      lambda_e = tmp_lambda_e,
                                                      lambda_m = tmp_lambda_m,
                                                      mu = tmp_mu) - 1, interval = c(-1.2, 4), tol = 1e-16)
    gamma_input_seq <- seq(-0.25, 1.5, by = 0.01)
    phi_of_gamma_values <- sapply(gamma_input_seq, function(x) phi_of_gamma(x, lambda_e = tmp_lambda_e,
                                                                            lambda_m = tmp_lambda_m,
                                                                            mu = tmp_mu))
    #R:= phi(0)
    #Here Rt_val is a value for R_0 as we assume time-independent birth processes
    Rt_val <- phi_of_gamma(0,
                           lambda_e = tmp_lambda_e,
                           lambda_m = tmp_lambda_m,
                           mu = tmp_mu)
    
    #data.table to store value of gamma, corresponding value of phi(gamma),
        # eventual values of R_0 and r,
        # and
        # the values of the birth+death processes used
    gamma_dt <- data.table(gamma = gamma_input_seq,
                           phi = phi_of_gamma_values,
                           R = rep(Rt_val, length(gamma_input_seq)),
                           r = rep(c(tmp_r_res$root), length(gamma_input_seq)),
                           lambda_e = rep(tmp_lambda_e, length(gamma_input_seq)),
                           lambda_m = rep(tmp_lambda_m, length(gamma_input_seq)),
                           mu = rep(tmp_mu, length(gamma_input_seq)),
                           IND = rep(i, length(gamma_input_seq)))
    
    
    
    ovr_gamma_dt <- rbind(ovr_gamma_dt, gamma_dt)
  }
  return(ovr_gamma_dt)
}

process_and_plot_R_r_relationship <- function(ovr_resulting_gamma_dt,
                                              green_pal){
  #Function for processing R-r results and plotting
  #Future generalisations are possible here
  
  
  #Set up data.table for drawing horizontal + vertical 
    # lines that intersect phi(gamma) to derive r as unique
      # root of phi(gamma) = 1
  intersection_dt <- unique(subset(ovr_resulting_gamma_dt, 
                       select = c("r", "IND")))
  setkeyv(intersection_dt, "IND")
  ovr_intersection_dt <- NULL
  for(i in 1:nrow(intersection_dt)){
    tmp_intersection_dt <- data.table(x = c(-Inf, intersection_dt$r[i], intersection_dt$r[i]),
                       y = c(1, 1, 0))
    tmp_intersection_dt[, IND:= factor(i)] #For colour scheme
    ovr_intersection_dt <- rbind(ovr_intersection_dt, tmp_intersection_dt)
  }
  ovr_intersection_dt[, IND:= factor(IND)]
  #Figure 1 (A)
  first_ovr_resulting_gamma_dt <- subset(ovr_resulting_gamma_dt, IND %in% c(1, 2, 3))
  first_ovr_tmp <- subset(ovr_intersection_dt, IND %in% c(1, 2, 3))
  
  #Just extract r values
  rt_seq <- unique(first_ovr_tmp[which(!is.infinite(x))]$x)
  rt_1 <- rt_seq[1]
  rt_2 <- rt_seq[2]
  rt_3 <- rt_seq[3]
  
  #Just extract R_0 values
  Rt_dt <- unique(subset(first_ovr_resulting_gamma_dt, gamma == 0 & IND %in% c(1, 2, 3)))
  Rt_seq <- Rt_dt$phi #R_0
  first_ovr_resulting_gamma_dt[, IND:= factor(IND, levels = c(1, 2, 3),labels = c(1, 2, 3))]
  first_ovr_tmp[, IND:= factor(IND, levels = c(1, 2, 3),labels = c(1, 2, 3))]
  Rt_dt[, IND:= factor(IND, levels = c(1, 2, 3),labels = c(1, 2, 3))]
  
  labels <- c("1" = expression(lambda["E"] ~ "= 0.5" * ", "* lambda["M"] ~ "= 0.5"),
              "2" = expression(lambda["E"] ~ "= 0.75" * ", "* lambda["M"] ~ "= 0.75"),
              "3" = expression(lambda["E"] ~ "= 1.0" * ", "* lambda["M"] ~ "= 1.0"))
  # first_ovr_resulting_gamma_dt[, IND:= factor(IND)]
  first_R_r_plot <- ggplot(first_ovr_resulting_gamma_dt)+
    geom_line(aes(x = gamma, y = phi, color = IND),
              linewidth = 1.4)+
    theme_bw()+
    geom_line(data = first_ovr_tmp, 
              aes(x = x, y = y, color = IND),
              linetype = 2,
              linewidth = 1.4,
              alpha = 0.4)+
    geom_point(data= Rt_dt, aes(x = gamma, y = phi, color = IND),
               size = 5)+
    geom_point(data= Rt_dt, aes(x = r, y = gamma, color = IND),
               size = 5)+
    #Add r and R labels
    annotate("text", x = rt_1, y= 0.025, size = 12, 
             label = expression(r) ,
             color = green_pal[1], fontface = 3, vjust = -0.12)+
    annotate("text", x = rt_2, y= 0.025, size = 12, 
             label = expression(r) ,
             color = green_pal[2], fontface = 3, vjust = -0.12)+
    annotate("text", x = rt_3, y= 0.025, size = 12, 
             label = expression(r) ,
             color = green_pal[3], fontface = 3, vjust = -0.12)+
    annotate("text", y = Rt_seq[1]+ 0.005, x= 0.025, size = 12, 
             label = expression(R[0]) ,
             color = green_pal[1], fontface = 3, vjust = -0.1)+
    annotate("text", y = Rt_seq[2] + 0.005, x= 0.025, size = 12, 
             label = expression(R[0]) ,
             color = green_pal[2], fontface = 3, vjust = -0.1)+
    annotate("text", y = Rt_seq[3]+ 0.005, x= 0.025, size = 12, 
             label = expression(R[0]) ,
             color = green_pal[3], fontface = 3, vjust = -0.1)+
    labs(x = expression(gamma),
         y = expression(phi ~ "(" * gamma ~ ")"))+
    scale_color_manual(values = c("1" = green_pal[1],
                                  "2" = green_pal[2],
                                  "3" = green_pal[3]),
                       labels = labels)+
    geom_hline(yintercept = 0, linetype = "solid", color = "black") +  # Add horizontal line at y=0
    geom_vline(xintercept = 0, linetype = "solid", color = "black") +  # Add vertical line at x=0+
    scale_y_continuous(limits = c(0,5.1), expand = c(0, 0),
                       breaks = seq(0, 5.0, by = 1.0),
                       labels = label_number(accuracy = 0.1)) +  
    annotate("text", x = 1.25, 4.0, size = 11, label = expression(mu ~ "= 0.025"),
             color = "darkblue", fontface = 3)+
    
    theme(text = element_text(size = 28),
          axis.text.x = element_text(size=28),
          axis.text.y = element_text(size=28),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.title = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_blank(),
          axis.title=element_text(size=28), 
          legend.text=element_text(size=28)+geom_text(size = 28),
          legend.position = "bottom",legend.key.size = unit(3.5, 'cm'))+
    guides(color = guide_legend(nrow =1))
  
  
  #Repeat for higher death rate (mu = 0.125)
  second_ovr_resulting_gamma_dt <- subset(ovr_resulting_gamma_dt, IND %in% c(7, 8, 9))
  second_ovr_tmp <- subset(ovr_intersection_dt, IND %in% c(7, 8, 9))
  second_ovr_resulting_gamma_dt[, IND:= factor(IND, levels = c(7, 8, 9),labels = c(7, 8, 9))]
  second_ovr_tmp[, IND:= factor(IND, levels = c(7, 8, 9),labels = c(7, 8, 9))]
  
  #Just extract r values
  rt_seq <- unique(second_ovr_tmp[which(!is.infinite(x))]$x)
  rt_1 <- rt_seq[1]
  rt_2 <- rt_seq[2]
  rt_3 <- rt_seq[3]
  
  #Just extract R_0 values
  Rt_dt <- unique(subset(second_ovr_resulting_gamma_dt, gamma == 0 & IND %in% c(7, 8, 9)))
  Rt_dt[, IND:= factor(IND, levels = c(7, 8, 9),labels = c(7, 8, 9))]
  Rt_seq <- Rt_dt$phi #R_0
  # second_ovr_resulting_gamma_dt[, IND:= factor(IND)]
  
  second_R_r_plot <- ggplot(second_ovr_resulting_gamma_dt)+
    geom_line(aes(x = gamma, y = phi, color = IND),
              linewidth = 1.4)+
    theme_bw()+
    geom_line(data = second_ovr_tmp, 
              aes(x = x, y = y, color = IND),
              linetype = 2,
              linewidth = 1.4,
              alpha = 0.4)+
    geom_point(data= Rt_dt, aes(x = gamma, y = phi, color = IND),
               size = 5)+
    geom_point(data= Rt_dt, aes(x = r, y = gamma, color = IND),
               size = 5)+
    #Add r and R labels
    annotate("text", x = rt_1, y= 0.025, size = 12, 
             label = expression(r) ,
             color = green_pal[1], fontface = 3, vjust = -0.12)+
    annotate("text", x = rt_2, y= 0.025, size = 12, 
             label = expression(r) ,
             color = green_pal[2], fontface = 3, vjust = -0.12)+
    annotate("text", x = rt_3, y= 0.025, size = 12, 
             label = expression(r) ,
             color = green_pal[3], fontface = 3, vjust = -0.12)+
    annotate("text", y = Rt_seq[1]+ 0.005, x= 0.025, size = 12, 
             label = expression(R[0]) ,
             color = green_pal[1], fontface = 3, vjust = -0.1)+
    annotate("text", y = Rt_seq[2] + 0.005, x= 0.025, size = 12, 
             label = expression(R[0]) ,
             color = green_pal[2], fontface = 3, vjust = -0.1)+
    annotate("text", y = Rt_seq[3]+ 0.005, x= 0.025, size = 12, 
             label = expression(R[0]) ,
             color = green_pal[3], fontface = 3, vjust = -0.1)+
    labs(x = expression(gamma),
         y = expression(phi ~ "(" * gamma ~ ")"))+
    scale_color_manual(values = c("7" = green_pal[1],
                                  "8" = green_pal[2],
                                  "9" = green_pal[3]),
                       labels = labels)+
    geom_hline(yintercept = 0, linetype = "solid", color = "black") +  # Add horizontal line at y=0
    geom_vline(xintercept = 0, linetype = "solid", color = "black") +  # Add vertical line at x=0+
    scale_y_continuous(limits = c(0,5.1), expand = c(0, 0),
                       breaks = seq(0, 5.0, by = 1.0),
                       labels = label_number(accuracy = 0.1)) +  
    annotate("text", x = 1.25, 4.0, size = 11, label = expression(mu ~ "= 0.125"),
             color = "darkblue", fontface = 3)+
    theme(text = element_text(size = 28),
          axis.text.x = element_text(size=28),
          axis.text.y = element_text(size=28),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.title = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_blank(),
          axis.title=element_text(size=28), 
          legend.text=element_text(size=28)+geom_text(size = 28),
          legend.position = "bottom",legend.key.size = unit(3.5, 'cm'))+
    guides(color = guide_legend(nrow =1))
  
  

  top_panel_plot <- ggarrange(plotlist = list(NULL, NULL,
                                              first_R_r_plot,
                                              second_R_r_plot), 
                              common.legend = TRUE,
                              legend = "bottom",
                              labels = c("(A)", "(B)", "", ""),
                              heights = c(0.1, 0.95),
                              font.label = list(size = 34, face = "bold"),
                              nrow = 2, ncol = 2)
  return(top_panel_plot)
}




#R0 as a function of r and the birth+death process parameters
#Appendix B.4 Derivation
#Figure 2 (C)-(D) ----
R_r_analytical_function <- function(mu,
                                    lambda_e,
                                    lambda_m,
                                    r){
  res <- (lambda_e + r)*(lambda_m + mu + r)/(lambda_e*(lambda_m + mu))
  return(res)
}


#Note for later: These would probably be better separated into two functions
#Function for calculating + plotting R as a function of r and the birth+death process parameters
#Mosquito death rate is varied and infection birth processes fixed
R_r_mu_varying_plot_function <- function(blue_pal){
  mu_seq <- c(0.025, 0.125, 0.5) #Vary death rate
  fixed_lambda_e <- 0.75
  fixed_lambda_m <- 0.75
  r_seq <- seq(-0.25, 0.5, by = 0.01)
  plot_mus_R_r_dt <- NULL #data.table to store results
  for(i in 1:length(mu_seq)){
    tmp_mu <- mu_seq[i]
    R_res <- sapply(r_seq, function(x){
      R_r_analytical_function(tmp_mu, fixed_lambda_e, fixed_lambda_m, x)})
    R_r_dt <- data.table(R_val = R_res,
                         r_val = r_seq)
    R_r_dt[, mu := tmp_mu]
    R_r_dt[, lambda_m := fixed_lambda_m]
    R_r_dt[, lambda_e := fixed_lambda_e]
    plot_mus_R_r_dt <- rbind(plot_mus_R_r_dt, R_r_dt)
  }
  
  mu_labels <- c("1" = expression(mu ~ "= 0.025  "),
                 "2" = expression(mu ~ "= 0.125  "),
                 "3" = expression(mu ~ "= 0.50  "))
  #Index as factor variable for colour scheme
  plot_mus_R_r_dt[, IND:= factor(mu, levels = c(0.025, 0.125, 0.5),
                                 labels = c(1, 2, 3))]
  mus_R_r_plot <- 
    ggplot(plot_mus_R_r_dt, 
           aes(x = R_val, y = r_val, color = IND)) + geom_line(linewidth = 1.5)+
    theme_bw()+
    coord_cartesian(ylim = c(-0.3, 0.55),
                    xlim = c(0, 2.8),
                    expand = F)+  
    scale_x_continuous(limits = c(0, 2.8),
                       breaks = seq(0, 2.8, by = 0.5))+
    labs(y = paste0(expression(r), " (per day)"),
         x = expression(R["0"]))+
    scale_color_manual(values = c("1" = blue_pal[1],
                                  "2" = blue_pal[2],
                                  "3" = blue_pal[3]),
                       labels = mu_labels)+
    geom_hline(yintercept = 0, linetype = "solid", color = "black") +  # Add horizontal line at y=0
    geom_vline(xintercept = 1, linetype = "solid", color = "black") +  # Add vertical line at x=0+
    annotate("text", x = 2.25, 0.15, size = 11, label = expression(lambda["M"] ~ "= 0.75, "),
             color = "darkblue", fontface = 3)+
    annotate("text", x = 2.25, 0.075, size = 11, label = expression(lambda["E"] ~ "= 0.75 "),
             color = "darkblue", fontface = 3)+
    theme(text = element_text(size = 28),
          axis.text.x = element_text(size=28),
          axis.text.y = element_text(size=28),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.title = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_blank(),
          axis.title=element_text(size=28), 
          legend.text=element_text(size=28)+geom_text(size = 28),
          legend.position = "bottom",legend.key.size = unit(2.0, 'cm'))+
    guides(color = guide_legend(nrow =1))
  return(mus_R_r_plot)
}

#Figure 2 (D)
#As above, but vary infection process rate and fix death rate ----
R_r_lambda_m_varying_plot_function <- function(blue_pal){
  lambda_m_seq <- c(0.5, 0.75, 1.0)
  fixed_lambda_e <- 0.75
  fixed_mu <- 0.075
  r_seq <- seq(-0.25, 0.5, by = 0.01)
  plot_lambda_m_R_r_dt <- NULL
  for(i in 1:length(lambda_m_seq)){
    tmp_lambda_m <- lambda_m_seq[i]
    R_res <- sapply(r_seq, function(x){
      R_r_analytical_function(fixed_mu, fixed_lambda_e, tmp_lambda_m, x)})
    R_r_dt <- data.table(R_val = R_res,
                         r_val = r_seq)
    R_r_dt[, mu := fixed_mu]
    R_r_dt[, lambda_m := tmp_lambda_m]
    R_r_dt[, lambda_e := fixed_lambda_e]
    plot_lambda_m_R_r_dt <- rbind(plot_lambda_m_R_r_dt, R_r_dt)
  }

  plot_lambda_m_R_r_dt[, lambda_m:= factor(lambda_m)]
  plot_lambda_m_R_r_dt[, log_R_val:= log(R_val)]
  
  lambda_m_labels <- c("1" = expression(lambda["M"] ~ "= 0.5  "),
                       "2" = expression(lambda["M"] ~ "= 0.75  "),
                       "3" = expression(lambda["M"] ~ "= 1.0  "))
  plot_lambda_m_R_r_dt[, IND:= factor(lambda_m, levels = c(0.5, 0.75, 1.0),
                                      labels = c(1, 2, 3))]
  R_r_relationship_lambda_m_plot <- 
    ggplot(plot_lambda_m_R_r_dt, 
           aes(x = R_val, y = r_val, color = IND)) + geom_line(linewidth = 1.5)+
    theme_bw()+
    coord_cartesian(ylim = c(-0.3, 0.55),
                    xlim = c(0, 2.8),
                    expand = F)+  
    scale_x_continuous(limits = c(0, 2.8),
                       breaks = seq(0, 2.8, by = 0.5))+
    labs(y = paste0(expression(r), " (per day)"),
         x = expression(R["0"]))+
    scale_color_manual(values = c("1" = blue_pal[1],
                                  "2" = blue_pal[2],
                                  "3" = blue_pal[3]),
                       labels = lambda_m_labels)+
    geom_hline(yintercept = 0, linetype = "solid", color = "black") +  # Add horizontal line at y=0
    geom_vline(xintercept = 1, linetype = "solid", color = "black") +  # Add vertical line at x=0+
    annotate("text", x = 2.2, 0.1, size = 11, label = expression(mu ~ "= 0.075"),
             color = "darkblue", fontface = 3)+
    theme(text = element_text(size = 28),
          axis.text.x = element_text(size=28),
          axis.text.y = element_text(size=28),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.title = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_blank(),
          axis.title=element_text(size=28), 
          legend.text=element_text(size=28)+geom_text(size = 28),
          legend.position = "bottom",legend.key.size = unit(2.0, 'cm'))+
    guides(color = guide_legend(nrow =1))
  return(R_r_relationship_lambda_m_plot)
}


plot_four_panel_figure_two <- function(top_panel_plot,
                                       mus_R_r_plot,
                                       R_r_relationship_lambda_m_plot){
  lower_panel_plot <- ggarrange(mus_R_r_plot,
                              R_r_relationship_lambda_m_plot, nrow = 1,
                              labels = c("(C)", "(D)"),
                              vjust = -0.2,
                              font.label = list(size = 34, face = "bold"))
  R_r_analytical_plot_two_panels <- ggarrange(top_panel_plot,
                                            lower_panel_plot, 
                                            nrow = 2)
  return(R_r_analytical_plot_two_panels) 
}




# Siraj et Al (2017) ----

#Following functions are from Siraj et Al (2017) 
# (just needed for toy temperature-dependent example)


#Temp-dependent death rate of mosquitoes
mortalityRT <- function(temp,fldcxn=fieldcorxn) {
  dd<-seq(0,120,length.out=(120*24+2))
  nwdd<-  data.frame(Days=dd,Temperature=rep(temp, (120*24+2)), Study_number=5, Feed_B=2, Feed_S=1) ## uses algam_85re.Rdata
  nwdd<-  cbind(nwdd,logDay=log(nwdd$Days+1), logTemp=log(nwdd$Temperature+1))  ## +1 avoids log(0) 
  prediction <-as.vector(unlist(stats::predict(algam,newdata = (nwdd), se.fit = TRUE, type = "response")$fit))
  prediction<- prediction[-1]
  prediction[1:24]<- prediction[1:24]/prediction[1]
  prediction[which(prediction>1)]<-1
  prediction[which(prediction<=0.001)]<-0
  diffDeath<- -diff(prediction)
  diffDeath<- diffDeath/sum(diffDeath)
  return(1/sum(dd[2:length(prediction)]*diffDeath) + fldcxn) 
}

#Intrinsic incubation period
iip_func <- function(iip_params = iip_params, 
                     age,
                     temp){
  iip_mu <- exp(exp(iip_params[1] + iip_params[2]*temp))
  iip_tau <- iip_params[3]
  tmp_iipdf = dlnorm(age, meanlog=log(iip_mu), sdlog=(1/sqrt(iip_tau)))
  if(age == 0){tmp_iipdf <- 0}
  return(tmp_iipdf)
}

#Human-to-mosquito period
hmtp_func <- function(age){
  barer = hinfectiousness(1)$histo ## histogram in Nishuara & Halstead (we are using 0:6, actually is -2:4 relative)
  stmle = mle.norm(barer)
  stpdf = dnorm(age,stmle[1],stmle[2])
  if(age == 0){stpdf <- 0}
  return(stpdf)
}

#Extrinsic incubation period
eip_func <- function(eip_params,
                     age,
                     temp){
  
  if(age != 0)
  {
    eip_mu <- exp(exp(eip_params[1] + eip_params[2]*temp))
    eip_tau <- eip_params[3]
    tmp_eipdf = dlnorm(age, meanlog=log(eip_mu), sdlog=(1/sqrt(eip_tau)))
  }
  else{
    tmp_eipdf <- 0
  }
  return(tmp_eipdf)
}

#Mosquito-to-human period
mhtp_function <- function(age,
                          temp){
  tmp_mortp <- mortalityRT(temp, fldcxn = fieldcorxn)
  mortp <- tmp_mortp
  tmp_mortdf = dexp(age, mortp)
  
  if(age == 0){tmp_mortdf <- 0}
  return(tmp_mortdf)
}

#Used for toy example of temperature-dependent GI
set_siraj_2017_params <- function(){
  siraj_2017_params <- c(0.56, 0, 13.7, 
              2.9, -0.08, 4.9) #Default Parameters from Siraj et al(2017)
  return(siraj_2017_params)
}



#Random number generator functions for distributions from
  # Siraj et al (2017)

rand_iip_func <- function(iip_params = iip_params, 
                          temp){
  iip_mu <- exp(exp(iip_params[1] + iip_params[2]*temp))
  iip_tau <- iip_params[3]
  sampled_iip = rlnorm(1, meanlog=log(iip_mu), sdlog=(1/sqrt(iip_tau)))
  return(sampled_iip)
}

rand_hmtp_func <- function(){
  barer = hinfectiousness(1)$histo ## histogram in Nishuara & Halstead (we are using 0:6, actually is -2:4 relative)
  stmle = mle.norm(barer)
  sampled_hmtp = rnorm(1,stmle[1],stmle[2])
  return(sampled_hmtp)
}

rand_eip_func <- function(eip_params,
                          temp){
  eip_mu <- exp(exp(eip_params[1] + eip_params[2]*temp))
  eip_tau <- eip_params[3]
  sampled_eip = rlnorm(1, meanlog=log(eip_mu), sdlog=(1/sqrt(eip_tau)))
  return(sampled_eip)
}

rand_mhtp_function <- function(temp){
  tmp_mortp <- mortalityRT(temp, fldcxn = fieldcorxn)
  sampled_mhtp = rexp(1, tmp_mortp)
  return(sampled_mhtp)
}

setup_simulated_data <- function(){
  sim_step_size <- 1
  sim_seq_times <- seq(1, 205, by = sim_step_size)
  sim_dt <- data.table(TIME = seq(1, max(sim_seq_times), 
                                  by = sim_step_size))
  fitting_sim_dt <- subset(sim_dt, TIME >= 150)
  fitting_sim_dt[, FITTING_TIME:= seq(1, nrow(fitting_sim_dt))]
  
  sim_temp <- rep(0, nrow(sim_dt))
  sim_temp[1] <- runif(1, 20, 45)
  for(i in 2:length(sim_temp)){
    tmp_noise <- rnorm(1, 0, 2)
    sim_temp[i] <- sim_temp[i-1] + tmp_noise
    
    #Reflecting boundaries
    if(sim_temp[i] < 20){
      sim_temp[i] <- sim_temp[i-1] + abs(tmp_noise)
    }
    else if(sim_temp[i] > 45){
      sim_temp[i] <- sim_temp[i-1] - (tmp_noise)
      
    }
    # }
  }
  sim_dt[, temp:= sim_temp]
  return(sim_dt)
}


read_all_times_data <- function(option){
  #Fortaleza (56-day)
  if(option == "fortaleza_season"){
    location_all_data <- data.table(read.csv("fortaleza_season_temps_dt.csv"))
  }
  #Fortaleza (Year)
  else if (option == "fortaleza_year"){
    location_all_data <- data.table(read.csv("fortaleza_2012_temps_dt.csv"))
    
  }  
  #Singapore
  else if (option == "singapore_season"){
    location_all_data <- data.table(read.csv("singapore_season_temps_dt.csv"))
  }
  else if (option == "foz_year"){
    location_all_data <- data.table(read.csv("foz_2012_temps_dt.csv"))
    
  }
  else if (option == "simulated"){
  #Simulated
    location_all_data <- data.table(read.csv("simulated_dt.csv"))
  }
}

read_fitting_period_data <- function(location_all_data, option){
  if(option == "fortaleza_season"){
    #Fortaleza (57 days)
    #Chosen 57 days to be length of short period with cases
    fitting_dt <- subset(location_all_data, WEEK >= 13 & WEEK <= 20 & YEAR == 2012) #Subset portion of dengue season
    fitting_dt[, FITTING_TIME:= seq(1, nrow(fitting_dt))]
  }
  else if (option == "fortaleza_year"){
    #Fortaleza (Year)
    fitting_dt <- subset(location_all_data, YEAR == 2012)
    fitting_dt[, FITTING_TIME:= seq(1, nrow(fitting_dt))]
  }  
  else if (option == "singapore_season"){
  #Singapore
  #Chosen 57 days to be length of short period with cases
    fitting_dt <- subset(location_all_data, TIME > (nrow(location_all_data) - 57))
    fitting_dt[, FITTING_TIME:= seq(1, nrow(fitting_dt))]
  }
  else if (option == "foz_year"){
    fitting_dt <- subset(location_all_data, YEAR == 2012)
    fitting_dt[, FITTING_TIME:= seq(1, nrow(fitting_dt))]
  }
  else if (option == "simulated"){
  #Simulated
    fitting_dt <- subset(location_all_data, TIME >= 150)
    fitting_dt[, FITTING_TIME:= seq(1, nrow(fitting_dt))]
  }

  return(fitting_dt)
}

#Estimating generation time distribution using a 
  #Monte Carlo sampler accounting for times at which stages of the transmission cycle occur
  # i.e. our "stage-specific" approach
time_varying_monte_carlo_generator_function_temps <- function(number_samples,
                                                              fixed_t,
                                                              temps_dt,
                                                              params){
  #Inputs: number of MC samples, calendar time t, temperature data.table,
    # parameters for assumed distributions between stages
  iip_params <- params[c(1:3)]
  eip_params <- params[c(4:6)]
  tau <- rep(NA, number_samples)
  
  #These are aV values at time t
  aV_samples <- replicate(number_samples, 
                          rand_mhtp_function(temps_dt[which(TIME == fixed_t)]$temp))
  #Empty vectors to store samples + densities
  aW_samples <- rep(NA, number_samples)
  aI_samples <- rep(NA, number_samples)
  aE_samples <- rep(NA, number_samples)
  
  aV_densities <- rep(NA, number_samples)
  aW_densities <- rep(NA, number_samples)
  aI_densities <- rep(NA, number_samples)
  aE_densities <- rep(NA, number_samples)
  unnormalised_gt_density <- rep(NA, number_samples)
  for(i in 1:number_samples){
    if(i %% 2500 ==0){gc()} #Memory!
    tmp_aV <- aV_samples[i] #This age will be used downstream
    aV_densities[i] <- mhtp_function(tmp_aV,
                                     temps_dt[which(TIME == fixed_t)]$temp)
    aW_samples[i] <- rand_eip_func(eip_params, 
                                   temps_dt[which(TIME == fixed_t-round(tmp_aV))]$temp)
    if(aW_samples[i] > (fixed_t - tmp_aV)){ #recall, we assume t>>a_W
      aW_samples[i] <- rand_eip_func(eip_params, 
                                     temps_dt[which(TIME == fixed_t-round(tmp_aV))]$temp)
    }
    tmp_aW <-  aW_samples[i]
    aW_densities[i] <- eip_func(eip_params, 
                                age = tmp_aW,
                                temps_dt[which(TIME == round(fixed_t-tmp_aV))]$temp)
    
    aI_samples[i] <- rand_hmtp_func()
    if(aI_samples[i] > (fixed_t - tmp_aV - tmp_aW)){
      aI_samples[i] <- rand_hmtp_func()
    }
    tmp_aI <- c(aI_samples[i])
    aI_densities[i] <- hmtp_func(age = tmp_aI)

    aE_samples[i] <-
      rand_iip_func(iip_params = iip_params,
                    temps_dt[which(TIME == max(round(fixed_t - tmp_aV - tmp_aW - tmp_aI), 1))]$temp)
    if(aE_samples[i] > (fixed_t - tmp_aV - tmp_aW - tmp_aI)){
      aE_samples[i]<-  rand_iip_func(iip_params = iip_params, 
                                     temps_dt[which(TIME == max(round(fixed_t - tmp_aV - tmp_aW - tmp_aI), 1))]$temp)
      
    }
    
    tmp_aE <- c(aE_samples[i])
    
    
    #Total delay tau = sample from w(t, tau)
    tau[i] <- tmp_aE + tmp_aW + tmp_aV + tmp_aI
    
    aE_densities[i] <- iip_func(iip_params = iip_params,
                                age = tmp_aE,
                                temp = temps_dt[which(TIME == max(round(fixed_t - tmp_aV - tmp_aW - tmp_aI), 1))]$temp)
    unnormalised_gt_density[i] <- aV_densities[i]*aW_densities[i]*aI_densities[i]*aE_densities[i]
  }
  
  gt_densities <-unnormalised_gt_density/sum(unnormalised_gt_density)
  gt_dt <- data.table(density = gt_densities,
                      tau = tau)
  
  return(list(tau, gt_densities,
              gt_dt, 
              aE_samples, aI_samples, aW_samples, aV_samples))
}


#Ignoring calendar times at which stages of transmission cycle occur
  # GT distribution based on temperature observed only at time t
time_varying_monte_carlo_generator_function_ignoring_transm_cycle_temps <- 
  function(number_samples,
           fixed_t,
           temps_dt,
           params){
  iip_params <- params[c(1:3)]
  eip_params <- params[c(4:6)]
  tau <- rep(NA, number_samples)
  
  aV_samples <- replicate(number_samples, 
                          rand_mhtp_function(temps_dt[which(TIME == fixed_t)]$temp))
  
  aW_samples <- rep(NA, number_samples)
  aI_samples <- rep(NA, number_samples)
  aE_samples <- rep(NA, number_samples)
  aV_densities <- rep(NA, number_samples)
  aW_densities <- rep(NA, number_samples)
  aI_densities <- rep(NA, number_samples)
  aE_densities <- rep(NA, number_samples)
  unnormalised_gt_density <- rep(NA, number_samples)
  for(i in 1:number_samples){
    if(i %% 2500 ==0){gc()}
    tmp_aV <- aV_samples[i]
    aV_densities[i] <- mhtp_function(tmp_aV,
                                     temps_dt[which(TIME == fixed_t)]$temp)
    
    aW_samples[i] <- rand_eip_func(eip_params, 
                                   temps_dt[which(TIME == fixed_t)]$temp)
    if(aW_samples[i] > (fixed_t - tmp_aV)){
      aW_samples[i] <- rand_eip_func(eip_params, 
                                     temps_dt[which(TIME == fixed_t)]$temp)
    }
    tmp_aW <-  aW_samples[i]
    aW_densities[i] <- eip_func(eip_params, 
                                age = tmp_aW,
                                temps_dt[which(TIME == round(fixed_t))]$temp)
    
    aI_samples[i] <- rand_hmtp_func()
    if(aI_samples[i] > (fixed_t - tmp_aV - tmp_aW)){
      aI_samples[i] <- rand_hmtp_func()
    }
    tmp_aI <- c(aI_samples[i])
    aI_densities[i] <- hmtp_func(age = tmp_aI)
    
    aE_samples[i] <- 
      rand_iip_func(iip_params = iip_params, 
                    temps_dt[which(time == round(fixed_t))]$temp)
    if(aE_samples[i] > (fixed_t - tmp_aV - tmp_aW - tmp_aI)){
      aE_samples[i] <- rand_iip_func(iip_params = iip_params, 
                                     temps_dt[which(time == round(fixed_t))]$temp)
    }
    
    tmp_aE <- c(aE_samples[i])
    
    tau[i] <- tmp_aE + tmp_aW + tmp_aV + tmp_aI
    
    aE_densities[i] <- iip_func(iip_params = iip_params,
                                age = tmp_aE,
                                temp = temps_dt[which(time == round(fixed_t))]$temp)
    unnormalised_gt_density[i] <- aV_densities[i]*aW_densities[i]*aI_densities[i]*aE_densities[i]
  }
  
  gt_densities <-unnormalised_gt_density/sum(unnormalised_gt_density)
  gt_dt <- data.table(density = gt_densities,
                      tau = tau)
  return(list(tau, gt_densities,
              gt_dt,
              aE_samples, aI_samples, aW_samples, aV_samples))
}


#Application to simulated + real-world temperature data ----

#Function for running the estimate of time-varying GI distribution by Monte Carlo sampling
run_time_varying_monte_carlo <- 
  function(number_mc_samples,
           all_data,
           fitting_data,
           transm_cycle = TRUE,
           siraj_2017_params){
    #all_data includes historical temperatures relevant for estimating current GI
    tau_samples_dt <- NULL
    for(i in 1:length(fitting_data$FITTING_TIME)){
      #FITTING_TIME is index from 1 to number of rows of fitting data.table
      
      
      #Loop over calendar time t at which estimate GI
      tmp_fixed_t <- fitting_data$TIME[i] 
      #Subset relevant temperatures      
      all_data[, TIME:= as.numeric(TIME)]
      fitting_data[, TIME:= as.numeric(TIME)]
      tmp_time_varying_temps_dt <- 
        subset(all_data, 
               TIME <= tmp_fixed_t)
      if(transm_cycle == TRUE){
        tmp_results <- 
          time_varying_monte_carlo_generator_function_temps(number_mc_samples,
                                                          tmp_fixed_t,
                                                          tmp_time_varying_temps_dt,
                                                          siraj_2017_params)
      }
      else{
        tmp_results <- 
          time_varying_monte_carlo_generator_function_ignoring_transm_cycle_temps(number_mc_samples,
                                                            tmp_fixed_t,
                                                            tmp_time_varying_temps_dt,
                                                            siraj_2017_params)
      }
      tau_samples <- tmp_results[[1]]
      tau_samples_dt <- rbind(tau_samples_dt, tau_samples)

    }
    tau_dt <- as.data.frame(tau_samples_dt)
    setnames(tau_dt, paste0("TIME_", seq(1, ncol(tau_dt))))
    tau_dt[, SAMPLE:= seq(1, nrow(tau_dt))] #SAMPLE = Index for MC samples
    tau_dt <- melt(tau_dt, id.var = "SAMPLE")
    tau_dt[, variable:= gsub("TIME_", "", variable)] 
    tau_dt[, variable:= as.numeric(variable)]
    setnames(tau_dt, "variable", "TIME") #TIME = Date at which we are evaluating gen time distribution
    return(tau_dt)
}



combine_and_compare_gen_time_distributions <- 
  function(tau_dt_1,
           tau_dt_2){
  tau_dt1 <- subset(tau_dt1, select = c("SAMPLE", "TIME", "value"))
  tau_dt1[, TYPE:= "Stage_Specific"]
  tau_dt2 <- subset(tau_dt2, select = c("SAMPLE", "TIME", "value"))
  tau_dt2[, TYPE:= "Existing"]
  tau_both_dt <- rbind(tau_dt, tau_dt2)
}









#Plotting for Figure 3 ----
plot_yearly_gi_dt <- function(dt){
  summary_dt <- dt[, list(CI_U = quantile(value, probs = 0.95),
                          CI_L = quantile(value, probs = 0.05),
                          MEDIAN = quantile(value, probs = 0.5))
                   , by = c("DATE", "TYPE")]
  summary_plot <- ggplot(summary_dt) + 
    geom_line(aes(x = DATE, y = MEDIAN, colour = TYPE),
              linewidth = 1.5)+
    geom_point(aes(x = DATE, y = MEDIAN, colour = TYPE),
               size = 2.5)+
    
    geom_ribbon(aes(x = DATE, ymin = CI_L, ymax = CI_U, fill = TYPE), alpha = 0.15)+
    theme_bw()+
    coord_cartesian(ylim = c(0, 125))+
    scale_y_continuous(breaks = seq(0, 125, by = 25))+
    labs(x = "Date", y = "GT Delay (Days)")+
    theme(legend.position = "bottom")+
    scale_fill_manual(name = "Method",
                      labels = c("Existing", "Stage-\nSpecific"),
                      values = c("forestgreen", "darkorange1"),
                      aesthetics = c("colour", "fill"))+
    theme(text = element_text(size = 28),
          axis.text.x = element_text(size=28),
          axis.text.y = element_text(size=28),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.title = element_text(face = "bold"),
          plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_blank(),
          axis.title=element_text(size=28), 
          legend.text=element_text(size=28)+geom_text(size = 28),
          legend.position = "bottom",legend.key.size = unit(4.5, 'cm'))
  return(summary_plot)
}

plot_yearly_temps_dt <- function(dt){
  temps_plot <- ggplot(dt) + 
    geom_line(aes(x = DATE, y = temp),
              linewidth = 1.15, color = "maroon4")+
    geom_point(aes(x = DATE, y = temp),
               size = 3.5, color = "maroon4", shape = 18)+
    theme_bw()+
    # scale_x_continuous(breaks = seq(0, 365, by = 30))+
    scale_y_continuous(breaks = seq(5, 35, by = 5))+
    coord_cartesian(ylim = c(5, 35))+
    labs(x = "Date", y = "Temperature (Degrees C)")+
    theme(legend.position = "bottom")+
    theme(text = element_text(size = 28),
          axis.text.x = element_text(size=28),
          axis.text.y = element_text(size=28),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.title = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_blank(),
          axis.title=element_text(size=28), 
          legend.text=element_text(size=28)+geom_text(size = 28),
          legend.position = "bottom",legend.key.size = unit(3.5, 'cm'))
  return(temps_plot)
}

arrange_yearly_gi_plots <- function(plot1, plot2){
  arranged_yearly_gi_plot <- ggarrange(NULL, NULL,
            plot1,
            plot2,
            nrow = 2, ncol = 2,
            heights = c(0.1, 1),
            common.legend = TRUE,
            legend = "bottom")
  return(arranged_yearly_gi_plot)
}

plot_temps_dt <- function(dt){
  temps_plot <- ggplot(dt) + 
    geom_line(aes(x = TIME, y = temp),
              linewidth = 1.15, color = "maroon4")+
    geom_point(aes(x = TIME, y = temp),
               size = 3.5, color = "maroon4", shape = 18)+
    theme_bw()+
    scale_x_continuous(breaks = seq(0, 55, by = 10))+
    scale_y_continuous(breaks = seq(20, 45, by = 5))+
    coord_cartesian(xlim = c(1, 56),
                    ylim = c(20, 40))+
    labs(x = "Time Index (Days)", y = "Temperature (Degrees C)")+
    theme(legend.position = "bottom")+
    theme(text = element_text(size = 28),
          axis.text.x = element_text(size=28),
          axis.text.y = element_text(size=28),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.title = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_blank(),
          axis.title=element_text(size=28), 
          legend.text=element_text(size=28)+geom_text(size = 28),
          legend.position = "bottom",legend.key.size = unit(3.5, 'cm'))
  return(temps_plot)
}

arrange_yearly_temps_plots <- function(plot1, plot2){
  arranged_yearly_temps_plot <- ggarrange(NULL, NULL,  #Leaving some space at top for labels
                                plot1,
                                plot2,
                                nrow = 2, ncol = 2,
                                labels = c("(D)", "(E)", "",""),
                                heights = c(0.175, 1),
                                hjust = -0.75,
                                font.label = list(size = 40, face = "bold"))
  return(arranged_yearly_temps_plot)
}

plot_gi_dt <- function(dt){
  summary_dt <- dt[, list(CI_U = quantile(value, probs = 0.95),
                          CI_L = quantile(value, probs = 0.05),
                          MEDIAN = quantile(value, probs = 0.5))
                   , by = c("TIME", "TYPE")]
  summary_plot <- ggplot(summary_dt) + 
    geom_line(aes(x = TIME, y = MEDIAN, colour = TYPE),
              linewidth = 1.5)+
    geom_point(aes(x = TIME, y = MEDIAN, colour = TYPE),
               size = 2.5)+
    
    geom_ribbon(aes(x = TIME, ymin = CI_L, ymax = CI_U, fill = TYPE), alpha = 0.15)+
    theme_bw()+
    scale_x_continuous(breaks = seq(0, 55, by = 10))+
    coord_cartesian(xlim = c(1, 56),
                    ylim = c(0, 80))+
    scale_y_continuous(breaks = seq(0, 80, by = 20))
    labs(x = "Time Index (Days)", y = "GT Delay (Days)")+
    theme(legend.position = "bottom")+
    scale_fill_manual(name = "Method",
                      labels = c("Existing", "Stage-\nSpecific"),
                      values = c("forestgreen", "darkorange1"),
                      aesthetics = c("colour", "fill"))+
    theme(text = element_text(size = 28),
          axis.text.x = element_text(size=28),
          axis.text.y = element_text(size=28),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.title = element_text(face = "bold"),
          plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_blank(),
          axis.title=element_text(size=28), 
          legend.text=element_text(size=28)+geom_text(size = 28),
          legend.position = "bottom",legend.key.size = unit(4.5, 'cm'))
  return(summary_plot)
}

arrange_temp_plots <- function(plot1, plot2, plot3){
  arranged_plot <- ggarrange(NULL, NULL, NULL, #Leaving some space at top for labels
                             plot1, plot2, plot3, 
                             ncol = 3, nrow = 2,
                             heights = c(0.15, 1),
                             labels = c("(A)", "(B)", "(C)", "", "", ""),
                             hjust = -0.75,
                             font.label = list(size = 40, face = "bold"),
                             common.legend = TRUE,
                             legend = "none")
  return(arranged_plot)
}

arrange_gi_plots <- function(plot1, plot2, plot3){
  summary_gi_plot <- ggarrange(plot1, 
                             plot2,
                             plot3, 
                             ncol = 3,
                             common.legend = TRUE,
                             legend = "none")
  return(summary_gi_plot)
}
setup_temps_dt <- function(fitting_dt){
  temps_dt <- subset(temps_dt, select = c("FITTING_TIME", "temp", "DATE"))
  setnames(temps_dt, "FITTING_TIME", "TIME")
  return(temps_dt)
}



arrange_final_five_panel_figure_3_plot <- 
  function(temps_plots, summary_gi_plot, 
           arranged_yearly_temps_plot, arranged_yearly_gi_plot){
  application_gis_plot <- ggarrange(temps_plots, summary_gi_plot,
                                    nrow = 2, heights = c(1, 1), common.legend = TRUE, legend = "none")
  
  yearly_plots <- ggarrange(arranged_yearly_temps_plots,
                            arranged_yearly_gi_plot, common.legend = TRUE,
                            legend = "bottom",
                            heights = c(0.75, 1.1),
                            nrow = 2)
  yearly_plots
}