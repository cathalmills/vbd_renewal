# Created by use_targets().
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
library(tarchetypes)
usethis::use_git()
tar_option_set(packages = c("readr", "dplyr", "ggplot2",
                            "stats", "pracma", "data.table",
                            "tidyverse", "gridExtra", "ggpubr",
                            "scales", "RColorBrewer"),
               format = "rds")
tar_source()
#Data = fortaleza, sao jao, singapore, fortaleza, simulated,
# and data from Siraj et al. (2017)
load("epidemicGrowth.RData") #Load data

list(
  #Figures 1 + 2
  tar_target(green_pal, green_palette_setup_function()),
  tar_target(ovr_gamma_dt, R_r_relationship_exponential_diff_lambdas()),
  tar_target(top_panel_plot, process_and_plot_R_r_relationship(ovr_gamma_dt, green_pal)),
  tar_target(blue_pal, blue_palette_setup_function()),
  tar_target(mus_R_r_plot, R_r_mu_varying_plot_function(blue_pal)),
  tar_target(R_r_relationship_lambda_m_plot, R_r_lambda_m_varying_plot_function(blue_pal)),
  tar_target(R_r_analytical_plot_two_panels, plot_four_panel_figure_two(top_panel_plot,
                                        mus_R_r_plot,
                                        R_r_relationship_lambda_m_plot)),
  
  
  
  #Application to real-world + simulated data
  tar_target(siraj_2017_params,
             set_siraj_2017_params()),
  #Load in temperature data for entire period
  tar_target(fortaleza_season_temps_dt, read_all_times_data(option = "fortaleza_season")),
  tar_target(fortaleza_2012_temps_dt, read_all_times_data(option = "fortaleza_year")),
  tar_target(singapore_season_temps_dt, read_all_times_data(option = "singapore_season")),
  tar_target(foz_2012_temps_dt, read_all_times_data(option = "foz_year")),
  tar_target(simulated_temps_dt, read_all_times_data(option = "simulated")),
  
  #Subset temperature data for fitting period
  tar_target(fitting_fortaleza_season_temps_dt, 
             read_fitting_period_data(location_all_data = fortaleza_season_temps_dt, 
                                        option = "fortaleza_season")),
  tar_target(fitting_fortaleza_2012_temps_dt, 
             read_fitting_period_data(location_all_data = fortaleza_2012_temps_dt, 
                                      option = "fortaleza_year")),
  tar_target(fitting_singapore_season_temps_dt, 
             read_fitting_period_data(location_all_data = singapore_season_temps_dt, 
                                      option = "singapore_season")),
  tar_target(fitting_foz_2012_temps_dt, 
             read_fitting_period_data(location_all_data = foz_2012_temps_dt, 
                                      option = "foz_year")),
  tar_target(fitting_simulated_temps_dt, 
             read_fitting_period_data(location_all_data = simulated_temps_dt, 
                                      option = "simulated")),
  
  #Estimate generation time distributions for each setting using 
    # i) no stage-specific approach
    # ii) stage-specific approach
  tar_target(fortaleza_season_ignoring_transm_cycle_taus, run_time_varying_monte_carlo(number_mc_samples = 2000,
                                                                                       fitting_data = fitting_fortaleza_season_temps_dt,
                                                                                       all_data = fortaleza_season_temps_dt,
                                                                                       transm_cycle = FALSE,
                                                                                       siraj_2017_params)),
  tar_target(fortaleza_season_transm_cycle_taus, run_time_varying_monte_carlo(number_mc_samples = 2000,
                                                                              all_data = fortaleza_season_temps_dt,
                                                                              fitting_data = fitting_fortaleza_season_temps_dt,
                                                                              transm_cycle = TRUE,
                                                                              siraj_2017_params)),
  tar_target(fortaleza_2012_ignoring_transm_cycle_taus, run_time_varying_monte_carlo(number_mc_samples = 2000, 
                                                                                       fitting_data = fitting_fortaleza_2012_temps_dt,
                                                                                       all_data = fortaleza_2012_temps_dt,
                                                                                       transm_cycle = FALSE,
                                                                                       siraj_2017_params)),
  tar_target(fortaleza_2012_transm_cycle_taus, run_time_varying_monte_carlo(number_mc_samples = 2000, 
                                                                            fitting_data = fitting_fortaleza_2012_temps_dt,
                                                                            all_data = fortaleza_2012_temps_dt,
                                                                            transm_cycle = TRUE,
                                                                            siraj_2017_params)),
  tar_target(singapore_season_ignoring_transm_cycle_taus, run_time_varying_monte_carlo(number_mc_samples = 2000, 
                                                                                       fitting_data = fitting_singapore_season_temps_dt,
                                                                                       all_data = singapore_season_temps_dt,
                                                                                       transm_cycle = FALSE,
                                                                                       siraj_2017_params)),
  tar_target(singapore_season_transm_cycle_taus, run_time_varying_monte_carlo(number_mc_samples = 2000, 
                                                                              fitting_data = fitting_singapore_season_temps_dt,
                                                                              all_data = singapore_season_temps_dt,
                                                                              transm_cycle = TRUE,
                                                                              siraj_2017_params)),
  tar_target(foz_2012_ignoring_transm_cycle_taus, run_time_varying_monte_carlo(number_mc_samples = 2000, 
                                                                                     fitting_data = fitting_foz_2012_temps_dt,
                                                                                     all_data = foz_2012_temps_dt,
                                                                                     transm_cycle = FALSE,
                                                                                     siraj_2017_params)),
  tar_target(foz_2012_transm_cycle_taus, run_time_varying_monte_carlo(number_mc_samples = 2000, 
                                                                            fitting_data = fitting_foz_2012_temps_dt,
                                                                            all_data = foz_2012_temps_dt,
                                                                            transm_cycle = TRUE,
                                                                            siraj_2017_params)),
  tar_target(simulated_ignoring_transm_cycle_taus, run_time_varying_monte_carlo(number_mc_samples = 2000, 
                                                                               fitting_data = fitting_simulated_temps_dt,
                                                                               all_data = simulated_temps_dt,
                                                                               transm_cycle = FALSE,
                                                                               siraj_2017_params)),
  tar_target(simulated_transm_cycle_taus, run_time_varying_monte_carlo(number_mc_samples = 2000, 
                                                                      fitting_data = fitting_simulated_temps_dt,
                                                                      all_data = simulated_temps_dt,
                                                                      transm_cycle = TRUE,
                                                                      siraj_2017_params))
)

