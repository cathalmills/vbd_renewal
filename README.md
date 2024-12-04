# Renewal equations for vector-borne diseases

The following repository supports the manuscript "Renewal Equations for vector-borne diseases".

Fully reproducible code has been structured in the format of targets, a Make-like pipeline tool/package in R.

## Notes on application:
  - The application to simulated and real-world data demonstrates calculation of the time-varying generation time distribution under simplifying assumptions of temperature-dependent generation time distributions between stages of the transmission cycle (see functions.R).
  - The objective of the applications is not to estimate necessarily realistic generation time distributions, but rather to demonstrate that even if we assume temperature-dependent conditions, a stage-specific approach is required.
  - Please note that the framework allows for more general time-varying generation time distributions by changing the birth process functions that control the generation time distribution. See manuscript for further details.
  - Studied locations: i) Ang Mo Kio (a neighbourhood in Singapore), ii) São João do Tauape (a neighbourhood in Fortaleza, Brazil), and iii) Foz do Iguaçu (a city in Brazil).
  - Temperature data sources: The ERA5 global reanalysis dataset and Meterological Service Singpaore.

## References: 
-  W. M. Landau. The targets R package: a dynamic Make-like function-oriented pipeline
  toolkit for reproducibility and high-performance computing. Journal of Open Source Soft-
  ware, 6(57):2959, 2021. URL https://doi.org/10.21105/joss.02959.
- Hersbach, B. Bell, P. Berrisford, S. Hirahara, A. Horányi, J. Muñoz-Sabater, J. Nicolas,
  C. Peubey, R. Radu, D. Schepers, A. Simmons, C. Soci, S. Abdalla, X. Abellan, G. Bal-
  samo, P. Bechtold, G. Biavati, J. Bidlot, M. Bonavita, G. De Chiara, P. Dahlgren, D. Dee,
  M. Diamantakis, R. Dragani, J. Flemming, R. Forbes, M. Fuentes, A. Geer, L. Haim-
  berger, S. Healy, R. J. Hogan, E. Hólm, M. Janisková, S. Keeley, P. Laloyaux, P. Lopez,
  C. Lupu, G. Radnoti, P. De Rosnay, I. Rozum, F. Vamborg, S. Villaume, and J. Thépaut.
  The ERA5 global reanalysis. Quarterly Journal of the Royal Meteorological Society, 146
  (730):1999–2049, July 2020. ISSN 0035-9009, 1477-870X. doi: 10.1002/qj.3803. URL
  https://rmets.onlinelibrary.wiley.com/doi/10.1002/qj.3803.
- M. S. Singapore. Historical Daily Records, Aug. 2024. URL https://www.weather.
  gov.sg/climate-historical-daily/.
