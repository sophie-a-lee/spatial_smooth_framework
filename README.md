# spatial_smooth_framework

Data and R code to support Lee S.A. et al., *(in press).* A Bayesian modelling framework to quantify multiple sources of spatial variation for disease mapping.

This work presents a novel Bayesian statistical modelling framework that allows multiple sources of spatial connectivity (for example, distance-based and human movement-based) into a spatial statistical model and quantifies the relative contribution of each term to the overall spatial structure of the model. This is achieved by applying penalised smoothing splines to coordinates that describe the relative connectivity between areas or observations to create 2-dimensional smooth surfaces describing the spatial structure of the data. Connectivity coordinates can be generated from any continuous measure of connectivity, including distance or the number of people travelling between cities. Smooth surfaces can be incorporated into Bayesian hierarchical models where their interpretation is similar to traditional random effects. 

### List of files in the repository:

* In data:

  * dengue_south_brazil.csv: Data table containing dengue case data for municipalities 2001 - 2020 in South Brazil (source: https://datasus.saude.gov.br/informacoes-de-saude-tabnet/). Other variables included: municipality names and codes, population, and the goegraphical coordinates of the centroid of the municipality.
  * regic_network_south.rds: Shapefile containing the Brazilian urban network, consisting of 6,637 connections between municipalities in South Brazil
  * south_brazil_shp: Shapefile for South Brazil

* In R:

  * 00_load_data_functions.R: load packages and data into the R session, and loads functions used in analysis
  * 00_human_movement_coords.R: generates connectivity coordinates based on the number of people expected to move between municipalities
  * Files beginning '01_' simulate fictional data 
  * Files beginning '02_' fits spatial statistical models
  * Files beginning '03_' generate visualisations and goodness-of-fit statistics from model output
  * 04_sensitivity_nimble_bym2.R: sensitivity analysis, fitting and comparing spatial statistical models
