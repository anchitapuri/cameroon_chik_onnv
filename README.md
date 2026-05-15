# Analysis for ONNV and CHIKV in Cameroon 

A statistical framework using multi-pathogen serology data to simultaneously reconstruct cross-reactivity patterns and infer individual- and population-level infection histories. In addition, this repository includes spatial modelling of CHIKV and ONNV transmission across Cameroon.


# Pre processing data 
District level geometery from: 
-  Caedistricts179_region.shp (179 districts, 183 geometeries (MANOKA == 5 geometries) 
-  cmr_admin3.shp (360 districts)
-  Districts not present in either were spatially assinged to their nearest shapefile polygon
-  Population weighted centroids calculated for each districts, and used for downstream analysis

# MultiSero Model 
- Model (MultiSero_Model.stan) used to jointly inferred pathogen-specific 
  prevalence estimates and between-pathogen cross-reactivity 
- Analysis in [`R/MultiSero_Fitting.R`](R/MultiSero_Fitting.R)


# Spatial Analysis
- To explore the variability of ONNV prevalence across Cameroon, applied a spatially explicit catalytic model with a Bayesian framework implemented using INLA
- Analysis for this part in [`R/INLA_SpatialPrediction.R`](R/INLA_SpatialPrediction.R)


