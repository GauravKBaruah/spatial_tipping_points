

# Abrupt transitions and its indicators in mutualistic meta-networks: effects of network topology, size of metacommunities and species dispersal
Author: G Baruah, email: gbaruahecoevo@gmail.com

# R scripts

1. `01_Code_spatial_ews.R` this script deals with the overall simulation of all the treatments and organizes the results into a dataframe.
2. `02_code_species_spatial_ews.R` - this script deals with simulation of all the treatments but takes data at the species level of each meta-network and organizes into a data frame.

3. `buildlandscape.R` - R script for building metacommunity landscape.
4. `dispKernels.R` - R script for the kernel of dispersal used in the manuscript.

5. `Functions.R` -  R script of all the functions used in the model simulations.

# R data

1. example_network_collapse.RData  -  this is an example data for a mutualistic network forced to collapse that could be used to produce figure S6 in the manuscript.

2. scale_of_ews_webs_all_data.RData - this is the data for all 56 webs used in a meta-network framework forced to collapse originated from using the `01_Code_spatial_ews.R` script and used to produce figures 2 and 3 in the main-text.

3. Species_tipping_points_spatial_data_apr15.RData - this is the data for all 56 webs used in a meta-network framework but statistics are measured at the level of the species used to produce figure 4 in the main-text and originated from using the `02_code_species_spatial_ews.R` script. 

# Paper

1. `scale_transition.Rmd` is the rmarkdown file. In this R markdown script, the entire paper is written, annotated with codes, comments, and citation formats.
