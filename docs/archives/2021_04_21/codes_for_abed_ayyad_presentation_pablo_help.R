#The following script explains how to repeat the steps which I presented during my talk to the Calgary R Users Group
# on Wednesday, 21 April 2021. The aim of the talk was to give an idea of how solar energy data can be used in R and 
# some easy ways to access it. I will do my best to explain the following code with comments, but it probably makes 
# sense to have a look at the slides (should be published by the Calgary R Users' Group soon enough) and especially the
# bibliography/references slide at the end to learn more. 

#Perhaps now is a good time to say: one major limitation of what I have placed down here is that I use monthly 
# values. I would argue that these are not great, but still valuable in some ways 

#Finally, let me thank Pablo Adames from the Calgary R Users' Group from testing the codes so thoroughly and 
# pointing out a few flaws which I acknowledge below with (PA)

#First, I am going to load all of the packages which we might need 
# since this script will end with a map, I will load quite a few graphics and GIS-related packages 
library(rgeos)
library(maptools)
library(ggplot2)
library(sp)
library(rgdal)
library(raster)
library(scales) # (all of the above) To make pretty maps 
library(plyr) # To bind spatial data to the maps 
library(nasapower) #To extract meteorological data etc 

# If you want to recreate the world cities/insolation map in the first part of the talk, please go 
# to the very end of this script 


#Here's some data that we will need to do the modeling with time 

silicon_cell_parameters = c("T_NOCT", "eta_STC", "beta_STC", "U_0" , "U_1")
mSi_values_Koehl = c(45, 18.4, -0.38, 30.02, 6.28)
pSi_values_Koehl = c(46, 14.1, -0.45, 30.02, 6.28)
aSi_values_Koehl = c(46, 6, -0.19, 25.73, 10.67)
uSi_values_Koehl = c(45, 9.5, -0.24, 30.02, 6.28)
CdTe_values_Koehl = c(45, 10.7, -0.25, 23.37, 5.44)

koehl_table = data.frame(cbind(mSi_values_Koehl, pSi_values_Koehl, aSi_values_Koehl, uSi_values_Koehl, CdTe_values_Koehl), stringsAsFactors = FALSE)
colnames(koehl_table) = c("mSi", "pSi", "aSi", "microSi", "CdTe") 
rownames(koehl_table) = silicon_cell_parameters



#I want to create the functions which take in the meteorological data and give us usable values
# I will make functions 
# As explained, PV cell efficiencies for a given technology are a function of cell temperature
# which itself is a function of ambient temperature (at 2 m height), irradiance and wind speed 


#The function below takes ambient temperatures, irradiances, windspeeds and the specification of the 
# type of technology 

koehl_model_monthly_cell_temperatures <- function(ambients_df, irradiances_df, windspeeds_df, techno_specification)
{
  selected_technology = dplyr::select(koehl_table, techno_specification)
  u_0 = selected_technology[4,]
  u_1 = selected_technology[5,]
  #This function will return the monthly cell temperatures 
  monthly_cell_temperatures = matrix(0, nrow = nrow(ambients_df), ncol = 12)
  
  for(i in 1:nrow(windspeeds_df))
  {
    for(j in 1:12)
    {
      monthly_cell_temperatures[i,j] = ambients_df[i,j] + (irradiances_df[i,j])/(windspeeds_df[i,j] + u_0 + u_1)
    }
  }
  
  return(monthly_cell_temperatures)
}

#The monthly cell temperatures feed into the efficiency of the cells 
# NB: this was fixed compared to the version presented at the talk 

create_monthly_etas <- function(cell_temperatures_df, irradiances_df)
{
  eta_STC = 0.154
  pv_gamma_ref = 0.12
  beta_ref = 0.0045
  eta_cellsDF = data.frame(matrix(0, nrow = nrow(irradiances_df), ncol = 12), stringsAsFactors = FALSE)
  
  for(i in 1:nrow(irradiances_df))
  {
    for(j in 1:12)
    {
      temp_prime = cell_temperatures_df[i,j] - 25
      eta_cellsDF[i,j] = eta_STC*(1 - beta_ref*(temp_prime)) + eta_STC*pv_gamma_ref*log(irradiances_df[i,j], base = 10)
    }
  }
  return(eta_cellsDF)  
}

#We need insolation, length of daylight hours, temperature and wind speed data to make this work 
# We may as well use a function for that, too, but note that we can only download 3 parameters at a time 
# from the NASA API 


use_nasapower_GHI_temps_dayhours <- function(locations_df)
{
  monthly_ghi = data.frame(matrix(0, nrow = nrow(locations_df), ncol = 13), stringsAsFactors = FALSE)
  monthly_t2m = data.frame(matrix(0, nrow = nrow(locations_df), ncol = 13), stringsAsFactors = FALSE)
  monthly_dayhours = data.frame(matrix(0, nrow = nrow(locations_df), ncol = 12), stringsAsFactors = FALSE)
  
  # "ALLSKY_SFC_SW_DWN" is the NASA POWER word for "global horizontal irradiance" 
  
  for(i in 1:nrow(locations_df))
  {
    full_request = nasapower::get_power(community = "SSE", pars = c("ALLSKY_SFC_SW_DWN", "T2M", "SG_DAY_HOUR_AVG"), 
                                        temporal_average = "CLIMATOLOGY", lonlat = c(locations_df$LON[i], locations_df$LAT[i]))
    
    ghi_request = full_request[which(full_request$PARAMETER == "ALLSKY_SFC_SW_DWN"),]
    t2m_request = full_request[which(full_request$PARAMETER == "T2M"),]
    dayhours_request = full_request[which(full_request$PARAMETER == "SG_DAY_HOUR_AVG"),]
    monthly_ghi[i,] = data.frame(ghi_request[,4:16])
    monthly_t2m[i,] = data.frame(t2m_request[,4:16])
    monthly_dayhours[i,] = data.frame(dayhours_request[,4:15])
  }
  
  list_of_dfs <- list(monthly_ghi, monthly_t2m, monthly_dayhours)
  print("This is a list which includes GHI, temperatures at 2 M and daylight hours, in that order")
  #Give ourselves a reminder when we use this code, so we extract things in the right order as the need arises 
  return(list_of_dfs)
}

use_nasapower_windspeeds <- function(locations_df)
{
  windspeed_df = data.frame(matrix(0, nrow = nrow(locations_df), ncol = 13))
  
  for(i in 1:nrow(windspeed_df))
  {
    wind_request = nasapower::get_power(community = "SSE", pars = "WS2M", temporal_average = "CLIMATOLOGY", lonlat = c(locations_df$LON[i], locations_df$LAT[i]))
    windspeed_df[i,] = data.frame(wind_request[,4:16])
  }
  return(windspeed_df)
}

#Please bear in mind: the units in which these values are given are important. I hope that the way
# this emerges in the script will be "natural" enough for most of you, but for now we plough on


#Okay, in order to use these functions, we will need to have data frames which hold the longitude and latitude data 
# How do we do that? I happen to know that there are 19 different districts in Alberta ...but first! 

#Canada has ISO country code "CA" or "CAN" sometimes; let's download a level-2 map from the Global Administrative Database
# Level 2 means 1 level below a province (level 0 is just the country, etc)
canada_map_level_2 = getData(name = "GADM", country = "CA", level = 2)
#You should see that there are 13 separate Level 1 provinces/territories in Canada, as found by: 
unique(canada_map_level_2$NAME_1)


#Note above that the "Spatial Polygons Data Frame" class of object, of which 
# the Canada level 2 and Alberta Level 2 (which we will create now) map are 
# examples, inherits the subsetting methods from data frame:

alberta_map_level_2 = canada_map_level_2[which(canada_map_level_2$NAME_1 == "Alberta"),]

#In Alberta, the names are sadly not so illustrative (run the next line and see): 
unique(alberta_map_level_2$NAME_2)
#How sad. Instead, we can create a vector of names which we can use in later display tables to show us what 
# the districts are, based on what the largest city in each district is (per Wikipedia)
alberta_district_names = c("Medicine Hat", "Lethbridge", "Claresholm", "Hanna", "Strathmore", "Calgary", "Wainwright", "Red Deer", "Rocky Mountain House", "Lloydminster", "Edmonton", "Cold Lake", "White Court", "Hinton", "Canmore", "Fort McMurray", "Slave Lake", "Grande Cache", "Grande Prairie" )

#We could just plot alberta_map_level_2, and in the end we will do that, colouring in the districts 
# with values from the amount of PV yield or LCOE (more later), but for now we just want to extract the 
# longitude and latitude coordinates of each district within Alberta. To do that, notice that 
# alberta_map_level_2 has multiple "polygons," each representing the shape file of a district. 
# to get the center of each such polygon, we do like so: 

## data.frame(gCentroid(alberta_map_level_2[i,]))
# Where "i" is our index running from 1 to 19, for the number of districts in Alberta
# gCentroid is from the "rgeos" package and allows us to find the center of a polygon
alberta_coordinates = data.frame(matrix(0, nrow = length(alberta_district_names), ncol = 2), 
                                  stringsAsFactors = FALSE)

#If you're wondering, having strings as factors set to false allows us to work more smoothly
# with numeric data in a data frame
for(i in 1:nrow(alberta_coordinates))
{
  alberta_coordinates[i,] = data.frame(gCentroid(alberta_map_level_2[i,]))
}

#We want to give the Alberta coordinates their correct names
rownames(alberta_coordinates) = alberta_district_names
colnames(alberta_coordinates) = c("LON", "LAT")
#Please watch out for the column names: they need to correspond to what we feed into the function for nasapower 
# data 

#Let's now download the raw data for these Alberta districts 
alberta_first_request = use_nasapower_GHI_temps_dayhours(alberta_coordinates)
alberta_windspeeds = use_nasapower_windspeeds(alberta_coordinates)

alberta_ghi_raw = alberta_first_request[[1]]
alberta_temps_2M_unadjusted = alberta_first_request[[2]]
alberta_daylength = alberta_first_request[[3]]

#It's worth doing a quick spot check on these: for all rows, does the length of day in each month (column) get 
# longer from January to July? Do the temperatures increase, too? Look at column 13 of the raw GHI table. It is
# given in units of kWh/(m^2 per day). The formulas for cell temperature and efficiency, however, use irradiance
# in units of Watts. Let's just convert that real quick, and assume that the insolation falls in neat monthly 
# "averages" for now. 

calgary.rug.alberta.irradiances_monthly = 1000*alberta_ghi_raw[,1:12]/alberta_daylength
calgary.rug.alberta.ambient.temperatures_adjusted_monthly = alberta_temps_2M_unadjusted[,1:12] + 0.4*abs(alberta_temps_2M_unadjusted[,1:12])

#It's also worth adjusting the GHI insolations into monthly values, which we will use when we calculate the yields 
days_in_months = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
#Make sure you got 12 there! 
length(days_in_months) == 12

#We will use this when we calculate yields in a minute 
calgary.rug.alberta.insolation_monthly = data.frame(matrix(0, nrow = nrow(calgary.rug.alberta.irradiances_monthly), ncol = 12), stringsAsFactors = FALSE)
for(i in 1:nrow(calgary.rug.alberta.insolation_monthly))
{
  calgary.rug.alberta.insolation_monthly[i,] = alberta_ghi_raw[i,1:12]*days_in_months
}

#Let me just apologise for my very inconsistent naming conventions, it's just that I wrote these codes in bursts in between doing other things
calgary.rug.alberta.cell_temps_monthly = koehl_model_monthly_cell_temperatures(ambients_df = calgary.rug.alberta.ambient.temperatures_adjusted_monthly, 
                                                                               irradiances_df = calgary.rug.alberta.irradiances_monthly, alberta_windspeeds, "mSi")

calgary.rug.alberta.etas_monthly = create_monthly_etas(cell_temperatures_df = calgary.rug.alberta.cell_temps_monthly, 
                                                       irradiances_df = calgary.rug.alberta.irradiances_monthly)

#You should see that, in Alberta, the efficiency of a solar cell should be swimming around 20%, and will increase in the colder months, as expected
# Now, these are for "cells" and there is a lot that we are leaving out in terms of yields, but let's plough on and we can deal with some of these
# simplifications later. 

calgary.rug.yields.monthly_1m_squared = data.frame(matrix(0, nrow = nrow(calgary.rug.alberta.etas_monthly), ncol = 12), stringsAsFactors = FALSE)
for(i in 1:nrow(calgary.rug.yields.monthly_1m_squared))
{
  for(j in 1:12)
  {
    calgary.rug.yields.monthly_1m_squared[i,j] = calgary.rug.alberta.etas_monthly[i,j]*calgary.rug.alberta.insolation_monthly[i,j]
  }
}

#I should point out, with data frames it's possible to simply multiply two DFs with just A * B, but I like doing this explicitly because it keeps 
# alert about multiplying matrices 

#So what we have now is what 1 meter squared of PV cells would yield over the course of a "typical" year. Do a quick check that the annual output of that last 
# DF is roughly 20% of the original display that we had at the begining of the slides--you should values ranging from 270 (270 kWh per month) to 224, use 
# range(rowSums(DF_name_here))

#Now, this is also kind of silly, because we don't really ever have "1 m squared of cells," instead we usually have panels of about 1.7 m^2 (or so). So, the first 
# thing we do is to build a data.frame which has the expected from one panel; but a panel's surface is not all just "cells," there is spacing between them. Let's 
# also knock off 5% for these losses. I will also make the area of a variable, so that it's easier for you to change this code when you feel like it. 

area_of_a_panel = 1.7 
calgary.rug.yields.monthly_1_panel = 0.95*area_of_a_panel*calgary.rug.yields.monthly_1m_squared

#Before anybody says anything: yes, if you're interested in helping build towards a CRAN package, do please email me. 

#Okay, so now that's good for panels, but what if you wanted to have a 1 MW (Megawatt) solar farm in Alberta? How much could we expect in PV yields, then?
# Those 1.7 m^2 panels have a power rating of 330 W (roughly). 

single_panel_power = 330 
panels_for_1MW_farm = (1*10^6)/single_panel_power

calgary.rug.yields.monthly_1MW_farm = panels_for_1MW_farm*calgary.rug.yields.monthly_1_panel

#These are expressed in terms of kWh however and since we're really fancy and own a 1 MW solar PV farm, we ought to switch to using MWh. We should also 
# then check to see what the capacity factor is now--if you remember from my talk, 1 MW solar farm has a maximum output of 8.76 GWh (8760 MWh), and you
# should see that Alberta performs better than the world average in 2018 (13.25%) for solar PV, but that, of course, there are limits for 
# solar 

calgary.rug.yields.monthly_1MW_farm_MWh = calgary.rug.yields.monthly_1MW_farm/1000
calgary.rug.capacity_factors = rowSums(calgary.rug.yields.monthly_1MW_farm_MWh)/8760

#...and I should probably not ignore this, it will come in handy later 
calgary.rug.annual.yields.1_MW_farm_MWh = rowSums(calgary.rug.yields.monthly_1MW_farm_MWh)


#Short of increasing the number of daylight hours per day, one way we have available to us to increase the capacity factor for solar PV in Alberta is 
# to use tilted panels and capture plane of array irradiance instead of GHI--see the slides for examples--but I'm leaving that as an exercise for the student. 
# Now we move on the fun part. 

##Calculating LCOE for a 1 MW solar farm in different districts of Alberta. As explained, the LCOE is the cost/(energy_yields), using a Net Present Value for 
# both (NPV). The higher the yields (capacity factor) at a given location, the lower the cost. Likewise, for solar farms, land area footprint is a huge deal 
# in relation to "thermal" generators. Stat Can has placed the value of an acre of land in Alberta at CA$ 2,842; that's actually quite favourable compared 
# to other provinces and territories. Let's start there, but I reckon we add a premium of 40% at least.

cost_of_acre_alberta = 1.4*2842 

#...but of course, there should be a difference across all the districts within Alberta. 


calgary.rug.alberta.land.costs = rnorm(n = length(calgary.rug.capacity_factors), cost_of_acre_alberta, sd = 0.25*cost_of_acre_alberta)
#again, I'm doing this because it should be easy to see where to make things modular

#...but we don't buy land in acres, we buy it in square meters needed for our farm. Let's say we need to expand the area covered by the panels by 15% to 
# cover inter-row spacing; note that we would likely need even more land if the panels were tilted, to prevent shading between rows 
calgary.rug.alberta.area.needed_1MW_farm_msquared = panels_for_1MW_farm*1.72
acre_of_land = 4046.86 #Google is your friend, remember 

calgary.rug.alberta.acres.needed_1MW_farm = calgary.rug.alberta.area.needed_1MW_farm_msquared/acre_of_land

#You should get something close to 1.3 acres; let's just round it up to 1.5 in our later calculations 
# Note that this kind of only works because we're using flat panels collecting GHI; we would 
# need much more land if we were tilting our panels. This is an important point!

calgary.rug.alberta.land.costs_1MW_farm = 1.5*calgary.rug.alberta.land.costs

#So now, we have to calculate two other costs for solar PV farms: the installation costs, which are usually given in units of $ per kW "peak" ($/kW_p); we also
# need to calculate the costs of operations and maintenance, usually given as a percentage of the capital costs. The Canada Energy Regulator lists different 
# capital costs for utility scale solar PV; I'm going to go with CA$ 1,400 per kW_p, or $1.4 million for our 1 MW of panels. I will then also fit in 7% O & M 
# costs per year. 

calgary.rug.capital.costs = 1.4*10^6

#If you've been paying attention this far, make sure to load package "FinCal" for financial calculations
# including our NPV

library(FinCal)

# We could spend a lot of time fidgeting around with the NPV calculations, but I prefer to just use the built-in functions
# in R. You could of course look up the more extensive definition. For now, it's important to know that the function "npv" 
# in FinCal takes three arguments: the "Year 0" cashflows, which are our investments/capital costs; a vector which contains 
# the expenditures going forward; and an interest rate or discount rate. We will also just treat units of electricity 
# generated--MWh--like they're # Canadian dollars and we move on from there. Finally, you would need to know how long
# your project is expected to last. I'd say 20 years is reasonable for land-based solar PV systems (and honestly, they're
# now 25 years)

project_lifetime = 25

#Let's say we get a discount rate of 2%
discount_rate = 0.018
#It's worth pointing out that there is a function "discount.rate" in-built in FinCal

#We want 7% of the capital outlay to be 
calgary.rug.operations_and_maintenance_costs = rep(0.05*calgary.rug.capital.costs, project_lifetime)
#Notice that in our case here, we're just going to assume that it's the same O & M cost in each district within Alberta
# although I have a suspicion this wouldn't be true in "the real world" 


#I should have also mentioned that the NPV takes a "Year 0" in which the investments/capital costs are made but there is
# no electricity generation 

generate_LCOE <- function(rate_of_money, capital_outlays, operating_costs, yields)
{
  npv_money = FinCal::npv(r = rate_of_money, cf = c(capital_outlays, operating_costs))
  npv_electricity = FinCal::npv(r = rate_of_money, cf = c(0, rep(yields, project_lifetime)))
  lcoe_to_return = npv_money/npv_electricity
  return(lcoe_to_return)
}

#Now we run this function for each district of Alberta 
#First, I make a vector to hold the LCOE value for each district in Alberta
calgary.rug.lcoe.per_district = numeric(length = nrow(alberta_coordinates))

for(i in 1:length(calgary.rug.lcoe.per_district))
{
  capital_costs_here = calgary.rug.capital.costs + calgary.rug.alberta.land.costs[i]
  #Obviously, the cost of land has become negligible in my example, going on Stat Can's very low cost of land
  # ...but you now have the tools to do this yourself at home 
  calgary.rug.lcoe.per_district[i] = generate_LCOE(discount_rate, capital_costs_here, calgary.rug.operations_and_maintenance_costs, calgary.rug.annual.yields.1_MW_farm_MWh[i])
  
}

#Remember, to get the more human friendly "cents per kWh," simply divide by 1000

#If you punch everything in, you actually get a number which might at first seem high. Remember a few things though. 
# The numbers are actually not that bad compared to what is mentioned by the US' DOE: https://www.energy.gov/eere/solar/sunshot-2030
# Factor in, of course, the 0.8 conversion from CAD to USD and you start getting values that are not entirely out of this world

calgary.rug.lcoe.per_district/1000

#This is all kind of contrived data; I think I was being super pessimistic in some ways and that we could 
# get a much better LCOE if we tried the right combinations of things. Let's suppose though that we wanted to run with
# this. Could we make a map of Alberta which shows us which areas are more likely to be better sites for solar PV?

#First, let's get everything into one data frame. I also need to make sure that one key in this data frame matches
# what we have 
calgary.rug.solar.lcoe.dataframe = data.frame(unique(alberta_map_level_2@data$NAME_2), calgary.rug.lcoe.per_district, alberta_coordinates$LON, alberta_coordinates$LAT  , stringsAsFactors = FALSE)
colnames(calgary.rug.solar.lcoe.dataframe) = c("NAME_2", "LCOE", "lon", "lat")
#rownames(calgary.rug.solar.lcoe.dataframe) = row.names(alberta_coordinates)
#You might think that this was a bit redunant ... maybe, but it's useful for later

#Now, let's take that important column of data and make sure that it sits nestled in the right place for the geographical
# data

#We have to go through a few weird steps at this stage; the aim here is to create a polygon which can 
# spread out the values of our contrived "LCOE" 

#Let's "stretch" Alberta out properly
WGS84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
#Thank you (PA) for pointing out that this line was missing in the original code
# PA also used the following, which might work for you: 
# proj4string(alberta_map_level_2)
# alberta_map_transformed = spTransform(alberta_map_level_2, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

#This is what I now did: 
alberta_map_transformed = spTransform(alberta_map_level_2, WGS84)
#Let's give ourselves an id column
alberta_map_transformed@data$id = rownames(alberta_map_transformed@data)
#Now we do the important step of attaching the LCOE values to each district 
alberta_map_transformed@data = join(alberta_map_transformed@data, calgary.rug.solar.lcoe.dataframe, by = "NAME_2")
#Now, we turn the map into a dataframe
alberta_df_lcoe = fortify(alberta_map_transformed)
#Finally, we make sure that the data frame covers all points on the map 
alberta_df_lcoe = join(alberta_df_lcoe, alberta_map_transformed@data, by = "id")
#Incidentally, you do NOT want to run the above commands multiple times, please follow the steps exactly and carefully

#Now make sure to run this line; it actually depends on your version of ggplot, but it's probably good practice anyway 
colnames(alberta_df_lcoe) <- make.unique(names(alberta_df_lcoe))



#Thanks PA

#Now we can plot a pretty map showing relative LCOE

ggplot() + geom_polygon(data = alberta_df_lcoe, aes(x = long, y = lat, group = group, fill = LCOE), size = 0.12) + scale_alpha_continuous(name="Contrived LCOE data", palette = "Spectral", breaks = 7)



##Clearly, there is so much that is "contrived" about the data I have presented above--and I stand by the claim that 
# I am seriously overestimating LCOE from solar PV. What's important to understand however is that meteorological factors
# and economic as well as technological factors come together to create a landscape where deploying solar PV becomes a 
# less or more favourable proposition based on where you are in the world. Remember: if you own your own home/roof, 
# go ahead and install those panels, it would be hard for you to end up ignoring them, but try to have the panels facing 
# south. 

##As a kind of bonus, the following code will produce the map of world cities which I showed early on in the
# slides 

#To get the longitude, latitude for the cities in question, you will need to install nominatim for Open Street Maps
#devtools::install_github("hrbrmstr/nominatim")
# You will ned an OSM API key which you can get for free from Map Quest 
# Then you just populate the data frame like so 

calgary.rug.cities = c("Calgary", "Vancouver", "Edmonton", "Yellowknife", "Winnipeg", "Utrecht", "Buenos Aires", "Ramallah", "Doha", "Berlin")
calgary.rug.countries = c("CA", "CA", "CA", "CA", "CA", "NL", "AR", "PS", "QA", "DE")
#You need to know the ISO country codes for the countries in question, I suppose 


calgary.rug.solar.df = data.frame(matrix(nrow = length(calgary.rug.cities), ncol = 2))
for(i in 1:length(calgary.rug.cities))
{
  place_holder = nominatim::osm_geocode(query = calgary.rug.cities[i], country_codes = calgary.rug.countries[i], key = myosmapikey)
  calgary.rug.solar.df[i,1] = place_holder$lon
  calgary.rug.solar.df[i,2] = place_holder$lat
  print("One more city")
}
#We will need 13 columns of insolation data, one for each month and one for the annual values
calgary.rug.solar.monthly = data.frame(matrix(nrow = nrow(calgary.rug.solar.df), ncol = 13))
for(i in 1:nrow(calgary.rug.solar.monthly))
{
  request_solar_ghi = nasapower::get_power(community = "SSE", pars = "ALLSKY_SFC_SW_DWN", temporal_average = "CLIMATOLOGY", lonlat = c(calgary.rug.solar.df[i,1], calgary.rug.solar.df[i,2]))
  #To see how this next line works, you really have to know what the structure of the returned response from NASAPOWER
  calgary.rug.solar.monthly[i,] = data.frame(request_solar_ghi[which(request_solar_ghi$PARAMETER == "ALLSKY_SFC_SW_DWN"),][,4:16])
}
colnames(calgary.rug.solar.monthly) = colnames(data.frame(request_solar_ghi[,4:16]))
#Check this to see what it is
calgary.rug.solar.monthly
rownames(calgary.rug.solar.monthly) = calgary.rug.cities
calgary.rug.solar.df = cbind(calgary.rug.solar.df, calgary.rug.solar.monthly$ANN)
calgary.rug.solar.df
row.names(calgary.rug.solar.df) = calgary.rug.cities
colnames(calgary.rug.solar.df) = c("LON", "LAT", "Annual GHI")
calgary.rug.solar.df

#The annual value is given as a "average day," so we need to fix it 
calgary.rug.solar.df$`Annual GHI` = calgary.rug.solar.df$`Annual GHI`*365
calgary.rug.solar.df

calgary.rug.cities.diffuse = vector(length = length(calgary.rug.cities))
for(i in 1:length(calgary.rug.cities.diffuse))
{
  diffuse_request = nasapower::get_power(community = "SSE", pars = "DIFF", temporal_average = "CLIMATOLOGY", lonlat = c(calgary.rug.cities.for.presenting$LON[i], calgary.rug.cities.for.presenting$LAT[i]))
  annual_diffuse = data.frame(diffuse_request)$ANN
  calgary.rug.cities.diffuse[i] = annual_diffuse
}

calgary.rug.cities.dnr = vector(length = length(calgary.rug.cities))
for(i in 1:length(calgary.rug.cities))
{
  dnr_request = nasapower::get_power(community = "SSE", pars = "DNR", temporal_average = "CLIMATOLOGY", lonlat = c(calgary.rug.cities.for.presenting$LON[i], calgary.rug.cities.for.presenting$LAT[i]))
  annual_dnr = data.frame(dnr_request)$ANN
  calgary.rug.cities.dnr[i] = annual_dnr
}

calgary.rug.cities.diffuse[is.na(calgary.rug.cities.diffuse)] <- 0
calgary.cities.diff_to_dnr = calgary.rug.cities.diffuse/calgary.rug.cities.dnr

#Yelloknife is a special case it seems; just not a lot of sunlight in some parts of the year 
calgary.rug.yellowknife_diffuse1 = data.frame(nasapower::get_power(community = "SSE", pars = "DIFF", temporal_average = "CLIMATOLOGY", lonlat = c(calgary.rug.cities.for.presenting$LON[4], calgary.rug.cities.for.presenting$LAT[4])))
calgary.rug.yellowknife_diffuse = sum((calgary.rug.yellowknife_diffuse1)[,4:16], na.rm = TRUE)
calgary.rug.yellowknife_dnr1 = data.frame(nasapower::get_power(community = "SSE", pars = "DNR", temporal_average = "CLIMATOLOGY", lonlat = c(calgary.rug.cities.for.presenting$LON[4], calgary.rug.cities.for.presenting$LAT[4])))
calgary.rug.yellowknife_dnr = sum(calgary.rug.yellowknife_dnr1[,4:16], na.rm = TRUE)
calgary.rug.yellowknife_diff_to_dnr = calgary.rug.yellowknife_diffuse/calgary.rug.yellowknife_dnr
calgary.cities.diff_to_dnr[4] = calgary.rug.yellowknife_diff_to_dnr


#Thanks once more to PA for mentioning that I left this out
calgary.rug.cities.for.presenting = cbind(calgary.rug.cities, calgary.rug.solar.df)

calgary.rug.cities.for.presenting$DiffDNR = calgary.cities.diff_to_dnr


#Before going any further, please have a quick check to see that the different data frames are reasonable before trying to plot them

myPaletter = colorRampPalette(rev(brewer.pal(11, "RdBu")))
scaled_map = scale_color_gradientn(colors = myPaletter(25), limits = c(min(calgary.rug.solar.df$DiffDNR2), max(calgary.rug.solar.df$DiffDNR2)), name = "Diffuse-to-DNR \n %-age")


#The map for diffuse to DNI 
ggplot2::ggplot(data = world_map) +
  geom_sf() +
  geom_point(data = calgary.rug.solar.df[1:10,],
             mapping = aes(x = LON, y = LAT,
                           color = (DiffDNR2)),
             size = 5) +   coord_sf(xlim = c(min(calgary.rug.cities.for.presenting$LON)-10,
                                             max(calgary.rug.cities.for.presenting$LON) + 10),
                                    ylim = c(min(calgary.rug.cities.for.presenting$LAT)-10,
                                             max(calgary.rug.cities.for.presenting$LAT)+10), expand = F) +
  scaled_map +
  ggtitle("Diffuse to DNR annually \n selected cities") +theme(panel.background = element_rect(fill = "lightblue",  color = "blue", size = 0.5))

