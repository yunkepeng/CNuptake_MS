
## A complete dataset of net minerlization from Gill and Finzi et al. (publicly available):

Gill AL, Finzi AC. 2016. Belowground carbon flux links biogeochemical cycles and resource-use efficiency at the global scale. Ecology Letters 19: 1419–1428.


### Variables
* lon (degree): longtitude 
* lat (degree): latitude 
* site: site name given in original data
* references: reference information of data
* z (m): elevation
* Begin_year: start year of measurement (from original data; if missing, using long-term period, e.g. 1984-2013) 
* End_year: end year of measurement (from original data; if missing, using long-term period, e.g. 1984-2013) 
* year_start: start year of measurement (from original data, but converting to 1980 if the measurement year is before 1980, for the use when interpolating site vcmax, fAPAR and climates) 
* year_end: end year of measurement (from original data, but converting to 1980 if the measurement year is before 1980, for the use when interpolating site vcmax, fAPAR and climates) 
* pft : biomes, either in forest or grassland
* file : file name as derived from whose collection
* Nmin (gN/m2/yr): net minerlization (all derived from Gill & Finzi 2016; Ecology Letters)
* Tg (degree celcius): growth temperature
* PPFD (umol/m2/s): average growing-season PPFD
* Total PPFD (mol/m2): Total growing-season PPFD of the year
* vpd (kPa): vapor-pressure-deficient averaged from growing season
* fAPAR: fAPAR interpolated from fAPAR3g map product
* mapped_age (years): age interpolated from map product
* CNrt: soil C/N interpolated from map product
* LMA (g/m2): LMA interpolated from map product
* vcmax25 (umol/m2/s): vcmax25 interpolated from map product
* nhx (nh3 concentration)
* noy (noy concentration)
* ndep: sum of nh3 + noy = nitrogen deposition from ingestr pacakege (gN/m2/yr) (originally Lamarque et al. 2011)
