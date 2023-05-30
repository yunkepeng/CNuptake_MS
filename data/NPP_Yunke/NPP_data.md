
## A complete dataset of BP, ANPP, BNPP, C/N in tissues, Net minerlization, and site-level predictors, for the use of Peng et al. (first chapter).
### Author: Yunke Peng
### Data-use policy: the data has received all contributed authors (co-authors of MS) permission. It will be publicly available once this MS is published. But now it is still for restricted use.

### Variables
* lon (degree): longtitude 
* lat (degree): latitude 
* site: site name given in original data
* references: reference information of data
* z (m): elevation
* Begin_year: start year of measurement (from original data; if missing, using long-term period, e.g. 1991-2010) 
* End_year: end year of measurement (from original data; if missing, using long-term period, e.g. 1991-2010) 
* year_start: start year of measurement (from original data, but converting to 1980 if the measurement year is before 1980, for the use when interpolating site vcmax, fAPAR and climates) 
* year_end: end year of measurement (from original data, but converting to 1980 if the measurement year is before 1980, for the use when interpolating site vcmax, fAPAR and climates) 
* NPP.foliage (gC/m2/yr) :  Net-primary-production within foliage
* NPP.stem (gC/m2/yr) :  Net-primary-production within stem
* NPP.wood (gC/m2/yr) :  Net-primary-production within wood
* NPP.fine (gC/m2/yr) :  Net-primary-production within fine-root
* NPP.coarse (gC/m2/yr) :  Net-primary-production within coarse-root
* ANPP_2 (gC/m2/yr) :  Net-primary-production in aboveground (ANPP2 = NPP.foliage + NPP.wood)
* BNPP_1 (gC/m2/yr) :  Net-primary-production in belowground 
* TNPP_1 (gC/m2/yr) :  Net-primary-production in total (TNPP_1 = ANPP_2 + BNPP_1) 
* age (years) : measured age in years
* LAI (m2/m2)  : leaf area index
* observedfAPAR  : "measured" fAPAR as derived from measured LAI 
* GPP (gC/m2/yr) : Gross-primary-production in measurement
* soilCN : measured soil C/N
* pft : biomes, either in forest or grassland
* file : file name as derived from whose collection
* CN_stem_final : Stem C/N
* CN_wood_final : Wood C/N
* CN_leaf_final : Leaf C/N
* CN_root_final : Root C/N
* lnf_obs_final (gN/m2/yr): Leaf N flux (= NPP.foliage / leaf C/N )
* bnf_obs_final (gN/m2/yr): root N flux (= BNPP_2 / root C/N )
* wnf_obs_final (gN/m2/yr): wood N flux (= NPP.wood / wood C/N )
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
