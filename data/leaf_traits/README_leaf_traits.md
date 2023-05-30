## The dataset has collected all publicly available vcmax data from multople references. Most of them was already published together with Peng et al. 2021 Communications Biology. But this dataset is the most updated version. 

## This dataset is free for further researchers to use them, since it is all already publicly avaialable.

### Variables
* lon (degree): longtitude 
* lat (degree): latitude 
* z (m): elevation
* start_yr: start year of measurement (from original data; if missing, using long-term period, e.g. 1991-2010) 
* end_yr: end year of measurement (from original data; if missing, using long-term period, e.g. 1991-2010) 
* site: site name given in original data
* sources: sources information of data
* species, family, genus: species, family and genus name
* Tleaf: leaf temperature
* narea: leaf nitrogen per area (g/m2)
* parea: leaf phosphorus per area (g/m2)
* lma: leaf mass per area (g/m2)
* vcmax: maximum rate of carboxylation (umol/m2/s)
* vcmax25: maximum rate of carboxylation at 25 degC (umol/m2/s)
* jmax: maximum rate of electron transport (umol/m2/s)
* jmax25: maximum rate of electron transport at 25 degC (umol/m2/s)
* Date: exact time for measurement collection
* Year: year of measurement
* carea: leaf carbon per area (g/m2)
* C_percent: leaf C percentage (%)
* Aarea: Light-saturated photosynthetic carbon assimilation per unit leaf area (umol/m2/s)

#below is reference information

   source                                                  number 
   <chr>                                                    <int>
 1 Atkin et al. 2015 New Phytologist                         1137
 2 Bahar et al 2017 New Phytologist                           301 (most of them was replicated to Atkin et al., therefore suggested to remove them when doing analysis)
 3 Bloomfield et al 2018 New Phytologist                      940
 4 Bloomfield; PhD unpublished                                244
 5 Cernusak et al 2011 Agricultural and Forest Meteorology     81
 6 Dong Ning collection                                       296
 7 Maire et al. 2015 GEB                                     2070 (vcmax and jmax were calculated basing on one-point method!)
 8 TROBIT; Thomas Domingues & Jon Lloyd                      2199
 9 Walker et al 2014  Ecology and Evolution                   194
10 Wang et al. 2017 Ecology                                   294
11 Xu et al. 2021 Tree Physiology                             427


#reference infrormation: 
Atkin O K, Bloomfield K J, Reich P B, Tjoelker M G, Asner G P, Bonal D, Bönisch G, Bradford M G, Cernusak L A, Cosio E G, Creek D, Crous K Y, Domingues T F, Dukes J S, Egerton J J G, Evans J R, Farquhar G D, Fyllas N M, Gauthier P P G, Gloor E, Gimeno T E, Griffin K L, Guerrieri R, Heskel M A, Huntingford C, Ishida F Y, Kattge J, Lambers H, Liddell M J, Lloyd J, Lusk C H, Martin R E, Maksimov A P, Maximov T C, Malhi Y, Medlyn B E, Meir P, Mercado L M, Mirotchnick N, Ng D, Niinemets Ü, O’Sullivan O S, Phillips O L, Poorter L, Poot P, Prentice I C, Salinas N, Rowland L M, Ryan M G, Sitch S, Slot M, Smith N G, Turnbull M H, Vanderwel M C, Valladares F, Veneklaas E J, Weerasinghe L K, Wirth C, Wright I J, Wythers K R, Xiang J, Xiang S and Zaragoza-Castells J 2015 Global variability in leaf respiration in relation to climate, plant functional types and leaf traits New Phytol. 206 614–36

Bahar N H A, Ishida F Y, Weerasinghe L K, Guerrieri R, O’Sullivan O S, Bloomfield K J, Asner G P, Martin R E, Lloyd J, Malhi Y, Phillips O L, Meir P, Salinas N, Cosio E G, Domingues T F, Quesada C A, Sinca F, Escudero Vega A, Zuloaga Ccorimanya P P, del Aguila-Pasquel J, Quispe Huaypar K, Cuba Torres I, Butrón Loayza R, Pelaez Tapia Y, Huaman Ovalle J, Long B M, Evans J R and Atkin O K 2017 Leaf-level photosynthetic capacity in lowland Amazonian and high-elevation Andean tropical moist forests of Peru New Phytol. 214 1002–18

Bloomfield K J, Prentice I C, Cernusak L A, Eamus D, Medlyn B E, Rumman R, Wright I J, Boer M M, Cale P, Cleverly J, Egerton J J G, Ellsworth D S, Evans B J, Hayes L S, Hutchinson M F, Liddell M J, Macfarlane C, Meyer W S, Togashi H F, Wardlaw T, Zhu L and Atkin O K 2019 The validity of optimal leaf traits modelled on environmental conditions New Phytol. 221 1409–23

Cernusak L A, Hutley L B, Beringer J, Holtum J A M and Turner B L 2011 Photosynthetic physiology of eucalypts along a sub-continental rainfall gradient in northern Australia Agric. For. Meteorol. 151 1462–70
Domingues T F, Meir P, Feldpausch T R, Saiz G, Veenendaal E M, Schrodt F, Bird M, Djagbletey G, Hien F, Compaore H, Diallo A, Grace J and Lloyd J 2010 Co-limitation of photosynthetic capacity by nitrogen and phosphorus in West Africa woodlands Plant, Cell Environ. 33 959–80

Dong N, Colin Prentice I, Evans B J, Caddy-Retalic S, Lowe A J and Wright I J 2017 Leaf nitrogen from first principles: Field evidence for adaptive variation with climate Biogeosciences 14 481–95
Ferreira Domingues T, Ishida F Y, Feldpausch T R, Grace J, Meir P, Saiz G, Sene O, Schrodt F, Sonké B, Taedoumg H, Veenendaal E M, Lewis S and Lloyd J 2015 Biome-specific effects of nitrogen and phosphorus on the photosynthetic characteristics of trees at a forest-savanna boundary in Cameroon Oecologia 178 659–72 Online: http://dx.doi.org/10.1007/s00442-015-3250-5

Maire V, Wright I J, Prentice I C, Batjes N H, Bhaskar R, van Bodegom P M, Cornwell W K, Ellsworth D, Niinemets Ü, Ordonez A, Reich P B and Santiago L S 2015 Global effects of soil and climate on leaf photosynthetic traits and rates Glob. Ecol. Biogeogr. 24 706–17

Meir P, Levy P E, Grace J, Jarvis P G, Ecology S P, Meir P, Levy P E, Grace J and Jarvis P G 2017 Photosynthetic parameters from two contrasting woody vegetation types in West Africa 192 277–87

Walker A P, Aranda I, Beckerman A P, Bown H, Cernusak L A, Dang Q L, Domingues T F, Gu L, Guo S, Han Q, Kattge J, Kubiske M, Manter D, Merilo E, Midgley G F, Porte A, Scales J C, Tissue D, Turnbull T, Warren C, Wohlfahrt G, Wood F I and Wullschleger S D 2014 A Global Data Set of Leaf Photosynthetic Rates, Leaf N and P, and Specific Leaf Area. ORNL DAAC, Oak Ridge, Tennessee, USA. https://doi.org/10.3334/ORNLDAAC/1224

Wang H, Harrison S P, Prentice I C, Yang Y, Bai F, Togashi H F, Wang M, Zhou S and Ni J 2018 The China Plant Trait Database: toward a comprehensive regional compilation of functional traits for land plants Ecology 99 500

Xu H, Wang H, Prentice I C, Harrison S P, Wang G and Sun X 2021 Predictability of leaf traits with climate and elevation: a case study in Gongga Mountain, China Tree Physiol. 41 1336–52

