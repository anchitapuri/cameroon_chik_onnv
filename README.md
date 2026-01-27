# cameroon_chik_onnv
Spatial Modelling of CHIK and ONNV in Cameroon

nrow(original data) = 6336
nrow(after dropping NA in CHIK / ONNV / MAY titre cols) = 5407 (Multisero model ran on this) 

Multisero model Labels: 
  - ONNV_pos: 
    - 0 = 4109
    - 1 = 1298
  - CHIK_pos: 
    - 0 = 5379
    - 1 = 28



# Pre processing data 
District level geometery from: (add links) 
-  Caedistricts179_region.shp (179 districts, 183 geometeries (MANOKA == 5 geometries) 
-  cmr_admin3.shp (360 districts)
  - only these district geometeries used from  cmr_admin3.shp
  - 'belabo','belel', 'evodoula', 'fotokol','gazawa', 'goulfey', 'nguelemendouka', 'njombe-penja', 'oku'
  - (these districts were in not Caedistricts179_region but were in the data) 

# Data 
208 unique districts 

# Mismatch between Shapefiles and Data: 
40 districts (875 rows in the data) 


```text
Mismatche mapping used: (second name is as it appears in the shapefile)
    district_lower == "njombe penja" ~ "njombe-penja",
    district_lower == "tchollire" ~ "tcholire",
    district_lower == "cite verte" ~ "cite vert",
    district_lower == "malentouen" ~ "malantouen",
    district_lower == "nkongsamba" ~ "nkonsamba",
    district_lower == "guidiguis" ~ "guidiguise",
    district_lower == "maroua 3" ~ "maroua rural",
    district_lower == "maroua 1" ~ "maroua urbain",
    district_lower == "kumba-north" ~ "kumba",
    district_lower == "bamenda 3" ~ "bamenda",
    district_lower == "garoua 1" ~ "garoua i",
    district_lower == "ngaoundal" ~ "ngaoundere rural",
    district_lower == "garoua urbain" ~ "garoua boulai",
    district_lower == "eyumodjock" ~ "eyumojock",
    district_lower == "garoua 2" ~ "garoua ii",
    district_lower == "ndikinimeki" ~ "ndikinimiki",
    district_lower == "mbandjock" ~ "mbanjock",
    district_lower == 'bandjoun' ~ "banjoun", 
    district_lower == 'bangangte' ~ "bangante",
    district_lower == 'bangourain' ~ "bangorain"
```


After using both shapefiles + mapping, these districts are still not matched (ie do have geometry information for these, 
and so rows with these districts are removed from downstream analysis) 

```text
 district_lower  n (= rows in the data  with these districts) 
1            boko 39
2            dang 33
3          mozogo 31
4        maroua 2 24
5          bangue 22
6          japoma 15
7            odza 15
8             abo 13
9      nkolbisson 13
10       mvog-ada  6
11           <NA>  3
12   garoua rural  2
13         abeche  1
14        biltine  1
```
    
# After preprocessing 
- Total rows in meta_data: 6336 
- Total rows after merge: 6336 
- Rows with geometry: 6118
  = 218 rows without geometry


# Remove NA from CHIK, ONNV and MAY + Remove duplicate samples 
Total remaining rows: 5185

(Multisero model + Spatial analysis on these 5185 rows with geometry info + valid CHIK / ONNV / MAY titres) 
