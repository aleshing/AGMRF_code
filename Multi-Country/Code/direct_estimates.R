library(SUMMER)
library(readstata13)
library(tidyverse)
library(rgdal)
library(spdep)
library(sf)

# Time Period: 2010-2014
# Countries:
# Burundi (2016, BUBR70DT, BUGE71FL), Ethiopia (2016, ETBR71DT, ETGE71FL), 
# Kenya (2014, KEBR72DT, KEGE71FL), Rwanda (2015, RWBR70DT, RWGE72FL), 
# Tanzania (2015, TZBR7ADT, TZGE7AFL), Uganda (2016, UGBR7BDT, UGGE7AFL)

#### Get direct estimates ####
get_direct_country <- function(file_names, country, survey_year, 
                               stratifying_vars = "v023", admin_info){
  print(paste("Getting direct estimates for", country, survey_year))
  if(admin_info == "DHSREGNA"){
    spatial_info <- 
      readOGR(file_names[2])@data %>%
      select(c(DHSCLUST, DHSREGNA)) %>%
      mutate(v001 = DHSCLUST, admin1 = DHSREGNA) %>% 
      select(v001, admin1)
  }
  else if(admin_info == "ADM1NAME"){
    spatial_info <- 
      readOGR(file_names[2])@data %>%
      select(c(DHSCLUST, ADM1NAME)) %>%
      mutate(v001 = DHSCLUST, admin1 = ADM1NAME) %>% 
      select(v001, admin1)
  }
  else if(admin_info == "Uganda"){
    spatial_info <- 
      readOGR(file_names[2])
    admin1_map <-
      readOGR("../Data/Maps/Uganda/sdr_subnational_boundaries_2022-06-01/shps")
    temp <- sp::over(spatial_info, admin1_map)
    spatial_info <- spatial_info@data %>%
      mutate(v001 = DHSCLUST, 
             admin1 = temp$DHSREGEN) 
    district_to_admin1 <- data.frame(ADM1NAME = c("MPIGI", "WAKISO", 
                                                  "KYANKWANZI", "KIBUKU", 
                                                  "BULAMBULI",  "SHEEMA", 
                                                  "MBARARA",  "BUVUMA" ),
                                     admin1 = c("Central", "Central",
                                                "Central","Eastern",
                                                "Eastern", "Western",
                                                "Western", "Central"))
    for(j in 1:nrow(spatial_info)){
      if(is.na(spatial_info$admin1[j])){
        to_replace <- which(district_to_admin1$ADM1NAME == 
                              spatial_info$ADM1NAME[j])
        spatial_info$admin1[j] <- district_to_admin1$admin1[to_replace]
      }
    }
    
    spatial_info <- spatial_info %>% 
      select(v001, admin1)
    
  }
  
  births_temp <- read.dta13(file_names[1], generate.factors = TRUE)
  births_merged <- merge(births_temp, spatial_info, by = "v001")
  if(country != "Ethiopia"){
    births <- getBirths(data = births_merged, 
                        variables = c("caseid", "v001", "v002", 
                                      "v004", "v005", "v021", "v022", 
                                      "v023", "v024", "v025", "admin1", 
                                      "bidx"),
                        strata = stratifying_vars,
                        surveyyear = survey_year, 
                        year.cut = seq(1990, 2020, by = 5))
  }
  else{
    births <- getBirths(data = births_merged, 
                        variables = c("caseid", "v001", "v002", 
                                      "v004", "v005", "v021", "v022", 
                                      "v023", "v024", "v025", "admin1", 
                                      "bidx"),
                        strata = stratifying_vars,
                        surveyyear = survey_year, 
                        cmc.adjust = 92, 
                        year.cut = seq(1990, 2020, by = 5))
  }
  
  
  
  births <- births %>%
    select(c(v001, v002, admin1, time, age, v005, strata, died))
  colnames(births) <- c("clustid", "id", "region", "time", "age", "weights", 
                        "strata", "died")
  
  years <- levels(births$time)
  direct_estimates <- getDirect(births = births, years = years, 
                                regionVar = "region", timeVar = "time", 
                                clusterVar = "~clustid + id", 
                                ageVar = "age", weightsVar = "weights",
                                national.only = FALSE) %>%
    mutate(region = paste(country, region, sep = "-"))
  
  return(direct_estimates)
}

num_countries <- 6
direct_estimates <- vector(mode = "list", length = num_countries)
file_names <- list(c("../Data/BUBR70DT/BUBR70FL.DTA", "../Data/BUGE71FL"),
                   c("../Data/ETBR71DT/ETBR71FL.DTA", "../Data/ETGE71FL"),
                   c("../Data/KEBR72DT/KEBR72FL.DTA", "../Data/KEGE71FL"),
                   c("../Data/RWBR70DT/RWBR70FL.DTA", "../Data/RWGE72FL"),
                   c("../Data/TZBR7ADT/TZBR7AFL.DTA", "../Data/TZGE7AFL"),
                   c("../Data/UGBR7BDT/UGBR7BFL.DTA", "../Data/UGGE7AFL"))
countries <- c("Burundi", "Ethiopia", "Kenya", "Rwanda", "Tanzania", "Uganda")
survey_years <- c(2016, 2016, 2014, 2015, 2015, 2016)
admin_info <- c("DHSREGNA", "DHSREGNA", "DHSREGNA", "ADM1NAME", "ADM1NAME", 
                "Uganda")

for(i in 1:num_countries){
  direct_estimates[[i]] <- get_direct_country(file_names[[i]], countries[i], 
                                              survey_years[i], 
                                              stratifying_vars = "v023", 
                                              admin_info = admin_info[i])
}
names(direct_estimates) <- countries
  
#### Get adjacency info ####
geos <- vector(mode = "list", length = num_countries)
sf::sf_use_s2(FALSE)
for(i in 1:num_countries){
  if(countries[[i]] != "Tanzania"){
    mapfile <- 
      paste0("../Data/Maps/", countries[i],
             "/sdr_subnational_boundaries_2022-06-01/shps/sdr_subnational_boundaries.shp")
  }
  else{
    mapfile <- 
      paste0("../Data/Maps/", countries[i],
             "/sdr_subnational_boundaries_2022-06-01/shps/sdr_subnational_boundaries2.shp")
  }
  geos[[i]] <- st_read(mapfile) %>%
    mutate(DHSREGEN = paste(countries[[i]], DHSREGEN, sep = "-"))
  
}
names(geos) <- countries

geo_combined <- bind_rows(geos) %>% 
  arrange(DHSREGEN)
amat_combined <- nb2mat(poly2nb(geo_combined), style = "B")
colnames(amat_combined) <- rownames(amat_combined) <- 
  geo_combined$DHSREGEN

# Fix border issues
# Ethiopia Kenya
amat_combined["Kenya-North Eastern", "Ethiopia-Somali"] <- 
  amat_combined["Ethiopia-Somali", "Kenya-North Eastern"] <- 1
amat_combined["Kenya-Eastern", "Ethiopia-Somali"] <- 
  amat_combined["Ethiopia-Somali", "Kenya-Eastern"] <- 1
amat_combined["Kenya-Eastern", "Ethiopia-Oromiya"] <- 
  amat_combined["Ethiopia-Oromiya", "Kenya-Eastern"] <- 1
amat_combined["Kenya-Eastern", "Ethiopia-SNNP"] <- 
  amat_combined["Ethiopia-SNNP", "Kenya-Eastern"] <- 1
amat_combined["Kenya-Eastern", "Ethiopia-SNNP"] <- 
  amat_combined["Ethiopia-SNNP", "Kenya-Eastern"] <- 1
amat_combined["Kenya-Rift Valley", "Ethiopia-SNNP"] <- 
  amat_combined["Ethiopia-SNNP", "Kenya-Rift Valley"] <- 1
# Kenya Uganda
amat_combined["Kenya-Rift Valley", "Uganda-Northern"] <- 
  amat_combined["Uganda-Northern", "Kenya-Rift Valley"] <- 1
amat_combined["Kenya-Rift Valley", "Uganda-Eastern"] <- 
  amat_combined["Uganda-Eastern", "Kenya-Rift Valley"] <- 1
amat_combined["Kenya-Nyanza", "Uganda-Eastern"] <- 
  amat_combined["Uganda-Eastern", "Kenya-Nyanza"] <- 1
# Rwanda Burundi
amat_combined["Rwanda-East", "Burundi-Kirundo"] <- 
  amat_combined["Burundi-Kirundo", "Rwanda-East"] <- 1
amat_combined["Rwanda-East", "Burundi-Muyinga"] <- 
  amat_combined["Burundi-Muyinga", "Rwanda-East"] <- 1
amat_combined["Rwanda-South", "Burundi-Kirundo"] <- 
  amat_combined["Burundi-Kirundo", "Rwanda-South"] <- 1
amat_combined["Rwanda-South", "Burundi-Ngozi"] <- 
  amat_combined["Burundi-Ngozi", "Rwanda-South"] <- 1
amat_combined["Rwanda-South", "Burundi-Kayanza"] <- 
  amat_combined["Burundi-Kayanza", "Rwanda-South"] <- 1
amat_combined["Rwanda-South", "Burundi-Cibitoke"] <- 
  amat_combined["Burundi-Cibitoke", "Rwanda-South"] <- 1
amat_combined["Rwanda-West", "Burundi-Cibitoke"] <- 
  amat_combined["Burundi-Cibitoke", "Rwanda-West"] <- 1

country_admin_table <- data.frame(Country = geo_combined$CNTRYNAMEE,
                                  Admin1 = geo_combined$DHSREGEN)

#### Save everything ####
# Remove islands in Tanzania
islands <- c("Tanzania-Kaskazini Pemba", "Tanzania-Kusini Pemba",
             "Tanzania-Mjini Magharibi", "Tanzania-Kaskazini Unguja", 
             "Tanzania-Kusini Unguja")
direct_estimates[["Tanzania"]] <- direct_estimates[["Tanzania"]] %>%
  filter(!(region %in% islands))
to_keep <- which(!(colnames(amat_combined) %in%  islands))
amat_combined <- amat_combined[to_keep, to_keep]
country_admin_table <- country_admin_table %>%
  filter(!(Admin1  %in% islands))
geo_combined <- geo_combined %>%
  filter(!(DHSREGEN  %in% islands))
# Fix some names
direct_estimates[["Ethiopia"]] <- direct_estimates[["Ethiopia"]] %>%
  mutate(region = recode(region,
                         "Ethiopia-Afar" = "Ethiopia-Affar",
                         "Ethiopia-Addis Adaba" = "Ethiopia-Addis Abeba",
                         "Ethiopia-SNNPR" = "Ethiopia-SNNP",
                         "Ethiopia-Benishangul" = "Ethiopia-Ben-Gumz",
                         "Ethiopia-Oromia" = "Ethiopia-Oromiya")) %>% 
  arrange(region, years)
direct_estimates[["Rwanda"]] <- direct_estimates[["Rwanda"]] %>%
  mutate(region = recode(region,
                         "Rwanda-Kigali City" = "Rwanda-City of Kigali")) %>% 
  arrange(region, years)

save(direct_estimates, file = "../Data/generated_data/direct_estimates.RData")
save(amat_combined, country_admin_table, geo_combined,
     file = "../Data/generated_data/admin_info.RData")
