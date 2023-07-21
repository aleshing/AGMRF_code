library(SUMMER)
library(readstata13)
library(dplyr)

#### Load helper functions ####
get_meta_analysis_est <- function(direct_estimates, years){
    means <- data.frame(years = years)
    means$est <- 0
    for(i in 1:length(years)){
        year <- years[i]
        temp <- direct_estimates %>% 
            filter(years == year & is.na(mean) == FALSE)
        means$est[i] <- SUMMER::expit(sum(temp$logit.est * temp$logit.prec) / 
                                          sum(temp$logit.prec))
    }
    means <- means %>% filter(is.nan(est) == FALSE)
    return(means)
}

#### Get births data ####
num_surveys <- 6
filenames <- c("../Data/RWBR21DT/RWBR21FL.dta",
               "../Data/RWBR41DT/RWBR41FL.DTA",
               "../Data/RWBR53DT/RWBR53FL.DTA",
               "../Data/RWBR5ADT/RWBR5AFL.DTA",
               "../Data/RWBR61DT/RWBR61FL.DTA",
               "../Data/RWBR70DT/RWBR70FL.DTA")
survey_years <- c(1992, 2000, 2005, 2008, 2010, 2015)
births_list <- vector(mode = "list", length = num_surveys)
for(i in 1:num_surveys){
    births_temp <- read.dta13(filenames[i], generate.factors = TRUE)
    if(i == num_surveys){
        births_list[[i]] <- getBirths(data = births_temp, strata = c("v023"),
                                      surveyyear = survey_years[i], 
                                      year.cut = seq(1985, 2020, by = 1))
    }
    else{
        births_list[[i]] <- getBirths(data = births_temp, 
                                      strata = c("v024", "v025"),
                                      surveyyear = survey_years[i], 
                                      year.cut = seq(1985, 2020, by = 1))
    }
    
    births_list[[i]] <- births_list[[i]][, c("v001", "v002", "v024", "time", 
                                             "age", "v005", "strata", "died")]
    colnames(births_list[[i]]) <- c("clustid", "id", "region", "time", "age", 
                                    "weights", "strata", "died")
}
names(births_list) <- survey_years

#### Get direct estimates ####
years <- levels(births_list[[1]]$time)
direct_estimates <- getDirectList(births = births_list, years = years, 
                                  regionVar = "region", timeVar = "time", 
                                  clusterVar = "~clustid + id", 
                                  ageVar = "age", weightsVar = "weights",
                                  national.only = TRUE)
head(direct_estimates)

#### Limit direct estimates to X years prior to survey ####
direct_estimates_all <- direct_estimates
limit_years <- 15
direct_estimates <- direct_estimates %>%
  filter(as.numeric(surveyYears) - as.numeric(years) <= limit_years)

#### Get meta analysis estimates ####
meta_analysis_estimates_all <- get_meta_analysis_est(direct_estimates_all, 
                                                     years)
meta_analysis_estimates <- get_meta_analysis_est(direct_estimates, years)

#### Get IGME estimates ####
# igme <- read.csv("../Data/IGME/IGME_U5MR_20191007.csv") 
# igme <- igme %>% filter(REF_AREA == "Rwanda", 
#                         INDICATOR == "Under-five mortality rate",
#                         SEX == "Total",
#                         SERIES_NAME == "UN IGME estimate 2019") 
# write.csv(igme, file = "../Data/IGME/IGME_U5MR_20191007_Rwanda.csv")
igme <- read.csv("../Data/IGME/IGME_U5MR_20191007_Rwanda.csv") 
igme <- igme %>%
    mutate(years = as.numeric(substr(as.character(TIME_PERIOD), 1, 4)),
           est = OBS_VALUE / 1000, lower = LOWER_BOUND / 1000, 
           upper = UPPER_BOUND / 1000) %>%
    filter(years >= 1985)

save(births_list, direct_estimates, meta_analysis_estimates, igme, 
     file = "../Data/generated_data/direct_estimates.RData")
