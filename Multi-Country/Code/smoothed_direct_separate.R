library(SUMMER)
library(readstata13)
library(dplyr)
library(rgdal)

#### Load helper functions and data ####
load(file = "../Data/generated_data/direct_estimates.RData")
load(file = "../Data/generated_data/admin_info.RData")
source("../../helper_files/smoothDirect_extra.R")

#### Get smoothed direct estimates, bym2 ####
direct_estimates_2010_2014 <- direct_estimates %>% 
  bind_rows() %>% filter(years == "10-14") %>%
  filter(!grepl("All", region))
amat <- amat_combined

countries <- unique(country_admin_table$Country)
smoothed_direct_separate <- vector(mode = "list", length = length(countries))
out_separate_list <- vector(mode = "list", length = length(countries))
for(i in 1:length(countries)){
  inds <- which(country_admin_table$Country == countries[i])
  temp_amat <- amat[inds, inds]
  temp_direct <- direct_estimates_2010_2014[inds, ]
  temp_smoothed <- smoothDirect_extra(data = temp_direct,
                                      Amat = temp_amat, 
                                      time.model = NULL,
                                      year_label = "10-14",
                                      year_range = c(2010, 2014), 
                                      is.yearly = FALSE, m = 5,
                                      verbose = FALSE)
  smoothed_direct_separate[[i]] <- temp_smoothed
  temp_out <- getSmoothed(temp_smoothed)
  out_separate_list[[i]] <- temp_out
}

out_separate <- out_separate_list %>% 
  bind_rows()

save(out_separate, smoothed_direct_separate, 
     file = "../Results/smoothed_direct_separate.RData")
