########################################################
#   Author: Isabel De Ramos                            #
#   Date Created: 14 April 2022                        #
#   Function: Life Expectancy CBSA Data Preparation    #
########################################################

### loading libraries
#library(dplyr)
library(tidyverse)
#library(stringr)
#library(purrr)

### Get mortality file names
# NOTE: individual level mortality files obtained from NCHS were cleaned and exported as .rdata's beforehand 
df_mort_files = tibble(file = list.files(path='../../NCHS Mortality Data/Clean/Raw', pattern = "mort[1|2]", full.names = TRUE)) %>% 
  mutate(year = str_sub(file,41,44) %>% parse_integer() ) %>% 
  arrange(year) 

#################### PREP WORK #################### 
# (1) create function that formats fips 
# (2) prepare population denominators 

#~~~~ (1) format fips ~~~~#
# create function that pads all fips (adds back leading zeroes)
format_5dig_fips <- function(fp) {
  ifelse(str_length(fp)==4,
         str_pad(fp, width=5, side="left", pad="0"),
         fp)
}
format_5digfips_v <- Vectorize(format_5dig_fips)

#~~~~ (2) preparing population denominators ~~~~#
# pop_county_asrh
# years 1990-2019 available 
# split by age group, gender, and race
# race is 4 categories = H, NHB, NHW, NHO

# bridged_pop_county_asrh
# years 1990-2019 available 
# split by age group (0,5,10,15...), gender, and race
# race is 5 categories = H, NHB, NHW, NHAIAN, NHAPI

# bridged_pop_county_asrhLT 
# years 1990-2019 available
# split by age group (0,1,5,10...), gender, and race
# race is 5 categories = H, NHB, NHW, NHAIAN, NHAPI

# using bridged_pop_county_asrhLT as more appropriate to build life tables with age groups 0,1,5,10...
pop_county <- read.csv("../../Population Data/Clean/bridged_pop_county_asrhLT.csv", header=TRUE, sep = ',') %>% 
  as_tibble() %>% 
  mutate(year=as.character(year),
         fips=as.character(fips),
         age_5yr_group=as.character(age_5yr_group),
         fips=format_5digfips_v(fips)) %>% 
  rename(`sex`=`male`,
         `race`=`hispanic_bridged`)


#################### CREATING FUNCTIONS TO CLEAN AND PREPARE NCHS MORTALITY AND POPULATION DATA #################### 
# (1) clean_nchs 
# (2) clean_nchs1
# (3) clean_popdenoms
# (4) clean_popdenoms1


#~~~~ (1) clean_nchs ~~~~#
# clean_nchs extracts variables of interest from mortality data - year, fips, age, sex, and race
clean_nchs <- function(file) {
  #file <- df_mort_files %>% filter(year%in%"1991") %>% pull(file) 
  load(file) 
  # after loading, object loaded into global environment is mortXXXX
  
  # NCHS DATASET - find counts per age group by race
  mort_dta <- mort_dta %>% 
    transmute(year=as.character(death_year),
              fips=res_fips_effective18,
              age=as.numeric(age_red),
              age_5yr_group=case_when(
                age == 0 ~ '0',
                age < 5 ~ '1',
                between(age, 5,84)~as.character(floor(age/5)*5),
                TRUE~"85"),
              sex=male,
              race=hispanic_bridged) %>%
    group_by(fips, year, age_5yr_group, sex, race) %>% 
    count(age_5yr_group) %>%
    rename(count=n) %>% 
    arrange(fips, year, age_5yr_group, sex, race) %>%
    ungroup()
  print(mort_dta) %>% 
    return()
}

#~~~~ (2) clean_nchs1 ~~~~#
clean_nchs1 <- function(nchs_dta_tmp) {
  #nchs_dta_tmp <- nchs
  
  # (1) filters out non-continental FIPS 
  nonUS_FIPS=c("^(02)","^(15)","^(60)","^(66)","^(69)","^(72)","^(78)") 
  # Alaska, American Samoa, Guam, Northern Mariana Islands, Puerto Rico, Virgin Islands
  nchs_dta_tmp <- nchs_dta_tmp %>% filter(!grepl(paste(nonUS_FIPS, collapse="|"), fips))
  
  # (2) filters out foreign residents (fips="00000", "00999)
  nchs_dta_tmp <- nchs_dta_tmp %>% filter(!fips%in%c("00000", "00999"))
  
  # (3) filters out missing fips (checked in RAW data, missing fips were res_fips from US territories (PR, GUA, VI , etc.))
  nchs_dta_tmp <- nchs_dta_tmp %>% filter(!is.na(fips)==TRUE)
  
  # (4) fixes FIPS that were renamed over time 
  nchs_dta_tmp <- nchs_dta_tmp %>% 
    mutate(fips=ifelse(fips=="12025", "12086", fips)) %>% # fixing 12025
    mutate(fips=ifelse(fips=="51560", "51005", fips)) %>% # fixing 51560
    mutate(fips=ifelse(fips=="51780", "51083", fips)) %>% # fixing 51780
    # special case: 30113 
    # pop_denoms %>% filter(fips%in%"30031") %>% count(fips, wt=pop_denom) >> population of 1863635
    # pop_denoms %>% filter(fips%in%"30067") %>% count(fips, wt=pop_denom) >> population of 332532
    # since FIPS 30031 has largest population, 30113, 30067 and 30031 will all merge into 30031 
    mutate(fips=ifelse(fips%in%c("30113", "30067"), "30031", fips)) %>%
    arrange(year, age_5yr_group, race) %>% 
    group_by(fips, year, age_5yr_group, sex, race) %>%
    summarise(count=sum(count)) %>% 
    ungroup()
  print(nchs_dta_tmp) %>% 
    return()
}


#~~~~ (3) clean_popdenoms ~~~~#
# clean_popdenoms extracts variables of interest - year, fips, age_5yr_group, sex, race, and pop_county
clean_popdenoms <- function(file) {
  #file <- df_mort_files %>% filter(year==2000) %>% pull(file) 
  
  # POP DENOMINATORS - isolate all fips / all groups for given year 
  pop_county_tmp <- pop_county %>% filter(year%in%substr(file, 18,21)) %>%
    group_by(fips, year, age_5yr_group, sex, race) %>% 
    summarise(pop_denom=sum(pop_county)) %>%
    ungroup()
  print(pop_county_tmp) %>%
    return()
}


#~~~~ (4) clean_popdenoms1 ~~~~#
## clean_popdenoms1  fixes problematic FIPS  
## (a) fixes FIPS that were renamed over time 
clean_popdenoms1 <- function(pop_county_tmp) {
  #pop_county_tmp <- pop_denoms
  
  # (a) fixes FIPS that were renamed over time 
  pop_county_tmp <- pop_county_tmp %>% mutate(fips=ifelse(fips=="46113", "46102", fips)) %>% # fixing 46102
    mutate(fips=ifelse(fips=="51560", "51005", fips)) %>%  # fixing 51560
    mutate(fips=ifelse(fips%in%c("30113", "30067"), "30031", fips)) %>% # fixing 30113 (see explanation above in clean_nchs1)
    mutate(sex=as.character(sex)) %>% 
    arrange(year, age_5yr_group, race) %>% 
    group_by(fips, year, age_5yr_group, sex, race) %>%
    summarise(pop_denom=sum(pop_denom)) %>% 
    ungroup()
  print(pop_county_tmp) %>%
    return()
}



#################### CREATING FUNCTIONS FOR STRATIFYING FUNCTIONS #################### 
# (1) crosswalk_cbsa_cz

#~~~~ (1) crosswalk_cbsa_cz ~~~~#
# crosswalk_cbsa_cz joins NCHS mortality data to population denoms, 
# adds CBSA's to FIPS
# adds CZ's to FIPS 

crosswalk_cbsa_cz <- function(nchs_dta_tmp, pop_county_tmp) {
  # join mortality (NCHS) to pops (POP DENOMS)
  county_dta_tmp <- nchs_dta_tmp %>% full_join(pop_county_tmp, by=c("fips", "year", "age_5yr_group", "sex", "race")) %>%
    arrange(fips, year, age_5yr_group, sex, race)
  # creating template to see full combination of fips / age groups / sexes/ races
  template <- expand_grid(fips=unique(county_dta_tmp$fips), 
                          age_5yr_group=unique(county_dta_tmp$age_5yr_group),
                          sex=unique(county_dta_tmp$sex),
                          race=unique(county_dta_tmp$race)) %>% 
    arrange(fips, age_5yr_group, sex, race)
  # replacing NA counts (deaths) with 0 
  county_dta_tmp <- template %>% full_join(county_dta_tmp, by=c("fips", "age_5yr_group", "sex", "race")) %>%
    fill(year) %>%
    mutate(count=replace_na(count,0)) %>% 
    arrange(fips, year, age_5yr_group, sex, race)
  # excluding NA data (caused by faulty fips)
  county_dta_tmp <- county_dta_tmp %>% filter(!is.na(year)==TRUE)
  # excluding fips 13999 (county of < 100k or Georgia HIV death)
      # # NOTE: For the years 1988 through 1991, if there were three or fewer deaths for a given Georgia county of residence 
      # (of deaths occurring in Georgia) with HIV infection (ICD codes *042-*044, 795.8) 
      # cited as a cause-of-death (underlying or non-underlying cause), these records were assigned a "missing" place of 
      # residence code (location code (FIPS code 13999).
      # These deaths do not appear in county death rates, but these deaths are included in the state and national death rates.
  county_dta_tmp <- county_dta_tmp %>% filter(!fips%in%"13999")
  # ADDING CBSA'S 0918 defintion
  # xwalk_cbsa_msa <- read.csv('../Crosswalks/US/Clean/fips_msa_xwalk.csv',
  #                            sep = ',', header =  T, colClasses = "character") %>% as_tibble() %>%
  #   select(fips, cbsa0918)
  # county_dta_tmp <- county_dta_tmp %>% left_join(xwalk_cbsa_msa, by="fips")
  # ADDING CBSA'S 0418 defintion
  xwalk_cbsa_msa <- read.csv('../Crosswalks/US/Clean/fips_msa_complex_xwalk.csv',  
                             sep = ',', header =  T, colClasses = "character") %>% as_tibble() %>% 
    select(fips, cbsa0418)
  county_dta_tmp <- county_dta_tmp %>% left_join(xwalk_cbsa_msa, by="fips")
  # ADDING CZ'S
  xwalk_cz_msa <- read.csv('../Crosswalks/US/Clean/county_FIPS_CZ_crosswalk.csv',  
                             sep = ',', header =  T, colClasses = "character") %>% as_tibble() %>% 
    select(fips, cz) %>% mutate(fips=format_5digfips_v(fips))
  county_dta_tmp <- county_dta_tmp %>% left_join(xwalk_cz_msa, by="fips")
  print(county_dta_tmp) %>%
    return()
}




#################### EXECUTING FUNCTIONS #################### 


### NEED YEARS 1990-2019
mort_files_tmp = df_mort_files %>% filter(between(year, 1990, 2019)) %>% pull(file)

# clean mortality data
nchs <- map_dfr(mort_files_tmp, ~clean_nchs(.x))
nchs <- clean_nchs1(nchs) 
#save(nchs, file='nchs.rdata')

# prepare population denominators
pop_denoms <- map_dfr(mort_files_tmp, ~clean_popdenoms(.x))
pop_denoms <- clean_popdenoms1(pop_denoms)
#save(pop_denoms, file='pop_denoms.rdata')

# load('../../NCHS Mortality Data/Clean/nchs.rdata')
# load('../../Population Data/Clean/pop_denoms.rdata')

# create master file by joining mortality and population denominators, and adding CBSA
master_dta_cbsa <- crosswalk_cbsa_cz(nchs, pop_denoms)
save(master_dta_cbsa, file='1990_2019_cbsa_mortality.rdata')



