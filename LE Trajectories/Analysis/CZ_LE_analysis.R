## PENDINGS for R1 EPID:
### think of doing resampling for the significance part?


#******************************************************#
#   Author: Isabel De Ramos                            #
#   Date Created: 14 April 2022                        #
#   Function: Life Expectancy CBSA Data Preparation    #
#******************************************************#

                    #******************************************************************************#
                    #       Enable document outline pane for better overview of R script           #
                    #                 and quicker navigation to any section                        #
                    #******************************************************************************#

# INSTALL LIBRARIES ----
library(tidyverse)
library(classInt)
library(grid)
library(gridExtra)
library(scales)
library(multcomp)
library(RColorBrewer)
library(scales)
library(ggpubr)
#install.packages("devtools")
#library(devtools)
#install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
#install.packages("remotes")
#remotes::install_github("timriffe/DemoTools")
#install.packages("DemoTools")
library(DemoTools)
library(RColorBrewer)
#remotes::install_github("coolbutuseless/ggpattern")
#install.packages("ggpattern")
library(ggpattern)
library(moments)
#install.packages("biscale")
library(biscale)
library(sf)
library(cowplot)
#devtools::install_github("ropenscilabs/rnaturalearth")
library(rnaturalearth)
# install.packages("rnaturalearthdata")
library(rnaturalearthdata)
library(sp)
# install.packages("spdep")
library(spdep)
library(rgeoda)
# install.packages("pals")
library(pals)



# PREP WORK ----
# (1) ggplot helpers
# (2) helper functions
#         (2.A) Camarda lifetable()
#         (2.B) le_lv()
#         (2.C) get_bivariate()
# (3) loading in .rdata
# (4) create 5-year pooled periods and factors
# (5) U.S. shapefiles 

#~~~~ (1) ggplot helpers ~~~~#
### prep for generating plots ###
select<-dplyr::select
fontsize<-16
isabel_theme <-theme_bw()+
  theme(axis.text=element_text(color="black", size=fontsize),
        axis.title=element_text(color="black", size=fontsize, face="bold"),
        plot.title=element_text(color="black", size=fontsize, face="bold", hjust = 0.5),
        strip.background=element_blank(),
        strip.text=element_text(color="black", size=fontsize, face="bold"),
        legend.text=element_text(color="black", size=fontsize),
        legend.title=element_text(color="black", size=fontsize, face="bold"))


map_theme <- isabel_theme+
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_rect(fill = "lightblue1"),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background = element_rect(fill = "lightblue1"),
        plot.margin = unit(c(0,0,0,0), "cm"),
        plot.tag = element_text(face="bold", size=40),
        legend.text = element_text(size=24),
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent')) #transparent legend panel

moran_colors<-c("lightgray", "dodgerblue4", "firebrick", "black", "white")

# bivariat3 palette
cols<-brewer.divdiv(n=9)
# change middle color for a darker gray
cols[5]<-"gray50"
cols[6]<- "skyblue3"

get_legend<-function(plot){
  grobs<-ggplotGrob(plot)$grobs
  legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]  
  return(legend)
}

#~~~~ (2) loading in helper functions ~~~~#

## (2.A) CAMARDA'S FUNCTIONS FOR LIFETABLES, LIFE EXPECTANCY 
# see https://sites.google.com/site/carlogiovannicamarda/r-stuff/life-expectancy-confidence-interval 
## lifetable() calculates life expectancy
## CIex() calculates 95% CI 
source('LifeTableFUN.R')

## (2.B) LE_LV FUNCTION 
le_lv<-function(.x, .y, age_num, sex="B"){
  ## .x -> the data with one race and CBSA
  # example: .x<- dta %>% filter(gender=="Men", cbsa=="10100", year5%in%'2015-2019')
  ## .y -> a dataframe wtih whichever CBSA or race category you are looking at
  # using lifetable() function from Camarda to create lifetable 
  lt <- lifetable(x=.x$age_5yr_group, Nx=.x$pop_denom, Dx=.x$count, sex=sex, ax=NULL) 
  # using CIex function from Camarda to obtain confidence intervals 
  ci_dta <- CIex(x=.x$age_5yr_group, Nx=.x$pop_denom, Dx=.x$count, sex=sex, ax=NULL, which.x=age_num, ns=1000, level=0.95)
  # adding variables needed to calculate lifespan variation 
  lt <- lt %>% mutate(xbar=NA_real_,
                      noname=NA_real_,
                      v=NA_real_,
                      sd=NA_real_)
  lt <- lt %>% mutate(xbar=x+ax,
                      noname=case_when(
                        x==0 ~ dx/lx*(xbar-ex)^2,
                        x!=0 ~ {
                          lx_forlv <- lt %>% filter(x==0) %>% pull(lx)
                          ex_forlv <- lt %>% filter(x==0) %>% pull(ex)
                          dx/lx_forlv*(xbar-ex_forlv)^2}),
                      v=rev(cumsum(rev(noname))),
                      sd=sqrt(v))
  # extracting LE
  le <- lt %>% filter(x==age_num) %>% pull(ex)
  # extracting 95% CI
  ci <- ci_dta %>% pluck("CIex") %>% unname() 
  # extracting LV 
  lv_sd <- lt %>% filter(x==age_num) %>% pull(sd) # this is standard deviation 
  lv_cv <- lv_sd/le # this is CoV = coefficient of variance (standard deviation / LE)
  # building df that summarizes results
  data.frame(le=le,
             lci=ci[1],
             uci=ci[2],
             lv_sd=lv_sd,
             lv_cv=lv_cv)
}

## (2.C) GET_BIVARIATE FUNCTION 

# creates bi_class by CBSA/CZ, gender, year3/year5
# NOTE: uses first period (1990-1992 for year, or 3 1990-1994 for year5) as baseline ref

get_bivarite <- function(.x, .y, gender_mw="Men", yr_35="year3", cbsa_cz="cbsa") { 
  
  dta <- results_cbsa_cz %>% filter(gender%in%gender_mw, year_type%in%yr_35, type%in%cbsa_cz)
  
  if(yr_35=="year3") {
    ## CATEGORIZATION OF LE'S, USING BASELINE (1990-1992) LE USING TERTILES
    le_bin <- dta %>% pull(le) %>% sort() %>% 
      quantile(probs=seq(0,1,1/3))
    
    ## CATEGORIZATION OF GAINS/LOSSES USING TERTILES
    dif_bin <- dta %>% select(year, id, gender, le) %>%
      group_by(year, id, gender) %>% 
      spread(year, le) %>% 
      mutate(abs_dif=`2017-2019`-`1990-1992`,
             rel_dif=`2017-2019`/`1990-1992`) %>% 
      ungroup() %>% 
      gather("year", "le", 3:12) %>% 
      pull(abs_dif) %>% 
      sort() %>% 
      quantile(probs=seq(0,1,1/3))
    
    df <- dta %>% select(year, id, gender, le) %>%
      group_by(year, id, gender) %>% 
      spread(year, le) %>% 
      mutate(abs_dif=`2017-2019`-`1990-1992`,
             rel_dif=`2017-2019`/`1990-1992`) %>% 
      ungroup() %>% 
      gather("year", "le", 3:12) %>% 
      filter(year%in%"1990-1992") %>% # using 1990-1992 as baseline
      mutate(cat=case_when(
        year%in%"1990-1992" & le>=le_bin[3] & abs_dif>dif_bin[3] ~ 9,
        year%in%"1990-1992" & le>=le_bin[3] & between(abs_dif, dif_bin[2], dif_bin[3]) ~ 8,
        year%in%"1990-1992" & le>=le_bin[3] &  abs_dif< dif_bin[2] ~ 7,
        year%in%"1990-1992" & between(le, le_bin[2], le_bin[3]) & abs_dif>dif_bin[3] ~ 6,
        year%in%"1990-1992" & between(le, le_bin[2], le_bin[3]) & between(abs_dif, dif_bin[2], dif_bin[3]) ~ 5,
        year%in%"1990-1992" & between(le, le_bin[2], le_bin[3]) & abs_dif< dif_bin[2] ~ 4,
        year%in%"1990-1992" & le<le_bin[2] & abs_dif>dif_bin[3] ~ 3, 
        year%in%"1990-1992" & le<le_bin[2] & between(abs_dif, dif_bin[2], dif_bin[3]) ~ 2,
        year%in%"1990-1992" & le<le_bin[2] & abs_dif< dif_bin[2] ~ 1)) 
    ## breakdown: abs_dif = difference from most recent period (2017-2019) to oldest (1990-1992)
    # 9 = high LE, most increase since 1990-1992
    # 8 = high LE, moderate increase since 1990-1992
    # 7 = high LE, most decrease since 1990-1992
    # 6 = moderate LE, most increase since 1990-1992
    # 5 = moderate LE, moderate increase since 1990-1992
    # 4 = moderate LE, most decrease since 1990-1992
    # 3 = low LE, most increase since 1990-1992
    # 2 = low LE, moderate increase since 1990-1992
    # 1 = low LE, most decrease since 1990-1992
    
    df <- bi_class(df, x=le, y=abs_dif, style="quantile", dim=3)
    print(df)
    
  } else if (yr_35=="year5") {
    ## CATEGORIZATION OF LE'S, USING BASELINE (1990-1992) LE USING TERTILES
    le_bin <- dta %>% pull(le) %>% sort() %>% 
      quantile(probs=seq(0,1,1/3))
    
    dif_bin <- dta %>% select(year, id, gender, le) %>%
      group_by(year, id, gender) %>% 
      spread(year, le) %>% 
      mutate(abs_dif=`2015-2019`-`1990-1994`,
             rel_dif=`2015-2019`/`1990-1994`) %>% 
      ungroup() %>% 
      gather("year", "le", 3:8) %>% 
      pull(abs_dif) %>% 
      sort() %>% 
      quantile(probs=seq(0,1,1/3))
    
    df <- dta %>% select(year, id, gender, le) %>%
      group_by(year, id, gender) %>% 
      spread(year, le) %>% 
      mutate(abs_dif=`2015-2019`-`1990-1994`,
             rel_dif=`2015-2019`/`1990-1994`) %>% 
      ungroup() %>% 
      gather("year", "le", 3:8) %>% 
      filter(year%in%"1990-1994") %>% # using 1990-1992 as baseline
      mutate(cat=case_when(
        year%in%"1990-1994" & le>=le_bin[3] & abs_dif>dif_bin[3] ~ 9,
        year%in%"1990-1994" & le>=le_bin[3] & between(abs_dif, dif_bin[2], dif_bin[3]) ~ 8,
        year%in%"1990-1994" & le>=le_bin[3] &  abs_dif< dif_bin[2] ~ 7,
        year%in%"1990-1994" & between(le, le_bin[2], le_bin[3]) & abs_dif>dif_bin[3] ~ 6,
        year%in%"1990-1994" & between(le, le_bin[2], le_bin[3]) & between(abs_dif, dif_bin[2], dif_bin[3]) ~ 5,
        year%in%"1990-1994" & between(le, le_bin[2], le_bin[3]) & abs_dif< dif_bin[2] ~ 4,
        year%in%"1990-1994" & le<le_bin[2] & abs_dif>dif_bin[3] ~ 3, 
        year%in%"1990-1994" & le<le_bin[2] & between(abs_dif, dif_bin[2], dif_bin[3]) ~ 2,
        year%in%"1990-1994" & le<le_bin[2] & abs_dif< dif_bin[2] ~ 1)) 
    
    df <- bi_class(df, x=le, y=abs_dif, style="quantile", dim=3)
    print(df)
  }
}


#~~~~ (3) loading in master datafile (see LE_data_prep.R) ~~~~#
load('1990_2019_cbsa_mortality.rdata')

#~~~~ (4) create year5 and year3 variable that pools periods from 2000-2019 ~~~~#
dta <- master_dta_cbsa %>% select(-fips) %>% 
  mutate(gender=as.numeric(sex),
         cbsa=cbsa0418,
         age_5yr_group=as.numeric(age_5yr_group),
         year5=case_when(
           year%in%c(2015:2019) ~ as.character('2015-2019'),
           year%in%c(2010:2014) ~ as.character('2010-2014'),
           year%in%c(2005:2009) ~ as.character('2005-2009'),
           year%in%c(2000:2004) ~ as.character('2000-2004'),
           year%in%c(1995:1999) ~ as.character('1995-1999'),
           year%in%c(1990:1994) ~ as.character('1990-1994')),
         year3=case_when(
           year%in%c(2017:2019) ~ as.character('2017-2019'),
           year%in%c(2014:2016) ~ as.character('2014-2016'),
           year%in%c(2011:2013) ~ as.character('2011-2013'),
           year%in%c(2008:2010) ~ as.character('2008-2010'),
           year%in%c(2005:2007) ~ as.character('2005-2007'),
           year%in%c(2002:2004) ~ as.character('2002-2004'),
           year%in%c(1999:2001) ~ as.character('1999-2001'),
           year%in%c(1996:1998) ~ as.character('1996-1998'),
           year%in%c(1993:1995) ~ as.character('1993-1995'),
           year%in%c(1990:1992) ~ as.character('1990-1992'))) %>%
  # factoring gender to help in ggplot later
  mutate(gender=factor(gender, levels=c(1,0),
                       labels=c("Men", "Women"))) %>% select(-cbsa0418) %>% 
  # dropping NA CZ's (this still keeps CZ's with no CBSA's!)
  filter(!is.na(cz)==TRUE) 



#~~~~ (5) U.S. shapefiles ~~~~#
# temp_shapefile <- tempfile()
# download.file("https://www2.census.gov/geo/tiger/GENZ2018/shp/cb_2018_us_cbsa_500k.zip", temp_shapefile)
# unzip(temp_shapefile, exdir="../Data/CBSA shapefile")
shp_cbsa <- read_sf("Shapefiles/CBSA shapefile/cb_2018_us_cbsa_500k.shp")

# temp_shapefile <- tempfile()
# download.file("https://sites.psu.edu/psucz/files/2018/09/ERS-2010-2fcqzlr.zip", temp_shapefile)
# unzip(temp_shapefile, exdir="../Data/CZ shapefile")
shp_cz <- read_sf("Shapefiles/CZ shapefile/ERS10.shp") %>% mutate(LM_Code=as.character(LM_Code))

# temp_shapefile <- tempfile()
# download.file("https://www2.census.gov/geo/tiger/GENZ2018/shp/cb_2018_us_state_500k.zip", temp_shapefile)
# unzip(temp_shapefile, exdir="../Data/state shapefile")
shp_state <- read_sf("Shapefiles/state shapefile/cb_2018_us_state_500k.shp")

# temp_shapefile <- tempfile()
# download.file("https://www2.census.gov/geo/tiger/GENZ2018/shp/cb_2018_us_region_500k.zip", temp_shapefile)
# unzip(temp_shapefile, exdir="Shapefiles/census region shapefile")
shp_census_region <- read_sf("Shapefiles/census region shapefile/cb_2018_us_region_500k.shp")

# temp_shapefile <- tempfile()
# download.file("https://www2.census.gov/geo/tiger/GENZ2018/shp/cb_2018_us_division_500k.zip", temp_shapefile)
# unzip(temp_shapefile, exdir="Shapefiles/census division shapefile")
shp_census_division <- read_sf("Shapefiles/census division shapefile/cb_2018_us_division_500k.shp")

df_state <- shp_state %>% filter(!GEOID%in%c("02", "15")) # only need continental U.S. 
# df_cbsa <- right_join(shp_cbsa, df, by=c("GEOID"="id"))
# df_cz <- right_join(shp_cz, df, by=c("LM_Code"="id"))
df_mexico <- ne_countries(country='mexico', scale=50, returnclass = "sf") %>% select(geometry)
df_canada <- ne_countries(country='canada', scale=50, returnclass = "sf") %>% select(geometry)



# LE CALCULATIONS ----

# IMPORTANT NOTE: ALL RESULTS BELOW HAVE ALREADY BEEN STORED IN LE_CBSA_RESULTS.RDATA (SKIP TO LINE 395)
run_LE<-F
if (run_LE){
  ## cz, 3 yr pooled ----
  results_czyr3 <-
    dta %>% group_by(year3, cz, age_5yr_group, gender) %>%  
    summarise(count=sum(count, na.rm=T),
              pop_denom=sum(pop_denom, na.rm=T)) %>%
    ungroup() %>%
    filter(gender=="Men") %>%
    # creating age-specific death rates variable, Mx
    mutate(mx=count/pop_denom,
           age_5yr_group=as.numeric(age_5yr_group)) %>%
    arrange(year3, cz, age_5yr_group, gender) %>%
    group_by(year3, cz, gender) %>%
    group_modify(~le_lv(., age_num=0, sex="M")) %>%
    bind_rows(
      dta %>% group_by(year3, cz, age_5yr_group, gender) %>%
        summarise(count=sum(count, na.rm=T),
                  pop_denom=sum(pop_denom, na.rm=T)) %>%
        ungroup() %>%
        filter(gender=="Women") %>%
        # creating age-specific death rates variable, Mx
        mutate(mx=count/pop_denom,
               age_5yr_group=as.numeric(age_5yr_group)) %>%
        arrange(year3, cz, age_5yr_group, gender) %>%
        group_by(year3, cz, gender) %>%
        group_modify(~le_lv(., age_num=0, sex="W"))
    )
  
  ## cz, 5 yr pooled ---- 
  results_czyr5 <-
    dta %>% group_by(year5, cz, age_5yr_group, gender) %>%  
    summarise(count=sum(count, na.rm=T),
              pop_denom=sum(pop_denom, na.rm=T)) %>%
    ungroup() %>%
    filter(gender=="Men") %>%
    # creating age-specific death rates variable, Mx
    mutate(mx=count/pop_denom,
           age_5yr_group=as.numeric(age_5yr_group)) %>%
    arrange(year5, cz, age_5yr_group, gender) %>%
    group_by(year5, cz, gender) %>%
    group_modify(~le_lv(., age_num=0, sex="M")) %>%
    bind_rows(
      dta %>% group_by(year5, cz, age_5yr_group, gender) %>%
        summarise(count=sum(count, na.rm=T),
                  pop_denom=sum(pop_denom, na.rm=T)) %>%
        ungroup() %>%
        filter(gender=="Women") %>%
        # creating age-specific death rates variable, Mx
        mutate(mx=count/pop_denom,
               age_5yr_group=as.numeric(age_5yr_group)) %>%
        arrange(year5, cz, age_5yr_group, gender) %>%
        group_by(year5, cz, gender) %>%
        group_modify(~le_lv(., age_num=0, sex="W"))
    )
  
  ## cbsa, 3 yr pooled ----
  results_cbsayr3 <-
    dta %>% filter(!cbsa%in%"") %>% filter(!is.na(cbsa)==TRUE) %>% 
    group_by(year3, cbsa, age_5yr_group, gender) %>%
    summarise(count=sum(count, na.rm=T),
              pop_denom=sum(pop_denom, na.rm=T)) %>%
    ungroup() %>%
    filter(gender=="Men") %>%
    # creating age-specific death rates variable, Mx
    mutate(mx=count/pop_denom,
           age_5yr_group=as.numeric(age_5yr_group)) %>%
    arrange(year3, cbsa, age_5yr_group, gender) %>%
    group_by(year3, cbsa, gender) %>%
    group_modify(~le_lv(., age_num=0, sex="M")) %>%
    bind_rows(
      dta %>% filter(!cbsa%in%"") %>% filter(!is.na(cbsa)==TRUE) %>% 
        group_by(year3, cbsa, age_5yr_group, gender) %>%
        summarise(count=sum(count, na.rm=T),
                  pop_denom=sum(pop_denom, na.rm=T)) %>%
        ungroup() %>%
        filter(gender=="Women") %>%
        # creating age-specific death rates variable, Mx
        mutate(mx=count/pop_denom,
               age_5yr_group=as.numeric(age_5yr_group)) %>%
        arrange(year3, cbsa, age_5yr_group, gender) %>%
        group_by(year3, cbsa, gender) %>%
        group_modify(~le_lv(., age_num=0, sex="W"))
    )
  
  
  ## cbsa, 5 yr pooled ----
  # FINDING LE FOR EACH CBSA BY SEX, 1990-2019
  results_cbsayr5 <-
    dta %>% filter(!cbsa%in%"") %>% filter(!is.na(cbsa)==TRUE) %>% 
    group_by(year5, cbsa, age_5yr_group, gender) %>%  
    summarise(count=sum(count, na.rm=T),
              pop_denom=sum(pop_denom, na.rm=T)) %>%
    ungroup() %>%
    filter(gender=="Men") %>%
    # creating age-specific death rates variable, Mx
    mutate(mx=count/pop_denom,
           age_5yr_group=as.numeric(age_5yr_group)) %>%
    arrange(year5, cbsa, age_5yr_group, gender) %>%
    group_by(year5, cbsa, gender) %>%
    group_modify(~le_lv(., age_num=0, sex="M")) %>%
    bind_rows(
      dta %>% filter(!cbsa%in%"") %>% filter(!is.na(cbsa)==TRUE) %>% 
        group_by(year5, cbsa, age_5yr_group, gender) %>%
        summarise(count=sum(count, na.rm=T),
                  pop_denom=sum(pop_denom, na.rm=T)) %>%
        ungroup() %>%
        filter(gender=="Women") %>%
        # creating age-specific death rates variable, Mx
        mutate(mx=count/pop_denom,
               age_5yr_group=as.numeric(age_5yr_group)) %>%
        arrange(year5, cbsa, age_5yr_group, gender) %>%
        group_by(year5, cbsa, gender) %>%
        group_modify(~le_lv(., age_num=0, sex="W"))
    )
  
  save(results_czyr3, results_czyr5, results_cbsayr3, results_cbsayr5, file='le_cbsa_results.rdata')
  
} else {
  load('le_cbsa_results.rdata')
}

# EXPLORATORY STUFF ----

# 1990-2019 le and lv trends # 
# LOOKING AT LE TRENDS (1990-2019) FOR MEN AND WOMEN FOR EACH CBSA
plots <- results_cbsayr3 %>% group_by(cbsa) %>% 
  group_map(~{
    #.x<-results_cbsa %>% filter(cbsa==10100)
    title<-paste0(unique(.y$cbsa))
    ggplot(.x, aes(x=year3, y=le, group=gender)) +
      geom_line(aes(color=gender)) +
      labs(title=title)
  })
pdf("test.pdf", width=10, height=7.5)
plots
dev.off()



# 2015-19 le & lv associations #
# MEN LE & LV ASSOCIATION, 2015-2019
ggplot(results_cbsayr5 %>% filter(year5%in%"2015-2019", gender%in%"Men"), aes(lv_cv, le)) +
  geom_point(size=1, color="blue") +
  geom_linerange(aes(ymin=lci, ymax=uci), color="blue")+
  geom_smooth(method=lm, level=0.95, color="black")+
  scale_y_continuous(limits=c(results_cbsayr5 %>% ungroup() %>% select(le) %>% min,
                              results_cbsayr5 %>% ungroup() %>% select(le) %>% max),
                     oob=rescale_none)+
  scale_x_continuous(limits=c(results_cbsayr5 %>% ungroup() %>% select(lv_cv) %>% min,
                              results_cbsayr5 %>% ungroup() %>% select(lv_cv) %>% max),
                     oob=rescale_none)+
  labs(x="Lifespan Variation in Men",
       y="Life Expectancy in Men",
       color="", fill="")+
  isabel_theme+
  theme(legend.position="none", legend.title=element_blank())

# WOMEN LE & LV ASSOCIATION, 2015-2019
ggplot(results_cbsayr5 %>% filter(year5%in%"2015-2019", gender%in%"Women"), aes(lv_cv, le)) +
  geom_point(size=1, color="red") +
  geom_linerange(aes(ymin=lci, ymax=uci), color="red")+
  geom_smooth(method=lm, level=0.95, color="black")+
  scale_y_continuous(limits=c(results_cbsayr5 %>% ungroup() %>% select(le) %>% min,
                              results_cbsayr5 %>% ungroup() %>% select(le) %>% max),
                     oob=rescale_none)+
  scale_x_continuous(limits=c(results_cbsayr5 %>% ungroup() %>% select(lv_cv) %>% min,
                              results_cbsayr5 %>% ungroup() %>% select(lv_cv) %>% max),
                     oob=rescale_none)+
  labs(x="Lifespan Variation in Women",
       y="Life Expectancy in Women",
       color="", fill="")+
  isabel_theme+
  theme(legend.position="none", legend.title=element_blank())




# 1990-2019 5yr men v women, le & lv #

## first, creating LE scatterplot
# spreading LE and 95% CIs by gender 
df <- results_cbsayr5 %>% select(cbsa, gender, le) %>% 
  spread(gender,le) %>% 
  rename(men_le=`Men`,
         women_le=`Women`) %>% 
  full_join(
    results_cbsayr5 %>% select(cbsa, gender, lci) %>% 
      spread(gender,lci) %>% 
      rename(men_lci=`Men`,
             women_lci=`Women`)) %>% 
  full_join(
    results_cbsayr5 %>% select(cbsa, gender, uci) %>% 
      spread(gender, uci) %>% 
      rename(men_uci=`Men`,
             women_uci=`Women`))
# generating scatterplot
ggplot(df, aes(x=men_le, y=women_le)) +
  geom_abline(intercept = 0, slope=1, lty=1)+
  geom_linerange(aes(ymin=women_lci, ymax=women_uci))+
  geom_linerange(aes(xmin=men_lci, xmax=men_uci))+
  geom_point(aes(fill=cbsa), size=1) +
  labs(x="Life Expectancy for Men (95% CI)",
       y="Life Expectancy for Women (95% CI)",
       color="", fill="", shape="")+
  #guides(color=F)+
  facet_wrap(~year5)+
  isabel_theme+
  theme(legend.position="none", legend.title=element_blank())


# computing RSE
results_czyr5 %>% 
  # back calculating SE (could do it in function above)
  mutate(width=uci-le,
         le_se=width/1.96,
         rse=le_se/le) %>% 
  ggplot(aes(x=year5, y=rse)) +
  geom_boxplot()+
  scale_y_continuous(limits=c(0,NA),
                     labels=percent_format(accuracy=1)) +
  facet_wrap(~gender) +
  labs(x="Year (5-year periods)",
       #title="Range (maximum-minimum)",
       y="Relative Standard Error (SE/LE)",
       color="", fill="", linetype="")+
  guides(linetype="none")+
  scale_x_discrete(labels=label_wrap(10))+
  isabel_theme+
  theme(legend.position = c(0.2, 0.2),
        legend.background = element_blank())
ggsave("../Tables & Figures/Appendix_Figure1.pdf", width=15, height=7.5)

# BIVARIATE MAPPING ---- 
# group based modeling trajectories 
# helps you identify what are the groups of trajectories that are here? 
#   aka commuting zones that started high then continued to go high, 
#   started low then went down, or started low then increased a lot! 
#   IDENTIFIES GROUPS OF TRAJECTORIES!

load('le_cbsa_results.rdata')

cz_hi <- c("134", "135") # CZ'S IN HAWAII
cz_ak <- seq(16,30) %>% as.character() # CZ'S IN ALASKA

results_cbsa_cz <- 
  results_cbsayr3 %>% ungroup() %>% mutate(type="cbsa", year_type="year3", year=year3, id=cbsa) %>% 
  select(-year3, -cbsa) %>% 
  bind_rows(
    results_cbsayr5 %>% ungroup() %>% mutate(type="cbsa", year_type="year5", year=year5, id=cbsa) %>% 
      select(-year5, -cbsa) 
  ) %>% bind_rows(
    results_czyr3 %>% ungroup() %>% mutate(type="cz", year_type="year3", year=year3, id=cz) %>% 
      select(-year3, -cz) 
  ) %>% bind_rows(
    results_czyr5 %>% ungroup() %>% mutate(type="cz", year_type="year5", year=year5, id=cz) %>% 
      select(-year5, -cz) 
  ) %>% 
  filter(!id%in%c("14720", "20780", "587")) %>% # QUARANTINING ROGUE CBSA/CZS FOR NOW 
  filter(!id%in%c(cz_ak, cz_hi)) # DROPPING CZ'S IN ALASKA AND HAWAII



## generating maps ----



# USING BIVARITE MAPPING
options<-expand.grid(gender=c("Men", "Women"), yr=c("year3", "year5"), unit=c("cz", "cbsa")) %>% as_tibble()
all_results_bivariate<-map(1:nrow(options), function(i){
  #i<-3
  opt_temp<-options %>% slice(i)
  # title<-paste0("LE among ", opt_temp$gender, ", ", opt_temp$unit, ", ", opt_temp$yr)
  
  if (opt_temp$unit=="cbsa") {
    temp <- right_join(shp_cbsa, results_cbsa_cz %>% get_bivarite(., gender_mw=opt_temp$gender, yr_35=opt_temp$yr, cbsa_cz=opt_temp$unit), by=c("GEOID"="id"))
    tmp_plot <- ggplot(temp, aes(geometry=geometry)) +
      geom_sf(mapping=aes(fill=bi_class), 
              color="white",
              size=0.1,
              show.legend=FALSE) + 
      bi_scale_fill(pal="DkViolet", dim=3) +
      geom_sf(data=st_transform(df_state, crs = st_crs(shp_cbsa)), size=0.1, color="black", fill=NA)+
      geom_sf(data=st_transform(df_mexico, crs=st_crs(shp_cbsa)), size=0.1, color="black", fill="darkgrey")+
      geom_sf(data=st_transform(df_canada, crs=st_crs(shp_cbsa)), size=0.1, color="black", fill="darkgrey")+
      geom_sf(data=st_transform(shp_census_region, crs = st_crs(shp_cbsa)), size=1.5, color="black", fill=NA)+
      geom_sf(data=st_transform(shp_census_division, crs = st_crs(shp_cbsa)), size=0.75, color="black", fill=NA)+
      coord_sf(xlim=st_bbox(temp)[c("xmin", "xmax")],
               ylim=st_bbox(temp)[c("ymin", "ymax")])+
      # labs(title=title)+
      map_theme
    tmp_legend <- bi_legend(pal = "DkViolet",
                            dim = 3,
                            xlab = "Higher LE",
                            ylab = "Most Increase Since 1990",
                            size = 8)+
      theme(panel.background = element_rect(fill = "lightblue1"),
            panel.border=element_blank(),
            plot.background = element_rect(fill = "lightblue1"),
            axis.ticks=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank())
    tmp_map <- ggdraw() +
      draw_plot(tmp_plot, 0, 0, 1, 1) +
      draw_plot(tmp_legend, 0.82, 0.2, 0.2, 0.2)
    tmp_map
    
  } else if (opt_temp$unit=="cz") {
    temp <- right_join(shp_cz, results_cbsa_cz %>% get_bivarite(., gender_mw=opt_temp$gender, yr_35=opt_temp$yr, cbsa_cz=opt_temp$unit), by=c("LM_Code"="id"))
    tmp_plot <- ggplot(temp, aes(geometry=geometry)) +
      geom_sf(mapping=aes(fill=bi_class), 
              color="white",
              size=0.1,
              show.legend=FALSE) + 
      bi_scale_fill(pal="DkViolet", dim=3) +
      geom_sf(data=st_transform(df_state, crs = st_crs(shp_cz)), size=0.1, color="black", fill=NA)+
      geom_sf(data=st_transform(shp_census_region, crs = st_crs(shp_cz)), size=1.5, color="black", fill=NA)+
      geom_sf(data=st_transform(shp_census_division, crs = st_crs(shp_cz)), size=0.75, color="black", fill=NA)+
      geom_sf(data=st_transform(df_mexico, crs=st_crs(shp_cz)), size=0.1, color="black", fill="darkgrey")+
      geom_sf(data=st_transform(df_canada, crs=st_crs(shp_cz)), size=0.1, color="black", fill="darkgrey")+
      geom_sf(data=temp %>% filter(LM_Code%in%"587") %>% select(geometry), size=0.1, color="black", fill="darkgrey")+
      coord_sf(xlim=st_bbox(temp)[c("xmin", "xmax")],
               ylim=st_bbox(temp)[c("ymin", "ymax")])+
      # labs(title=title)+
      map_theme
    tmp_legend <- bi_legend(pal = "DkViolet",
                            dim = 3,
                            xlab = "Higher LE",
                            ylab = "Most Increase Since 1990",
                            size = 8)+
      theme(panel.background = element_rect(fill = "lightblue1"),
            panel.border=element_blank(),
            plot.background = element_rect(fill = "lightblue1"),
            axis.ticks=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank())
    tmp_map <- ggdraw() +
      draw_plot(tmp_plot, 0, 0, 1, 1) +
      draw_plot(tmp_legend, 0.82, 0.2, 0.2, 0.2)
    tmp_map
  }
})


# USING BREWER PALETTE FOR MAPPING
options<-expand.grid(gender=c("Men", "Women"), yr=c("year3", "year5"), unit=c("cz", "cbsa")) %>% as_tibble()
all_results_brewer<-map(1:nrow(options), function(i){
  #i<-2
  opt_temp<-options %>% slice(i)
  # title<-paste0("LE among ", opt_temp$gender, ", ", opt_temp$unit, ", ", opt_temp$yr)
  
  if (opt_temp$unit=="cbsa") {
    temp <- right_join(shp_cbsa, results_cbsa_cz %>% get_bivarite(., gender_mw=opt_temp$gender, yr_35=opt_temp$yr, cbsa_cz=opt_temp$unit), by=c("GEOID"="id"))
    tmp_plot <- ggplot(temp, aes(geometry=geometry)) +
      geom_sf(mapping=aes(fill=bi_class), 
              color="white",
              size=0.1,
              show.legend=TRUE) + 
      geom_sf(data=st_transform(df_state, crs = st_crs(shp_cbsa)), size=0.1, color="black", fill=NA)+
      geom_sf(data=st_transform(df_mexico, crs=st_crs(shp_cbsa)), size=0.1, color="black", fill="darkgrey")+
      geom_sf(data=st_transform(df_canada, crs=st_crs(shp_cbsa)), size=0.1, color="black", fill="darkgrey")+
      coord_sf(xlim=st_bbox(temp)[c("xmin", "xmax")],
               ylim=st_bbox(temp)[c("ymin", "ymax")])+
      scale_fill_brewer(palette="BrBG")+
      # labs(title=title)+
      guides(color=guide_legend(reverse=TRUE))+
      map_theme + 
      theme(panel.background = element_rect(fill = "blue"),
            plot.background = element_rect(fill = "blue"))
    tmp_plot
    
  } else if (opt_temp$unit=="cz") {
    temp <- right_join(shp_cz, results_cbsa_cz %>% get_bivarite(., gender_mw=opt_temp$gender, yr_35=opt_temp$yr, cbsa_cz=opt_temp$unit), by=c("LM_Code"="id"))
    tmp_plot <- ggplot(temp, aes(geometry=geometry)) +
      geom_sf(mapping=aes(fill=bi_class), 
              color="white",
              size=0.1,
              show.legend=TRUE) + 
      geom_sf(data=st_transform(df_state, crs = st_crs(shp_cz)), size=0.1, color="black", fill=NA)+
      geom_sf(data=st_transform(df_mexico, crs=st_crs(shp_cz)), size=0.1, color="black", fill="darkgrey")+
      geom_sf(data=st_transform(df_canada, crs=st_crs(shp_cz)), size=0.1, color="black", fill="darkgrey")+
      coord_sf(xlim=st_bbox(temp)[c("xmin", "xmax")],
               ylim=st_bbox(temp)[c("ymin", "ymax")])+
      scale_fill_brewer(palette="BrBG")+
      # labs(title=title)+
      guides(color=guide_legend(reverse=TRUE))+
      map_theme
    tmp_plot
  }
})




## saving desired maps as PDF images ----

men_cz_5yr_map <- all_results_bivariate[[3]]
ggsave("../Tables & Figures/men_cz_5yr_map.pdf", men_cz_5yr_map, width=15, height=10)

women_cz_5yr_map <- all_results_bivariate[[4]]
ggsave("../Tables & Figures/women_cz_5yr_map.pdf", women_cz_5yr_map, width=15, height=10)



# SUMMARY STATISTICS OF LE ACROSS CZ OVER TIME ----
# What were the summary statistics for life expectancy over time within commuting zones: 
# mean, standard deviation, minimum, maximum? 

## overall CZ over time----
results_czyr5 %>% ungroup() %>% 
  select(year5, le) %>% 
  group_by(year5) %>% 
  summarise(mean_le=mean(le),
            sd_le=sd(le),
            min=min(le),
            max=max(le),
            median=median(le),
            q1=quantile(le, probs=0.25),
            q3=quantile(le, probs=0.75))
## overall CZ by gender over time-----
sumstat_yr_gender <- results_czyr5 %>% ungroup() %>% 
  select(year5, gender, le) %>% 
  group_by(year5, gender) %>% 
  summarise(median_le=median(le),
            sd_le=sd(le),
            min=min(le),
            max=max(le),
            q1=quantile(le, probs=0.25),
            q3=quantile(le, probs=0.75))


sumstat_table <- sumstat_yr_gender %>%  mutate(le_ci=paste0(median_le=format(median_le, digits=1, nsmall=1), 
                                           " ± (",
                                           format(sd_le, digits=1, nsmall=1),
                                           ")")) %>% # LE (95% CI) column 
  select(year5, gender, le_ci) %>% 
  spread(year5, le_ci) %>% full_join(
sumstat_yr_gender %>%  mutate(min_max=paste0("(",
                                             format(min, digits=1, nsmall=1), 
                                            ", ",
                                             format(max, digits=1, nsmall=1),
                                           ")")) %>%
  select(year5, gender, min_max) %>% 
  spread(year5, min_max) %>% 
  rename(col1=`1990-1994`,
         col2=`1995-1999`,
         col3=`2000-2004`,
         col4=`2005-2009`,
         col5=`2010-2014`,
         col6=`2015-2019`)) %>% 
  select(gender,
         `1990-1994`, col1,
         `1995-1999`, col2,
         `2000-2004`, col3,
         `2005-2009`, col4,
         `2010-2014`, col5,
         `2015-2019`, col6)
write.csv(sumstat_table, "../Tables & Figures/sumstat_table.csv", row.names=FALSE)

sumstat_table2 <- sumstat_yr_gender %>%  
  mutate(le_ci=paste0(median_le=format(median_le, digits=1, nsmall=1), 
                                                            " [",
                                                            format(q1, digits=1, nsmall=1), "-",
                      format(q3, digits=1, nsmall=1),
                                                            "]")) %>% 
  select(year5, gender, le_ci) %>% 
  pivot_wider(id_cols=year5, values_from=le_ci, names_from=gender)
write.csv(sumstat_table2, "../Tables & Figures/sumstat_table_v2.csv", row.names=FALSE)

# adding change statistics
sumstat_change<-results_czyr5 %>% ungroup() %>% 
  select(cz, year5, gender, le) %>% 
  pivot_wider(id_cols=c(gender, cz), values_from=le, names_from=year5) %>% 
  mutate(change=`2015-2019`-`1990-1994`) %>% 
  group_by(gender) %>% 
  summarise(median=median(change),
            q1=quantile(change, probs=0.25),
            q3=quantile(change, probs=0.75)) %>% 
  mutate(le_ci=paste0(median=format(median, digits=1, nsmall=1), 
                      " [",
                      format(q1, digits=1, nsmall=1), "-",
                      format(q3, digits=1, nsmall=1),
                      "]")) %>% 
  select(gender, le_ci) %>% 
  mutate(year5="Change") %>% 
  pivot_wider(id_cols=year5, names_from=gender, values_from=le_ci)
sumstat_table3<-sumstat_table2 %>% bind_rows(sumstat_change)
write.csv(sumstat_table3, "../Tables & Figures/sumstat_table_v3.csv", row.names=FALSE)

# VARIBILITY IN LE ACROSS CBSA'S OVER TIME ----

# NOTE: MAKE SURE TO RUN STEPS (1) - (5) IN PREP WORK PRIOR TO THIS 

# FINDING AVERAGE LE + ST_DEV + RANGE ACROSS CZ BY SEX, 1990-2019
results_coeffvar <- results_czyr5 %>% select(year5, cz, gender, le) %>% 
  group_by(year5, gender) %>%
  summarise(mean_le=mean(le),
            sd_le=sd(le),
            min=min(le),
            max=max(le)) %>% 
  mutate(coeff_var=sd_le/mean_le,
         range=max-min) %>% 
  mutate(year5=substr(year5, 1, 4)) %>% 
  mutate(year5=factor(year5,
                      levels=c(1990, 1995, 2000, 2005, 2010, 2015),
                      labels=c("1990- 1994",
                               "1995- 1999",
                               "2000- 2004",
                               "2005- 2009",
                               "2010- 2014",
                               "2015- 2019"))) 

coeffvar_men <- results_coeffvar %>% filter(gender%in%"Men") %>% 
  select(year5, coeff_var) %>% 
  spread(year5, coeff_var)
write.csv(coeffvar_men, "../Tables & Figures/coeffvar_men.csv", row.names=FALSE)

coeffvar_women <- results_coeffvar %>% filter(gender%in%"Women") %>% 
  select(year5, coeff_var) %>% 
  spread(year5, coeff_var)
write.csv(coeffvar_women, "../Tables & Figures/coeffvar_women.csv", row.names=FALSE)

coeffvar_trends <- 
  ggplot(results_coeffvar, aes(x=year5, y=coeff_var, group=gender))+
  geom_line(aes(linetype=gender))+
  scale_y_continuous(labels=percent_format(accuracy=1),
                     limits=c(0,NA)) +
  #facet_wrap(~gender) +
  geom_text(aes(label=format(coeff_var, digits=2, nsmall=2)), vjust=-1.5) +
  labs(x="5-year period",
       y="Coefficient of variation (%)",
       title="Coefficient of Variation (SD/mean)",
       color="", fill="", linetype="")+
  scale_x_discrete(labels=label_wrap(10))+
  isabel_theme+
  theme(legend.position = c(0.2, 0.2),
        legend.background = element_blank())

range_trends <- 
  ggplot(results_coeffvar, aes(x=year5, y=range, group=gender))+
  geom_line(aes(linetype=gender))+
  scale_y_continuous(limits=c(0,NA)) +
  #facet_wrap(~gender) +
  geom_text(aes(label=format(range, digits=1, nsmall=1)), vjust=-1.5) +
  labs(x="5-year period",
       title="Range (maximum-minimum)",
       y="Range (years)",
       color="", fill="", linetype="")+
  guides(linetype="none")+
  scale_x_discrete(labels=label_wrap(10))+
  isabel_theme+
  theme(legend.position = c(0.2, 0.2),
        legend.background = element_blank())
variability_trends<-arrangeGrob(grobs=list(coeffvar_trends,
                                           range_trends), ncol=2)
ggsave("../Tables & Figures/Figure2.pdf", variability_trends, width=15, height=7.5)



# FINDING OVERALL LE TREND FOR MEN AND WOMEN, CZ, 5YR PERIODS 

df <- dta %>% group_by(year5, age_5yr_group, gender) %>%  
  summarise(count=sum(count, na.rm=T),
            pop_denom=sum(pop_denom, na.rm=T)) %>%
  ungroup() %>%
  filter(gender=="Men") %>%
  # creating age-specific death rates variable, Mx
  mutate(mx=count/pop_denom,
         age_5yr_group=as.numeric(age_5yr_group)) %>%
  arrange(year5, age_5yr_group, gender) %>%
  group_by(year5, gender) %>%
  group_modify(~le_lv(., age_num=0, sex="M")) %>%
  bind_rows(
    dta %>% group_by(year5, age_5yr_group, gender) %>%
      summarise(count=sum(count, na.rm=T),
                pop_denom=sum(pop_denom, na.rm=T)) %>%
      ungroup() %>%
      filter(gender=="Women") %>%
      # creating age-specific death rates variable, Mx
      mutate(mx=count/pop_denom,
             age_5yr_group=as.numeric(age_5yr_group)) %>%
      arrange(year5, age_5yr_group, gender) %>%
      group_by(year5, gender) %>%
      group_modify(~le_lv(., age_num=0, sex="W"))
  ) %>% 
  ungroup() %>% 
  mutate(type="cz", year_type="year5", year=year5, id="overall") %>% 
  select(-year5) %>% 
  bind_rows(., results_cbsa_cz) %>% 
  mutate(year=substr(year, 1, 4)) %>% 
  mutate(year=factor(year,
                     levels=c(1990, 1995, 2000, 2005, 2010, 2015),
                     labels=c("1990- 1994",
                              "1995- 1999",
                              "2000- 2004",
                              "2005- 2009",
                              "2010- 2014",
                              "2015- 2019"))) 

cz_le_trends_men <- df %>% filter(gender%in%"Men", id%in%"overall") %>% 
  mutate(le_ci=paste0(le=format(le, digits=1, nsmall=1), 
                      " (",
                      format(lci, digits=1, nsmall=1),
                      ", ",
                      format(uci, digits=1, nsmall=1),
                      ")")) %>% # LE (95% CI) column
  select(year, le_ci) %>% 
  spread(year, le_ci)
write.csv(cz_le_trends_men, "../Tables & Figures/cz_le_trends_men.csv", row.names=FALSE)

cz_le_trends_women <- df %>% filter(gender%in%"Women", id%in%"overall") %>% 
  mutate(le_ci=paste0(le=format(le, digits=1, nsmall=1), 
                      " (",
                      format(lci, digits=1, nsmall=1),
                      ", ",
                      format(uci, digits=1, nsmall=1),
                      ")")) %>% # LE (95% CI) column
  select(year, le_ci) %>% 
  spread(year, le_ci)
write.csv(cz_le_trends_women, "../Tables & Figures/cz_le_trends_women.csv", row.names=FALSE)


czyr5_trends <- 
  ggplot(df %>% filter(type%in%"cz", year_type%in%"year5") %>% mutate(cols=ifelse(id%in%"overall", "2", "1")), 
         aes(x=year, y=le, color=cols, group=id))+
  geom_line(aes(color=cols, alpha=cols, size=cols))+
  facet_wrap(~gender) +
  scale_alpha_manual(values=c(0.25, 1))+
  scale_size_manual(values = c(0.5, 1))+
  scale_color_manual(labels=c("CZ", "Overall"), values=c("grey55", "black"))+
  labs(x="5-year period",
       y="Life expectancy (years)",
       color="", fill="")+
  scale_x_discrete(labels=label_wrap(10))+
  isabel_theme+
  guides(size=FALSE, alpha=FALSE) +
  theme(legend.position = c(0.8, 0.2),
        legend.background = element_blank())

ggsave("../Tables & Figures/figure1.pdf", czyr5_trends, width=15, height=10)


# GLOBAL MORAN'S I ----
# gives us a p-value and tells us if there's clustering / spatial autocorrelation
load('le_cbsa_results.rdata')
shp_cz <- read_sf("Shapefiles/CZ shapefile/ERS10.shp") %>% mutate(LM_Code=as.character(LM_Code))

cz_hi <- c("134", "135") # CZ'S IN HAWAII
cz_ak <- seq(16,30) %>% as.character() # CZ'S IN ALASKA

results_cbsa_cz <- 
  results_cbsayr3 %>% ungroup() %>% mutate(type="cbsa", year_type="year3", year=year3, id=cbsa) %>% 
  select(-year3, -cbsa) %>% 
  bind_rows(
    results_cbsayr5 %>% ungroup() %>% mutate(type="cbsa", year_type="year5", year=year5, id=cbsa) %>% 
      select(-year5, -cbsa) 
  ) %>% bind_rows(
    results_czyr3 %>% ungroup() %>% mutate(type="cz", year_type="year3", year=year3, id=cz) %>% 
      select(-year3, -cz) 
  ) %>% bind_rows(
    results_czyr5 %>% ungroup() %>% mutate(type="cz", year_type="year5", year=year5, id=cz) %>% 
      select(-year5, -cz) 
  ) %>% 
  filter(!id%in%c(cz_ak, cz_hi)) # DROPPING CZ'S IN ALASKA AND HAWAII


## Looping to do all periods/sexes
# add change first
results_cbsa_cz_moran<-results_cbsa_cz %>% bind_rows(
  bind_rows(results_cbsa_cz %>% get_bivarite(., gender_mw="Women", yr_35="year5", cbsa_cz="cz") %>% 
            mutate(year="Change") %>% 
            select(id, gender, year, le=abs_dif), 
          results_cbsa_cz %>% get_bivarite(., gender_mw="Men", yr_35="year5", cbsa_cz="cz") %>% 
            mutate(year="Change") %>% 
            select(id,gender, year, le=abs_dif)) %>% 
    mutate(type="cz", year_type="year5")) 

shp_clusters <- right_join(shp_cz, results_cbsa_cz_moran %>% 
                             filter(type%in%"cz",
                                    year_type%in%"year5"), 
                           by=c("LM_Code"="id")) %>% arrange(gender, LM_Code)
globalmoran_table<-shp_clusters %>% 
  group_by(gender, year) %>% 
  group_modify(~{
    #.x<-shp_clusters %>% filter(gender=="Men", year=="1990-1994")
    print(.y)
    neighbors <- poly2nb(.x %>% group_by(LM_Code))
    # fixing San Juan County, WA, CZ 604 (island, but manually input King County CZ 602 as neighbor)
    neighbors[[which(.x$LM_Code==604)]] <- which(.x$LM_Code%in%c(602))
    neighbors[[which(.x$LM_Code==602)]] <- c(neighbors[[which(.x$LM_Code==602)]], which(.x$LM_Code%in%c(604)))
    
    # it goes by indice in .x dataframe! for example: 
    # neighbors[[601]] is .x[601,] == CZ 92 (Philadelphia County, PA) 
    # neighbors of CZ 92 (Philly County) are 
    # 182 == .x[182,] == CZ 275 (Baltimore County, MD)
    # 307 == .x[307,] == CZ 390 (Camden County, NJ)
    # 308 == .x[308,] == CZ 391 (Bergen County, NJ)
    # 415 == .x[415,] == CZ 488 (Lancaster County, PA)
    # 418 == .x[418,] == CZ 490 (Berks County, PA)
    # 600 == .x[600,] == CZ 91  (Sussex County, DE)
    
    # use the moran.mc function
    temp <- moran.mc(x = .x %>% pull(le), 
                     # listw = list of neighbors, obtained from the nb adjacency matrix by using the nb2listw spdep function
                     listw = nb2listw(neighbors, style = "B"), 
                     # last, number of permutations for the permutation-based moran test
                     nsim = 9999)
    # extract moran statistic and p value
    globalmoran <- data.frame(statistic=temp$statistic[1], pval=temp$p.value) %>% 
      mutate(out=paste0(round(statistic, digits=3), " (", ifelse(pval<0.001, "<0.001", round(pval, digits=3)), ")")) 
    globalmoran
  })
# now create table
globalmoran_table_formatted<-globalmoran_table %>% select(year, gender, out) %>% 
  pivot_wider(id_cols=year, values_from=out, names_from=gender) %>% 
  rename(year5=year, men_moran=Men, women_moran=Women)
table1_new<-full_join(sumstat_table3, globalmoran_table_formatted) %>% 
  select(year5, Men, men_moran, Women, women_moran)
write.csv(table1_new, "../Tables & Figures/Table1_new.csv", row.names=FALSE)

# Map with Baseline and change in LEs
## Men
### Baseline

shp_rawle <- right_join(shp_cz, results_cbsa_cz_moran %>%
                             filter(type%in%"cz",
                                    year_type%in%"year5",
                                    year%in%c("1990-1994","Change"),
                                    gender%in%c("Men", "Women")),
                           by=c("LM_Code"="id")) %>% arrange(gender, LM_Code) %>% 
  filter(!LM_Code%in%"587") %>% 
  group_by(year, gender) %>% 
  mutate(jenks=cut(le, breaks=classIntervals(le,n=10,style="jenks")$brks, include.lowest = T))
  
le_colors<-rev(colorRampPalette(c("green","yellow", "orange","red"))(10))
maps_rawle<-shp_rawle %>% group_by(year, gender) %>% 
  group_map(~{
    #.x<-shp_rawle %>% filter(gender=="Men", year=="1990-1994")
    tag<-case_when(
      .y$gender=="Men" & .y$year=="1990-1994" ~ "A",
      .y$gender=="Men" & .y$year=="Change" ~ "C",
      .y$gender=="Women" & .y$year=="1990-1994" ~ "B",
      .y$gender=="Women" & .y$year=="Change" ~ "D",
    )
    ggplot()+
      geom_sf(data=.x, size=0,color=NA,
              aes(geometry=geometry, fill=(jenks)))+
      geom_sf(data=.x, size=.1,fill=NA,color="black",
              aes(geometry=geometry))+
      geom_sf(data=st_transform(df_state, crs = st_crs(shp_cz)), size=0.1, color="black", fill=NA)+
      geom_sf(data=st_transform(df_mexico, crs=st_crs(shp_cz)), size=0.1, color="black", fill="darkgrey")+
      geom_sf(data=st_transform(df_canada, crs=st_crs(shp_cz)), size=0.1, color="black", fill="darkgrey")+
      geom_sf(data=shp_cz %>% filter(LM_Code%in%"587") %>% select(geometry), size=0.1, color="black", fill="black")+
      geom_sf(data=st_transform(shp_census_region, crs = st_crs(shp_cz)), size=1.5, color="black", fill=NA)+
      geom_sf(data=st_transform(shp_census_division, crs = st_crs(shp_cz)), size=0.75, color="black", fill=NA)+
      annotate("text", x=st_bbox(shp_rawle)$xmin, 
               y=st_bbox(shp_rawle)$ymax, 
               hjust=0, vjust=0,
               size=10,
               label=tag)+
      scale_fill_manual(values=le_colors, name="")+
      coord_sf(xlim=st_bbox(shp_rawle)[c("xmin", "xmax")],
               ylim=st_bbox(shp_rawle)[c("ymin", "ymax")])+
      guides(size=F, alpha=F,
             fill=guide_legend(nrow=2)) +
      #labs(title="", tag="A")+
      map_theme +
      theme(legend.background = element_rect(fill='white'),
            legend.box.background = element_blank(),
            legend.text = element_text(size=19),
            legend.position="bottom")
    
  })
pall<-arrangeGrob(grobs=maps_rawle, 
                  ncol=2)
ggsave("../Tables & Figures/AppendixFigure2_rawle.pdf", pall, width=27.75, height=19)



# GETIS ORD ----

# USING GETIS ORD METHOD

## MEN ----

### baseline LE ----
shp_clusters <- right_join(shp_cz, results_cbsa_cz %>%
                             filter(type%in%"cz",
                                    year_type%in%"year5",
                                    year%in%"1990-1994",
                                    gender%in%"Men"),
                           by=c("LM_Code"="id")) %>% arrange(gender, LM_Code)

shp_clusters <- shp_clusters %>% filter(!LM_Code%in%"587")

queen_w <- queen_weights(shp_clusters)
gstar <- local_gstar(queen_w,shp_clusters %>% select(le))
shp_clusters$gstar_cluster <- lisa_clusters(gstar)
shp_clusters <- shp_clusters %>%
  mutate(gstar_cluster=factor(gstar_cluster, levels=0:4, labels=lisa_labels(gstar)))
table(shp_clusters$gstar_cluster)

# CREATING MAP FOR GETIS ORD: BASELINE LE FOR MEN
getisord_baselineLE_men <- ggplot()+
  geom_sf(data=shp_clusters %>% filter(!is.na(gstar_cluster)), size=0,color=NA,
          aes(geometry=geometry, fill=(gstar_cluster)))+
  geom_sf(data=shp_clusters, size=.1,fill=NA,color="black",
          aes(geometry=geometry))+
  geom_sf(data=st_transform(df_state, crs = st_crs(shp_cz)), size=0.1, color="black", fill=NA)+
  geom_sf(data=st_transform(df_mexico, crs=st_crs(shp_cz)), size=0.1, color="black", fill="darkgrey")+
  geom_sf(data=st_transform(df_canada, crs=st_crs(shp_cz)), size=0.1, color="black", fill="darkgrey")+
  geom_sf(data=shp_cz %>% filter(LM_Code%in%"587") %>% select(geometry), size=0.1, color="black", fill="black")+
  geom_sf(data=st_transform(shp_census_region, crs = st_crs(shp_cz)), size=1.5, color="black", fill=NA)+
  geom_sf(data=st_transform(shp_census_division, crs = st_crs(shp_cz)), size=0.75, color="black", fill=NA)+
  annotate("text", x=st_bbox(shp_clusters)$xmin, 
           y=st_bbox(shp_clusters)$ymax, 
           hjust=0, vjust=0,
           size=10,
           label="A")+
  scale_fill_manual(values=moran_colors, name="",
                    labels=c("Not Significant", "High or Increase", "Low or Decrease", "N/A"))+
  coord_sf(xlim=st_bbox(shp_clusters)[c("xmin", "xmax")],
           ylim=st_bbox(shp_clusters)[c("ymin", "ymax")])+
  guides(size=F, alpha=F,
         fill=guide_legend(nrow=1)) +
  #labs(title="", tag="A")+
  map_theme +
  theme(legend.background = element_blank(),
        legend.box.background = element_blank())
legend_sep<-get_legend(getisord_baselineLE_men);plot(legend_sep)
getisord_baselineLE_men<-getisord_baselineLE_men+guides(fill="none")
#ggsave("test.pdf", getisord_baselineLE_men, width=15, height=7.5)

### change in LE ----
shp_clusters <- right_join(shp_cz, results_cbsa_cz %>% get_bivarite(., gender_mw="Men", yr_35="year5", cbsa_cz="cz"), by=c("LM_Code"="id")) 
shp_clusters <- shp_clusters %>% filter(!LM_Code%in%"587")

queen_w <- queen_weights(shp_clusters)
gstar <- local_gstar(queen_w,shp_clusters %>% select(abs_dif))
shp_clusters$gstar_cluster <- lisa_clusters(gstar)
shp_clusters <- shp_clusters %>%
  mutate(gstar_cluster=factor(gstar_cluster, levels=0:4, labels=lisa_labels(gstar)))
table(shp_clusters$gstar_cluster)


# CREATING MAP FOR GETIS ORD: DIFFERENCE IN LIFE EXPECTANCY B/W 2015-2019 AND 1990-1994 FOR MEN 
getisord_diffLE_men <- ggplot()+
  geom_sf(data=shp_clusters %>% filter(!is.na(gstar_cluster)), size=0,color=NA,
          aes(geometry=geometry, fill=(gstar_cluster)))+
  geom_sf(data=shp_clusters, size=.1,fill=NA,color="black",
          aes(geometry=geometry))+
  geom_sf(data=st_transform(df_state, crs = st_crs(shp_cz)), size=0.1, color="black", fill=NA)+
  geom_sf(data=st_transform(df_mexico, crs=st_crs(shp_cz)), size=0.1, color="black", fill="darkgrey")+
  geom_sf(data=st_transform(df_canada, crs=st_crs(shp_cz)), size=0.1, color="black", fill="darkgrey")+
  geom_sf(data=shp_cz %>% filter(LM_Code%in%"587") %>% select(geometry), size=0.1, color="black", fill="black")+
  geom_sf(data=st_transform(shp_census_region, crs = st_crs(shp_cz)), size=1.5, color="black", fill=NA)+
  geom_sf(data=st_transform(shp_census_division, crs = st_crs(shp_cz)), size=0.75, color="black", fill=NA)+
  annotate("text", x=st_bbox(shp_clusters)$xmin, 
           y=st_bbox(shp_clusters)$ymax, 
           hjust=0, vjust=0,
           size=10,
           label="C")+
  scale_fill_manual(values=moran_colors, name="")+
  coord_sf(xlim=st_bbox(shp_clusters)[c("xmin", "xmax")],
           ylim=st_bbox(shp_clusters)[c("ymin", "ymax")])+
  guides(size=F, alpha=F, fill = guide_legend(override.aes = list(alpha=0))) +
  #labs(title=title) +
  #labs(title="", tag="C")+
  map_theme + 
  theme(legend.text=element_blank()) +
  theme(legend.position="none")


### biscale 3x3 plot: baseline LE & change in LE ----

# baseline df
shp_clusters <- right_join(shp_cz, results_cbsa_cz %>%
                             filter(type%in%"cz",
                                    year_type%in%"year5",
                                    year%in%"1990-1994",
                                    gender%in%"Men"),
                           by=c("LM_Code"="id")) %>% arrange(gender, LM_Code)
shp_clusters <- shp_clusters %>% filter(!LM_Code%in%"587")

queen_w <- queen_weights(shp_clusters)
gstar <- local_gstar(queen_w,shp_clusters %>% select(le))
shp_clusters$gstar_cluster <- lisa_clusters(gstar)
shp_clusters <- shp_clusters %>%
  mutate(gstar_cluster=factor(gstar_cluster, levels=0:4, labels=lisa_labels(gstar)))
table(shp_clusters$gstar_cluster)
baseline_clusters <- shp_clusters %>% 
  mutate(baseline_gstar_cluster=gstar_cluster) %>% 
  select(-gstar_cluster) %>% 
  st_drop_geometry()

# change in LE df 
shp_clusters <- right_join(shp_cz, results_cbsa_cz %>% get_bivarite(., gender_mw="Men", yr_35="year5", cbsa_cz="cz"), by=c("LM_Code"="id")) 
shp_clusters <- shp_clusters %>% filter(!LM_Code%in%"587")

queen_w <- queen_weights(shp_clusters)
gstar <- local_gstar(queen_w,shp_clusters %>% select(abs_dif))
shp_clusters$gstar_cluster <- lisa_clusters(gstar)
shp_clusters <- shp_clusters %>%
  mutate(gstar_cluster=factor(gstar_cluster, levels=0:4, labels=lisa_labels(gstar)))
table(shp_clusters$gstar_cluster)
diffLE_clusters <- shp_clusters %>% 
  mutate(diffLE_gstar_cluster=gstar_cluster) %>%
  select(-gstar_cluster) %>% 
  st_drop_geometry()

# combining both df's and droplevel unused factor levels (Undefined, Isolated)
df <- baseline_clusters %>% right_join(diffLE_clusters, by="LM_Code") %>% 
  select(LM_Code, baseline_gstar_cluster, diffLE_gstar_cluster) %>% 
  filter(!baseline_gstar_cluster%in%c("Undefined", "Isolated")) %>% 
  filter(!diffLE_gstar_cluster%in%c("Undefined", "Isolated"))
df$baseline_gstar_cluster <- droplevels(df$baseline_gstar_cluster)
df$diffLE_gstar_cluster <- droplevels(df$diffLE_gstar_cluster)
df <- df %>% mutate(baseline_gstar_cluster=case_when(baseline_gstar_cluster%in%"Low-Low" ~ "Low",
                                                     baseline_gstar_cluster%in%"High-High"~ "High",
                                                     baseline_gstar_cluster%in%"Not significant" ~ "Not significant"),
                    diffLE_gstar_cluster=case_when(diffLE_gstar_cluster%in%"Low-Low" ~ "Low",
                                                   diffLE_gstar_cluster%in%"High-High"~ "High",
                                                   diffLE_gstar_cluster%in%"Not significant" ~ "Not significant"))


# checking unique combos, should be 9
crossing(df$baseline_gstar_cluster, df$diffLE_gstar_cluster)
df %>% count(baseline_gstar_cluster, diffLE_gstar_cluster) 


df <- df %>% mutate(type=case_when(
  baseline_gstar_cluster=="Low" & diffLE_gstar_cluster=="Low" ~ "Low Baseline - Decreased",
  baseline_gstar_cluster=="Low" & diffLE_gstar_cluster=="Not significant" ~ "Low Baseline - NS Change",
  baseline_gstar_cluster=="Low" & diffLE_gstar_cluster=="High" ~ "Low Baseline - Increased",
  baseline_gstar_cluster=="Not significant" & diffLE_gstar_cluster=="Low" ~ "NS Baseline - Decreased",
  baseline_gstar_cluster=="Not significant" & diffLE_gstar_cluster=="Not significant" ~ "NS Baseline - NS Change",
  baseline_gstar_cluster=="Not significant" & diffLE_gstar_cluster=="High" ~ "NS Baseline - Increased",
  baseline_gstar_cluster=="High" & diffLE_gstar_cluster=="Low" ~ "High Baseline - Decreased",
  baseline_gstar_cluster=="High" & diffLE_gstar_cluster=="Not significant" ~ "High Baseline - NS Change",
  baseline_gstar_cluster=="High" & diffLE_gstar_cluster=="High" ~ "High Baseline - Increased"))

  
df <- df %>% mutate(type=factor(type,
                           levels=c("Low Baseline - Decreased",
                                    "Low Baseline - NS Change",
                                    "Low Baseline - Increased",
                                    "NS Baseline - Decreased",
                                    "NS Baseline - NS Change",
                                    "NS Baseline - Increased", 
                                    "High Baseline - Decreased",
                                    "High Baseline - NS Change",
                                    "High Baseline - Increased"))) %>% 
  arrange(type)
cluster_count <- df %>% count(type)

biscale_men_df <- right_join(shp_cz, df, by="LM_Code") %>% glimpse()

# CREATING MAP FOR BASELINE LE X CHANGE IN LE 
biscale_men_map <- ggplot()+
  geom_sf(data=biscale_men_df, size=0,color=NA,
          aes(geometry=geometry, fill=(type)))+
  geom_sf(data=biscale_men_df, size=.1,fill=NA,color="black",
          aes(geometry=geometry))+
  geom_sf(data=st_transform(df_state, crs = st_crs(shp_cz)), size=0.1, color="black", fill=NA)+
  geom_sf(data=st_transform(df_mexico, crs=st_crs(shp_cz)), size=0.1, color="black", fill="darkgrey")+
  geom_sf(data=st_transform(df_canada, crs=st_crs(shp_cz)), size=0.1, color="black", fill="darkgrey")+
  geom_sf(data=shp_cz %>% filter(LM_Code%in%"587") %>% select(geometry), size=0.1, color="black", fill="black")+
  geom_sf(data=st_transform(shp_census_region, crs = st_crs(shp_cz)), size=1.5, color="black", fill=NA)+
  geom_sf(data=st_transform(shp_census_division, crs = st_crs(shp_cz)), size=0.75, color="black", fill=NA)+
  annotate("text", x=st_bbox(biscale_men_df)$xmin, 
           y=st_bbox(biscale_men_df)$ymax, 
           hjust=0, vjust=0,
           size=10,
           label="A")+
  scale_fill_manual(values=cols)+
  coord_sf(xlim=st_bbox(biscale_men_df)[c("xmin", "xmax")],
           ylim=st_bbox(biscale_men_df)[c("ymin", "ymax")])+
  guides(color=guide_legend(reverse=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE))+
  #labs(title="Men")+
  #labs(tag="A")+
  map_theme+
  theme(legend.position="bottom", legend.title=element_blank())+
  theme(legend.background = element_blank(),
        legend.box.background = element_blank())


## WOMEN ----

### baseline LE ----
shp_clusters <- right_join(shp_cz, results_cbsa_cz %>%
                             filter(type%in%"cz",
                                    year_type%in%"year5",
                                    year%in%"1990-1994",
                                    gender%in%"Women"),
                           by=c("LM_Code"="id")) %>% arrange(gender, LM_Code)
shp_clusters <- shp_clusters %>% filter(!LM_Code%in%"587")

queen_w <- queen_weights(shp_clusters)
gstar <- local_gstar(queen_w,shp_clusters %>% select(le))
shp_clusters$gstar_cluster <- lisa_clusters(gstar)
shp_clusters <- shp_clusters %>%
  mutate(gstar_cluster=factor(gstar_cluster, levels=0:4, labels=lisa_labels(gstar)))
table(shp_clusters$gstar_cluster)

# CREATING MAP FOR GETIS ORD: BASELINE LE FOR WOMEN 
getisord_baselineLE_women <- ggplot()+
  geom_sf(data=shp_clusters %>% filter(!is.na(gstar_cluster)), size=0,color=NA,
          aes(geometry=geometry, fill=(gstar_cluster)))+
  geom_sf(data=shp_clusters, size=.1,fill=NA,color="black",
          aes(geometry=geometry))+
  geom_sf(data=st_transform(df_state, crs = st_crs(shp_cz)), size=0.1, color="black", fill=NA)+
  geom_sf(data=st_transform(df_mexico, crs=st_crs(shp_cz)), size=0.1, color="black", fill="darkgrey")+
  geom_sf(data=st_transform(df_canada, crs=st_crs(shp_cz)), size=0.1, color="black", fill="darkgrey")+
  geom_sf(data=shp_cz %>% filter(LM_Code%in%"587") %>% select(geometry), size=0.1, color="black", fill="black")+
  geom_sf(data=st_transform(shp_census_region, crs = st_crs(shp_cz)), size=1.5, color="black", fill=NA)+
  geom_sf(data=st_transform(shp_census_division, crs = st_crs(shp_cz)), size=0.75, color="black", fill=NA)+
  annotate("text", x=st_bbox(shp_clusters)$xmin, 
           y=st_bbox(shp_clusters)$ymax, 
           hjust=0, vjust=0,
           size=10,
           label="B")+
  scale_fill_manual(values=moran_colors, name="")+
  coord_sf(xlim=st_bbox(shp_clusters)[c("xmin", "xmax")],
           ylim=st_bbox(shp_clusters)[c("ymin", "ymax")])+
  guides(size=F, alpha=F, fill = guide_legend(override.aes = list(alpha=0))) +
  #labs(title=title) +
  #labs(title="", tag="B")+
  map_theme + 
  theme(legend.text=element_blank()) +
  theme(legend.position="none")


### change in LE ----

shp_clusters <- right_join(shp_cz, results_cbsa_cz %>% get_bivarite(., gender_mw="Women", yr_35="year5", cbsa_cz="cz"), by=c("LM_Code"="id"))
shp_clusters <- shp_clusters %>% filter(!LM_Code%in%"587")

queen_w <- queen_weights(shp_clusters)
gstar <- local_gstar(queen_w,shp_clusters %>% select(abs_dif))
shp_clusters$gstar_cluster <- lisa_clusters(gstar)
shp_clusters <- shp_clusters %>%
  mutate(gstar_cluster=factor(gstar_cluster, levels=0:4, labels=lisa_labels(gstar)))
table(shp_clusters$gstar_cluster)

# CREATING MAP FOR GETIS ORD: DIFFERENCE IN LIFE EXPECTANCY B/W 2015-2019 AND 1990-1994 FOR WOMEN 
getisord_diffLE_women <- ggplot()+
  geom_sf(data=shp_clusters %>% filter(!is.na(gstar_cluster)), size=0,color=NA,
          aes(geometry=geometry, fill=(gstar_cluster)))+
  geom_sf(data=shp_clusters, size=.1,fill=NA,color="black",
          aes(geometry=geometry))+
  geom_sf(data=st_transform(df_state, crs = st_crs(shp_cz)), size=0.1, color="black", fill=NA)+
  geom_sf(data=st_transform(df_mexico, crs=st_crs(shp_cz)), size=0.1, color="black", fill="darkgrey")+
  geom_sf(data=st_transform(df_canada, crs=st_crs(shp_cz)), size=0.1, color="black", fill="darkgrey")+
  geom_sf(data=shp_cz %>% filter(LM_Code%in%"587") %>% select(geometry), size=0.1, color="black", fill="black")+
  geom_sf(data=st_transform(shp_census_region, crs = st_crs(shp_cz)), size=1.5, color="black", fill=NA)+
  geom_sf(data=st_transform(shp_census_division, crs = st_crs(shp_cz)), size=0.75, color="black", fill=NA)+
  annotate("text", x=st_bbox(shp_clusters)$xmin, 
           y=st_bbox(shp_clusters)$ymax, 
           hjust=0, vjust=0,
           size=10,
           label="D")+
  scale_fill_manual(values=moran_colors, name="")+
  coord_sf(xlim=st_bbox(shp_clusters)[c("xmin", "xmax")],
           ylim=st_bbox(shp_clusters)[c("ymin", "ymax")])+
  guides(size=F, alpha=F, fill = guide_legend(override.aes = list(alpha=0))) +
  #labs(title=title) +
  #labs(title="", tag="D")+
  map_theme + 
  theme(legend.text=element_blank()) +
  theme(legend.position="none")

### biscale 3x3 plot: baseline LE & change in LE ----

# baseline df
shp_clusters <- right_join(shp_cz, results_cbsa_cz %>%
                             filter(type%in%"cz",
                                    year_type%in%"year5",
                                    year%in%"1990-1994",
                                    gender%in%"Women"),
                           by=c("LM_Code"="id")) %>% arrange(gender, LM_Code)
shp_clusters <- shp_clusters %>% filter(!LM_Code%in%"587")

queen_w <- queen_weights(shp_clusters)
gstar <- local_gstar(queen_w,shp_clusters %>% select(le))
shp_clusters$gstar_cluster <- lisa_clusters(gstar)
shp_clusters <- shp_clusters %>%
  mutate(gstar_cluster=factor(gstar_cluster, levels=0:4, labels=lisa_labels(gstar)))
table(shp_clusters$gstar_cluster)
baseline_clusters <- shp_clusters %>% 
  mutate(baseline_gstar_cluster=gstar_cluster) %>% 
  select(-gstar_cluster) %>% 
  st_drop_geometry()

# change in LE df 
shp_clusters <- right_join(shp_cz, results_cbsa_cz %>% get_bivarite(., gender_mw="Women", yr_35="year5", cbsa_cz="cz"), by=c("LM_Code"="id")) 
shp_clusters <- shp_clusters %>% filter(!LM_Code%in%"587")

queen_w <- queen_weights(shp_clusters)
gstar <- local_gstar(queen_w,shp_clusters %>% select(abs_dif))
shp_clusters$gstar_cluster <- lisa_clusters(gstar)
shp_clusters <- shp_clusters %>%
  mutate(gstar_cluster=factor(gstar_cluster, levels=0:4, labels=lisa_labels(gstar)))
table(shp_clusters$gstar_cluster)
diffLE_clusters <- shp_clusters %>% 
  mutate(diffLE_gstar_cluster=gstar_cluster) %>%
  select(-gstar_cluster) %>% 
  st_drop_geometry()

# combining both df's and droplevel unused factor levels (Undefined, Isolated)
df <- baseline_clusters %>% right_join(diffLE_clusters, by="LM_Code") %>% 
  select(LM_Code, baseline_gstar_cluster, diffLE_gstar_cluster) %>% 
  filter(!baseline_gstar_cluster%in%c("Undefined", "Isolated")) %>% 
  filter(!diffLE_gstar_cluster%in%c("Undefined", "Isolated"))
df$baseline_gstar_cluster <- droplevels(df$baseline_gstar_cluster)
df$diffLE_gstar_cluster <- droplevels(df$diffLE_gstar_cluster)
df <- df %>% mutate(baseline_gstar_cluster=case_when(baseline_gstar_cluster%in%"Low-Low" ~ "Low",
                                                     baseline_gstar_cluster%in%"High-High"~ "High",
                                                     baseline_gstar_cluster%in%"Not significant" ~ "Not significant"),
                    diffLE_gstar_cluster=case_when(diffLE_gstar_cluster%in%"Low-Low" ~ "Low",
                                                   diffLE_gstar_cluster%in%"High-High"~ "High",
                                                   diffLE_gstar_cluster%in%"Not significant" ~ "Not significant"))


# checking unique combos, should be 9
crossing(df$baseline_gstar_cluster, df$diffLE_gstar_cluster)
df %>% count(baseline_gstar_cluster, diffLE_gstar_cluster) 


df <- df %>% mutate(type=case_when(
  baseline_gstar_cluster=="Low" & diffLE_gstar_cluster=="Low" ~ "Low Baseline - Decreased",
  baseline_gstar_cluster=="Low" & diffLE_gstar_cluster=="Not significant" ~ "Low Baseline - NS Change",
  baseline_gstar_cluster=="Low" & diffLE_gstar_cluster=="High" ~ "Low Baseline - Increased",
  baseline_gstar_cluster=="Not significant" & diffLE_gstar_cluster=="Low" ~ "NS Baseline - Decreased",
  baseline_gstar_cluster=="Not significant" & diffLE_gstar_cluster=="Not significant" ~ "NS Baseline - NS Change",
  baseline_gstar_cluster=="Not significant" & diffLE_gstar_cluster=="High" ~ "NS Baseline - Increased",
  baseline_gstar_cluster=="High" & diffLE_gstar_cluster=="Low" ~ "High Baseline - Decreased",
  baseline_gstar_cluster=="High" & diffLE_gstar_cluster=="Not significant" ~ "High Baseline - NS Change",
  baseline_gstar_cluster=="High" & diffLE_gstar_cluster=="High" ~ "High Baseline - Increased"))


df <- df %>% mutate(type=factor(type,
                                levels=c("Low Baseline - Decreased",
                                         "Low Baseline - NS Change",
                                         "Low Baseline - Increased",
                                         "NS Baseline - Decreased",
                                         "NS Baseline - NS Change",
                                         "NS Baseline - Increased", 
                                         "High Baseline - Decreased",
                                         "High Baseline - NS Change",
                                         "High Baseline - Increased"))) %>% 
  arrange(type)
cluster_count <- df %>% count(type)

biscale_women_df <- right_join(shp_cz, df, by="LM_Code") 

# CREATING MAP FOR BASELINE LE X CHANGE IN LE 
biscale_women_map <- ggplot()+
  geom_sf(data=biscale_women_df, size=0,color=NA,
          aes(geometry=geometry, fill=(type)))+
  geom_sf(data=biscale_women_df, size=.1,fill=NA,color="black",
          aes(geometry=geometry))+
  geom_sf(data=st_transform(df_state, crs = st_crs(shp_cz)), size=0.1, color="black", fill=NA)+
  geom_sf(data=st_transform(df_mexico, crs=st_crs(shp_cz)), size=0.1, color="black", fill="darkgrey")+
  geom_sf(data=st_transform(df_canada, crs=st_crs(shp_cz)), size=0.1, color="black", fill="darkgrey")+
  geom_sf(data=shp_cz %>% filter(LM_Code%in%"587") %>% select(geometry), size=0.1, color="black", fill="black")+
  geom_sf(data=st_transform(shp_census_region, crs = st_crs(shp_cz)), size=1.5, color="black", fill=NA)+
  geom_sf(data=st_transform(shp_census_division, crs = st_crs(shp_cz)), size=0.75, color="black", fill=NA)+
  annotate("text", x=st_bbox(biscale_men_df)$xmin, 
           y=st_bbox(biscale_men_df)$ymax, 
           hjust=0, vjust=0,
           size=10,
           label="B")+
  scale_fill_manual(values=cols)+
  coord_sf(xlim=st_bbox(biscale_women_df)[c("xmin", "xmax")],
           ylim=st_bbox(biscale_women_df)[c("ymin", "ymax")])+
  guides(color=guide_legend(reverse=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE))+
  #labs(subtitle="Women")+
  #labs(tag="B")+
  map_theme+
  theme(legend.position="bottom", legend.title=element_blank())+
  theme(legend.background = element_blank(),
        legend.box.background = element_blank())



pall<-arrangeGrob(grobs=list(getisord_baselineLE_men, 
                             getisord_baselineLE_women, 
                             getisord_diffLE_men, 
                             getisord_diffLE_women), 
                  ncol=2)
pall<-arrangeGrob(grobs=list(pall, legend_sep), heights=c(20, 1), ncol=1)
ggsave("../Tables & Figures/Figure3_sep.pdf", pall, width=29, height=19)


legend<-get_legend(biscale_men_map)
biscale_men_map<-biscale_men_map+guides(fill="none") 
biscale_women_map<-biscale_women_map+guides(fill="none")
map_all<-arrangeGrob(grobs=list(biscale_men_map, biscale_women_map), ncol=2)
map_all_forreal<-arrangeGrob(grobs=list(map_all, legend), ncol=1, heights=c(15, 2))
ggsave("../Tables & Figures/Figure3_combined.pdf", map_all_forreal, width=30, height=11)


# exploring size of CZs by type of clustering
# first, create average pop across periods
pop<-dta %>% group_by(cz, year) %>% 
  summarize(pop=sum(pop_denom)) %>% 
  group_by(cz) %>% 
  summarise(avg_pop=mean(pop)) %>% 
  rename(LM_Code=cz)
biscale_pop<-bind_rows(biscale_women_df %>% as_tibble() %>% left_join(pop) %>% 
            select(LM_Code, type, avg_pop) %>% mutate(gender="Women"),
          biscale_men_df %>% as_tibble() %>% left_join(pop) %>% 
            select(LM_Code, type, avg_pop) %>% mutate(gender="Men")) %>% 
  mutate(type=sub(" - ", " -\n", type))
ggplot(biscale_pop, aes(x=type, y=avg_pop)) +
  geom_boxplot()+
  scale_y_continuous(trans="log10", breaks=10^(4:7),
                     labels=c("10K", "100K", "1M", "10M")) +
  facet_wrap(~gender) +
  labs(x="Year (5-year periods)",
       #title="Range (maximum-minimum)",
       y="Relative Standard Error (SE/LE)",
       color="", fill="", linetype="")+
  guides(linetype="none")+
  scale_x_discrete(labels=label_wrap(40))+
  isabel_theme+
  theme(legend.position = c(0.2, 0.2),
        legend.background = element_blank(),
        axis.text.x=element_text(angle=90))
  
  


# USING GETIS ORD METHOD WITH P<0.01 (OR TRY RESAMPLING LE?? AND THEN SHOW CLUSTERS WITH >95% prob of being there?)

## MEN ----

### baseline LE ----
shp_clusters <- right_join(shp_cz, results_cbsa_cz %>%
                             filter(type%in%"cz",
                                    year_type%in%"year5",
                                    year%in%"1990-1994",
                                    gender%in%"Men"),
                           by=c("LM_Code"="id")) %>% arrange(gender, LM_Code)

shp_clusters <- shp_clusters %>% filter(!LM_Code%in%"587")

queen_w <- queen_weights(shp_clusters)
gstar <- local_gstar(queen_w,shp_clusters %>% select(le), significance_cutoff = 0.01)
shp_clusters$gstar_cluster <- lisa_clusters(gstar)
shp_clusters <- shp_clusters %>%
  mutate(gstar_cluster=factor(gstar_cluster, levels=0:4, labels=lisa_labels(gstar)))
table(shp_clusters$gstar_cluster)

# CREATING MAP FOR GETIS ORD: BASELINE LE FOR MEN
getisord_baselineLE_men <- ggplot()+
  geom_sf(data=shp_clusters %>% filter(!is.na(gstar_cluster)), size=0,color=NA,
          aes(geometry=geometry, fill=(gstar_cluster)))+
  geom_sf(data=shp_clusters, size=.1,fill=NA,color="black",
          aes(geometry=geometry))+
  geom_sf(data=st_transform(df_state, crs = st_crs(shp_cz)), size=0.1, color="black", fill=NA)+
  geom_sf(data=st_transform(df_mexico, crs=st_crs(shp_cz)), size=0.1, color="black", fill="darkgrey")+
  geom_sf(data=st_transform(df_canada, crs=st_crs(shp_cz)), size=0.1, color="black", fill="darkgrey")+
  geom_sf(data=shp_cz %>% filter(LM_Code%in%"587") %>% select(geometry), size=0.1, color="black", fill="black")+
  geom_sf(data=st_transform(shp_census_region, crs = st_crs(shp_cz)), size=1.5, color="black", fill=NA)+
  geom_sf(data=st_transform(shp_census_division, crs = st_crs(shp_cz)), size=0.75, color="black", fill=NA)+
  annotate("text", x=st_bbox(shp_clusters)$xmin, 
           y=st_bbox(shp_clusters)$ymax, 
           hjust=0, vjust=0,
           size=10,
           label="A")+
  scale_fill_manual(values=moran_colors, name="",
                    labels=c("Not Significant", "High or Increase", "Low or Decrease", "N/A"))+
  coord_sf(xlim=st_bbox(shp_clusters)[c("xmin", "xmax")],
           ylim=st_bbox(shp_clusters)[c("ymin", "ymax")])+
  guides(size=F, alpha=F,
         fill=guide_legend(nrow=1)) +
  #labs(title="", tag="A")+
  map_theme +
  theme(legend.background = element_blank(),
        legend.box.background = element_blank())
legend_sep<-get_legend(getisord_baselineLE_men);plot(legend_sep)
getisord_baselineLE_men<-getisord_baselineLE_men+guides(fill="none")
#ggsave("test.pdf", getisord_baselineLE_men, width=15, height=7.5)

### change in LE ----
shp_clusters <- right_join(shp_cz, results_cbsa_cz %>% get_bivarite(., gender_mw="Men", yr_35="year5", cbsa_cz="cz"), by=c("LM_Code"="id")) 
shp_clusters <- shp_clusters %>% filter(!LM_Code%in%"587")

queen_w <- queen_weights(shp_clusters)
gstar <- local_gstar(queen_w,shp_clusters %>% select(abs_dif))
shp_clusters$gstar_cluster <- lisa_clusters(gstar)
shp_clusters <- shp_clusters %>%
  mutate(gstar_cluster=factor(gstar_cluster, levels=0:4, labels=lisa_labels(gstar)))
table(shp_clusters$gstar_cluster)


# CREATING MAP FOR GETIS ORD: DIFFERENCE IN LIFE EXPECTANCY B/W 2015-2019 AND 1990-1994 FOR MEN 
getisord_diffLE_men <- ggplot()+
  geom_sf(data=shp_clusters %>% filter(!is.na(gstar_cluster)), size=0,color=NA,
          aes(geometry=geometry, fill=(gstar_cluster)))+
  geom_sf(data=shp_clusters, size=.1,fill=NA,color="black",
          aes(geometry=geometry))+
  geom_sf(data=st_transform(df_state, crs = st_crs(shp_cz)), size=0.1, color="black", fill=NA)+
  geom_sf(data=st_transform(df_mexico, crs=st_crs(shp_cz)), size=0.1, color="black", fill="darkgrey")+
  geom_sf(data=st_transform(df_canada, crs=st_crs(shp_cz)), size=0.1, color="black", fill="darkgrey")+
  geom_sf(data=shp_cz %>% filter(LM_Code%in%"587") %>% select(geometry), size=0.1, color="black", fill="black")+
  geom_sf(data=st_transform(shp_census_region, crs = st_crs(shp_cz)), size=1.5, color="black", fill=NA)+
  geom_sf(data=st_transform(shp_census_division, crs = st_crs(shp_cz)), size=0.75, color="black", fill=NA)+
  annotate("text", x=st_bbox(shp_clusters)$xmin, 
           y=st_bbox(shp_clusters)$ymax, 
           hjust=0, vjust=0,
           size=10,
           label="C")+
  scale_fill_manual(values=moran_colors, name="")+
  coord_sf(xlim=st_bbox(shp_clusters)[c("xmin", "xmax")],
           ylim=st_bbox(shp_clusters)[c("ymin", "ymax")])+
  guides(size=F, alpha=F, fill = guide_legend(override.aes = list(alpha=0))) +
  #labs(title=title) +
  #labs(title="", tag="C")+
  map_theme + 
  theme(legend.text=element_blank()) +
  theme(legend.position="none")


### biscale 3x3 plot: baseline LE & change in LE ----

# baseline df
shp_clusters <- right_join(shp_cz, results_cbsa_cz %>%
                             filter(type%in%"cz",
                                    year_type%in%"year5",
                                    year%in%"1990-1994",
                                    gender%in%"Men"),
                           by=c("LM_Code"="id")) %>% arrange(gender, LM_Code)
shp_clusters <- shp_clusters %>% filter(!LM_Code%in%"587")

queen_w <- queen_weights(shp_clusters)
gstar <- local_gstar(queen_w,shp_clusters %>% select(le))
shp_clusters$gstar_cluster <- lisa_clusters(gstar)
shp_clusters <- shp_clusters %>%
  mutate(gstar_cluster=factor(gstar_cluster, levels=0:4, labels=lisa_labels(gstar)))
table(shp_clusters$gstar_cluster)
baseline_clusters <- shp_clusters %>% 
  mutate(baseline_gstar_cluster=gstar_cluster) %>% 
  select(-gstar_cluster) %>% 
  st_drop_geometry()

# change in LE df 
shp_clusters <- right_join(shp_cz, results_cbsa_cz %>% get_bivarite(., gender_mw="Men", yr_35="year5", cbsa_cz="cz"), by=c("LM_Code"="id")) 
shp_clusters <- shp_clusters %>% filter(!LM_Code%in%"587")

queen_w <- queen_weights(shp_clusters)
gstar <- local_gstar(queen_w,shp_clusters %>% select(abs_dif))
shp_clusters$gstar_cluster <- lisa_clusters(gstar)
shp_clusters <- shp_clusters %>%
  mutate(gstar_cluster=factor(gstar_cluster, levels=0:4, labels=lisa_labels(gstar)))
table(shp_clusters$gstar_cluster)
diffLE_clusters <- shp_clusters %>% 
  mutate(diffLE_gstar_cluster=gstar_cluster) %>%
  select(-gstar_cluster) %>% 
  st_drop_geometry()

# combining both df's and droplevel unused factor levels (Undefined, Isolated)
df <- baseline_clusters %>% right_join(diffLE_clusters, by="LM_Code") %>% 
  select(LM_Code, baseline_gstar_cluster, diffLE_gstar_cluster) %>% 
  filter(!baseline_gstar_cluster%in%c("Undefined", "Isolated")) %>% 
  filter(!diffLE_gstar_cluster%in%c("Undefined", "Isolated"))
df$baseline_gstar_cluster <- droplevels(df$baseline_gstar_cluster)
df$diffLE_gstar_cluster <- droplevels(df$diffLE_gstar_cluster)
df <- df %>% mutate(baseline_gstar_cluster=case_when(baseline_gstar_cluster%in%"Low-Low" ~ "Low",
                                                     baseline_gstar_cluster%in%"High-High"~ "High",
                                                     baseline_gstar_cluster%in%"Not significant" ~ "Not significant"),
                    diffLE_gstar_cluster=case_when(diffLE_gstar_cluster%in%"Low-Low" ~ "Low",
                                                   diffLE_gstar_cluster%in%"High-High"~ "High",
                                                   diffLE_gstar_cluster%in%"Not significant" ~ "Not significant"))


# checking unique combos, should be 9
crossing(df$baseline_gstar_cluster, df$diffLE_gstar_cluster)
df %>% count(baseline_gstar_cluster, diffLE_gstar_cluster) 


df <- df %>% mutate(type=case_when(
  baseline_gstar_cluster=="Low" & diffLE_gstar_cluster=="Low" ~ "Low Baseline - Decreased",
  baseline_gstar_cluster=="Low" & diffLE_gstar_cluster=="Not significant" ~ "Low Baseline - NS Change",
  baseline_gstar_cluster=="Low" & diffLE_gstar_cluster=="High" ~ "Low Baseline - Increased",
  baseline_gstar_cluster=="Not significant" & diffLE_gstar_cluster=="Low" ~ "NS Baseline - Decreased",
  baseline_gstar_cluster=="Not significant" & diffLE_gstar_cluster=="Not significant" ~ "NS Baseline - NS Change",
  baseline_gstar_cluster=="Not significant" & diffLE_gstar_cluster=="High" ~ "NS Baseline - Increased",
  baseline_gstar_cluster=="High" & diffLE_gstar_cluster=="Low" ~ "High Baseline - Decreased",
  baseline_gstar_cluster=="High" & diffLE_gstar_cluster=="Not significant" ~ "High Baseline - NS Change",
  baseline_gstar_cluster=="High" & diffLE_gstar_cluster=="High" ~ "High Baseline - Increased"))


df <- df %>% mutate(type=factor(type,
                                levels=c("Low Baseline - Decreased",
                                         "Low Baseline - NS Change",
                                         "Low Baseline - Increased",
                                         "NS Baseline - Decreased",
                                         "NS Baseline - NS Change",
                                         "NS Baseline - Increased", 
                                         "High Baseline - Decreased",
                                         "High Baseline - NS Change",
                                         "High Baseline - Increased"))) %>% 
  arrange(type)
cluster_count <- df %>% count(type)

biscale_men_df <- right_join(shp_cz, df, by="LM_Code") %>% glimpse()

# CREATING MAP FOR BASELINE LE X CHANGE IN LE 
biscale_men_map <- ggplot()+
  geom_sf(data=biscale_men_df, size=0,color=NA,
          aes(geometry=geometry, fill=(type)))+
  geom_sf(data=biscale_men_df, size=.1,fill=NA,color="black",
          aes(geometry=geometry))+
  geom_sf(data=st_transform(df_state, crs = st_crs(shp_cz)), size=0.1, color="black", fill=NA)+
  geom_sf(data=st_transform(df_mexico, crs=st_crs(shp_cz)), size=0.1, color="black", fill="darkgrey")+
  geom_sf(data=st_transform(df_canada, crs=st_crs(shp_cz)), size=0.1, color="black", fill="darkgrey")+
  geom_sf(data=shp_cz %>% filter(LM_Code%in%"587") %>% select(geometry), size=0.1, color="black", fill="black")+
  geom_sf(data=st_transform(shp_census_region, crs = st_crs(shp_cz)), size=1.5, color="black", fill=NA)+
  geom_sf(data=st_transform(shp_census_division, crs = st_crs(shp_cz)), size=0.75, color="black", fill=NA)+
  annotate("text", x=st_bbox(biscale_men_df)$xmin, 
           y=st_bbox(biscale_men_df)$ymax, 
           hjust=0, vjust=0,
           size=10,
           label="A")+
  scale_fill_manual(values=cols)+
  coord_sf(xlim=st_bbox(biscale_men_df)[c("xmin", "xmax")],
           ylim=st_bbox(biscale_men_df)[c("ymin", "ymax")])+
  guides(color=guide_legend(reverse=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE))+
  #labs(title="Men")+
  #labs(tag="A")+
  map_theme+
  theme(legend.position="bottom", legend.title=element_blank())+
  theme(legend.background = element_blank(),
        legend.box.background = element_blank())


## WOMEN ----

### baseline LE ----
shp_clusters <- right_join(shp_cz, results_cbsa_cz %>%
                             filter(type%in%"cz",
                                    year_type%in%"year5",
                                    year%in%"1990-1994",
                                    gender%in%"Women"),
                           by=c("LM_Code"="id")) %>% arrange(gender, LM_Code)
shp_clusters <- shp_clusters %>% filter(!LM_Code%in%"587")

queen_w <- queen_weights(shp_clusters)
gstar <- local_gstar(queen_w,shp_clusters %>% select(le))
shp_clusters$gstar_cluster <- lisa_clusters(gstar)
shp_clusters <- shp_clusters %>%
  mutate(gstar_cluster=factor(gstar_cluster, levels=0:4, labels=lisa_labels(gstar)))
table(shp_clusters$gstar_cluster)

# CREATING MAP FOR GETIS ORD: BASELINE LE FOR WOMEN 
getisord_baselineLE_women <- ggplot()+
  geom_sf(data=shp_clusters %>% filter(!is.na(gstar_cluster)), size=0,color=NA,
          aes(geometry=geometry, fill=(gstar_cluster)))+
  geom_sf(data=shp_clusters, size=.1,fill=NA,color="black",
          aes(geometry=geometry))+
  geom_sf(data=st_transform(df_state, crs = st_crs(shp_cz)), size=0.1, color="black", fill=NA)+
  geom_sf(data=st_transform(df_mexico, crs=st_crs(shp_cz)), size=0.1, color="black", fill="darkgrey")+
  geom_sf(data=st_transform(df_canada, crs=st_crs(shp_cz)), size=0.1, color="black", fill="darkgrey")+
  geom_sf(data=shp_cz %>% filter(LM_Code%in%"587") %>% select(geometry), size=0.1, color="black", fill="black")+
  geom_sf(data=st_transform(shp_census_region, crs = st_crs(shp_cz)), size=1.5, color="black", fill=NA)+
  geom_sf(data=st_transform(shp_census_division, crs = st_crs(shp_cz)), size=0.75, color="black", fill=NA)+
  annotate("text", x=st_bbox(shp_clusters)$xmin, 
           y=st_bbox(shp_clusters)$ymax, 
           hjust=0, vjust=0,
           size=10,
           label="B")+
  scale_fill_manual(values=moran_colors, name="")+
  coord_sf(xlim=st_bbox(shp_clusters)[c("xmin", "xmax")],
           ylim=st_bbox(shp_clusters)[c("ymin", "ymax")])+
  guides(size=F, alpha=F, fill = guide_legend(override.aes = list(alpha=0))) +
  #labs(title=title) +
  #labs(title="", tag="B")+
  map_theme + 
  theme(legend.text=element_blank()) +
  theme(legend.position="none")


### change in LE ----

shp_clusters <- right_join(shp_cz, results_cbsa_cz %>% get_bivarite(., gender_mw="Women", yr_35="year5", cbsa_cz="cz"), by=c("LM_Code"="id"))
shp_clusters <- shp_clusters %>% filter(!LM_Code%in%"587")

queen_w <- queen_weights(shp_clusters)
gstar <- local_gstar(queen_w,shp_clusters %>% select(abs_dif))
shp_clusters$gstar_cluster <- lisa_clusters(gstar)
shp_clusters <- shp_clusters %>%
  mutate(gstar_cluster=factor(gstar_cluster, levels=0:4, labels=lisa_labels(gstar)))
table(shp_clusters$gstar_cluster)

# CREATING MAP FOR GETIS ORD: DIFFERENCE IN LIFE EXPECTANCY B/W 2015-2019 AND 1990-1994 FOR WOMEN 
getisord_diffLE_women <- ggplot()+
  geom_sf(data=shp_clusters %>% filter(!is.na(gstar_cluster)), size=0,color=NA,
          aes(geometry=geometry, fill=(gstar_cluster)))+
  geom_sf(data=shp_clusters, size=.1,fill=NA,color="black",
          aes(geometry=geometry))+
  geom_sf(data=st_transform(df_state, crs = st_crs(shp_cz)), size=0.1, color="black", fill=NA)+
  geom_sf(data=st_transform(df_mexico, crs=st_crs(shp_cz)), size=0.1, color="black", fill="darkgrey")+
  geom_sf(data=st_transform(df_canada, crs=st_crs(shp_cz)), size=0.1, color="black", fill="darkgrey")+
  geom_sf(data=shp_cz %>% filter(LM_Code%in%"587") %>% select(geometry), size=0.1, color="black", fill="black")+
  geom_sf(data=st_transform(shp_census_region, crs = st_crs(shp_cz)), size=1.5, color="black", fill=NA)+
  geom_sf(data=st_transform(shp_census_division, crs = st_crs(shp_cz)), size=0.75, color="black", fill=NA)+
  annotate("text", x=st_bbox(shp_clusters)$xmin, 
           y=st_bbox(shp_clusters)$ymax, 
           hjust=0, vjust=0,
           size=10,
           label="D")+
  scale_fill_manual(values=moran_colors, name="")+
  coord_sf(xlim=st_bbox(shp_clusters)[c("xmin", "xmax")],
           ylim=st_bbox(shp_clusters)[c("ymin", "ymax")])+
  guides(size=F, alpha=F, fill = guide_legend(override.aes = list(alpha=0))) +
  #labs(title=title) +
  #labs(title="", tag="D")+
  map_theme + 
  theme(legend.text=element_blank()) +
  theme(legend.position="none")

### biscale 3x3 plot: baseline LE & change in LE ----

# baseline df
shp_clusters <- right_join(shp_cz, results_cbsa_cz %>%
                             filter(type%in%"cz",
                                    year_type%in%"year5",
                                    year%in%"1990-1994",
                                    gender%in%"Women"),
                           by=c("LM_Code"="id")) %>% arrange(gender, LM_Code)
shp_clusters <- shp_clusters %>% filter(!LM_Code%in%"587")

queen_w <- queen_weights(shp_clusters)
gstar <- local_gstar(queen_w,shp_clusters %>% select(le))
shp_clusters$gstar_cluster <- lisa_clusters(gstar)
shp_clusters <- shp_clusters %>%
  mutate(gstar_cluster=factor(gstar_cluster, levels=0:4, labels=lisa_labels(gstar)))
table(shp_clusters$gstar_cluster)
baseline_clusters <- shp_clusters %>% 
  mutate(baseline_gstar_cluster=gstar_cluster) %>% 
  select(-gstar_cluster) %>% 
  st_drop_geometry()

# change in LE df 
shp_clusters <- right_join(shp_cz, results_cbsa_cz %>% get_bivarite(., gender_mw="Women", yr_35="year5", cbsa_cz="cz"), by=c("LM_Code"="id")) 
shp_clusters <- shp_clusters %>% filter(!LM_Code%in%"587")

queen_w <- queen_weights(shp_clusters)
gstar <- local_gstar(queen_w,shp_clusters %>% select(abs_dif))
shp_clusters$gstar_cluster <- lisa_clusters(gstar)
shp_clusters <- shp_clusters %>%
  mutate(gstar_cluster=factor(gstar_cluster, levels=0:4, labels=lisa_labels(gstar)))
table(shp_clusters$gstar_cluster)
diffLE_clusters <- shp_clusters %>% 
  mutate(diffLE_gstar_cluster=gstar_cluster) %>%
  select(-gstar_cluster) %>% 
  st_drop_geometry()

# combining both df's and droplevel unused factor levels (Undefined, Isolated)
df <- baseline_clusters %>% right_join(diffLE_clusters, by="LM_Code") %>% 
  select(LM_Code, baseline_gstar_cluster, diffLE_gstar_cluster) %>% 
  filter(!baseline_gstar_cluster%in%c("Undefined", "Isolated")) %>% 
  filter(!diffLE_gstar_cluster%in%c("Undefined", "Isolated"))
df$baseline_gstar_cluster <- droplevels(df$baseline_gstar_cluster)
df$diffLE_gstar_cluster <- droplevels(df$diffLE_gstar_cluster)
df <- df %>% mutate(baseline_gstar_cluster=case_when(baseline_gstar_cluster%in%"Low-Low" ~ "Low",
                                                     baseline_gstar_cluster%in%"High-High"~ "High",
                                                     baseline_gstar_cluster%in%"Not significant" ~ "Not significant"),
                    diffLE_gstar_cluster=case_when(diffLE_gstar_cluster%in%"Low-Low" ~ "Low",
                                                   diffLE_gstar_cluster%in%"High-High"~ "High",
                                                   diffLE_gstar_cluster%in%"Not significant" ~ "Not significant"))


# checking unique combos, should be 9
crossing(df$baseline_gstar_cluster, df$diffLE_gstar_cluster)
df %>% count(baseline_gstar_cluster, diffLE_gstar_cluster) 


df <- df %>% mutate(type=case_when(
  baseline_gstar_cluster=="Low" & diffLE_gstar_cluster=="Low" ~ "Low Baseline - Decreased",
  baseline_gstar_cluster=="Low" & diffLE_gstar_cluster=="Not significant" ~ "Low Baseline - NS Change",
  baseline_gstar_cluster=="Low" & diffLE_gstar_cluster=="High" ~ "Low Baseline - Increased",
  baseline_gstar_cluster=="Not significant" & diffLE_gstar_cluster=="Low" ~ "NS Baseline - Decreased",
  baseline_gstar_cluster=="Not significant" & diffLE_gstar_cluster=="Not significant" ~ "NS Baseline - NS Change",
  baseline_gstar_cluster=="Not significant" & diffLE_gstar_cluster=="High" ~ "NS Baseline - Increased",
  baseline_gstar_cluster=="High" & diffLE_gstar_cluster=="Low" ~ "High Baseline - Decreased",
  baseline_gstar_cluster=="High" & diffLE_gstar_cluster=="Not significant" ~ "High Baseline - NS Change",
  baseline_gstar_cluster=="High" & diffLE_gstar_cluster=="High" ~ "High Baseline - Increased"))


df <- df %>% mutate(type=factor(type,
                                levels=c("Low Baseline - Decreased",
                                         "Low Baseline - NS Change",
                                         "Low Baseline - Increased",
                                         "NS Baseline - Decreased",
                                         "NS Baseline - NS Change",
                                         "NS Baseline - Increased", 
                                         "High Baseline - Decreased",
                                         "High Baseline - NS Change",
                                         "High Baseline - Increased"))) %>% 
  arrange(type)
cluster_count <- df %>% count(type)

biscale_women_df <- right_join(shp_cz, df, by="LM_Code") 

# CREATING MAP FOR BASELINE LE X CHANGE IN LE 
biscale_women_map <- ggplot()+
  geom_sf(data=biscale_women_df, size=0,color=NA,
          aes(geometry=geometry, fill=(type)))+
  geom_sf(data=biscale_women_df, size=.1,fill=NA,color="black",
          aes(geometry=geometry))+
  geom_sf(data=st_transform(df_state, crs = st_crs(shp_cz)), size=0.1, color="black", fill=NA)+
  geom_sf(data=st_transform(df_mexico, crs=st_crs(shp_cz)), size=0.1, color="black", fill="darkgrey")+
  geom_sf(data=st_transform(df_canada, crs=st_crs(shp_cz)), size=0.1, color="black", fill="darkgrey")+
  geom_sf(data=shp_cz %>% filter(LM_Code%in%"587") %>% select(geometry), size=0.1, color="black", fill="black")+
  geom_sf(data=st_transform(shp_census_region, crs = st_crs(shp_cz)), size=1.5, color="black", fill=NA)+
  geom_sf(data=st_transform(shp_census_division, crs = st_crs(shp_cz)), size=0.75, color="black", fill=NA)+
  annotate("text", x=st_bbox(biscale_men_df)$xmin, 
           y=st_bbox(biscale_men_df)$ymax, 
           hjust=0, vjust=0,
           size=10,
           label="B")+
  scale_fill_manual(values=cols)+
  coord_sf(xlim=st_bbox(biscale_women_df)[c("xmin", "xmax")],
           ylim=st_bbox(biscale_women_df)[c("ymin", "ymax")])+
  guides(color=guide_legend(reverse=TRUE),
         fill=guide_legend(nrow=3, byrow=TRUE))+
  #labs(subtitle="Women")+
  #labs(tag="B")+
  map_theme+
  theme(legend.position="bottom", legend.title=element_blank())+
  theme(legend.background = element_blank(),
        legend.box.background = element_blank())



pall<-arrangeGrob(grobs=list(getisord_baselineLE_men, 
                             getisord_baselineLE_women, 
                             getisord_diffLE_men, 
                             getisord_diffLE_women), 
                  ncol=2)
pall<-arrangeGrob(grobs=list(pall, legend_sep), heights=c(20, 1), ncol=1)
ggsave("../Tables & Figures/Figure3_sep.pdf", pall, width=29, height=19)


legend<-get_legend(biscale_men_map)
biscale_men_map<-biscale_men_map+guides(fill="none") 
biscale_women_map<-biscale_women_map+guides(fill="none")
map_all<-arrangeGrob(grobs=list(biscale_men_map, biscale_women_map), ncol=2)
map_all_forreal<-arrangeGrob(grobs=list(map_all, legend), ncol=1, heights=c(15, 2))
ggsave("../Tables & Figures/Figure3_combined.pdf", map_all_forreal, width=30, height=11)


## OLD CODE:



## MEN ----

### baseline LE ----
shp_clusters <- right_join(shp_cz, results_cbsa_cz %>% 
                             filter(type%in%"cz",
                                    year_type%in%"year5", 
                                    year%in%"1990-1994", 
                                    gender%in%"Men"), 
                           by=c("LM_Code"="id")) %>% arrange(gender, LM_Code)

neighbors <- poly2nb(shp_clusters %>% group_by(LM_Code))
# fixing San Juan County, WA, CZ 604 (island, but manually input King County CZ 602 as neighbor)
neighbors[[which(shp_clusters$LM_Code==604)]] <- which(shp_clusters$LM_Code%in%c(602))
neighbors[[which(shp_clusters$LM_Code==602)]] <- c(neighbors[[which(shp_clusters$LM_Code==602)]], which(shp_clusters$LM_Code%in%c(604)))

# it goes by indice in shp_clusters dataframe! for example: 
# neighbors[[601]] is shp_clusters[601,] == CZ 92 (Philadelphia County, PA) 
# neighbors of CZ 92 (Philly County) are 
# 182 == shp_clusters[182,] == CZ 275 (Baltimore County, MD)
# 307 == shp_clusters[307,] == CZ 390 (Camden County, NJ)
# 308 == shp_clusters[308,] == CZ 391 (Bergen County, NJ)
# 415 == shp_clusters[415,] == CZ 488 (Lancaster County, PA)
# 418 == shp_clusters[418,] == CZ 490 (Berks County, PA)
# 600 == shp_clusters[600,] == CZ 91  (Sussex County, DE)

# use the moran.mc function
temp <- moran.mc(x = shp_clusters %>% pull(le), 
                 # listw = list of neighbors, obtained from the nb adjacency matrix by using the nb2listw spdep function
                 listw = nb2listw(neighbors, style = "B"), 
                 # last, number of permutations for the permutation-based moran test
                 nsim = 9999)
# extract moran statistic and p value
globalmoran <- data.frame(statistic=temp$statistic[1], pval=temp$p.value) %>% 
  mutate(out=paste0(round(statistic, digits=3), " (", ifelse(pval<0.001, "<0.001", round(pval, digits=3)), ")")) 
globalmoran
## BASELINE LE FOR MEN: 0.622 (<0.001)


### change in LE ----
shp_clusters <- right_join(shp_cz, results_cbsa_cz %>% get_bivarite(., gender_mw="Men", yr_35="year5", cbsa_cz="cz"), by=c("LM_Code"="id"))

neighbors <- poly2nb(shp_clusters %>% group_by(LM_Code))
# fixing San Juan County, WA, CZ 604 (island, but manually input King County CZ 602 as neighbor)
neighbors[[which(shp_clusters$LM_Code==604)]] <- which(shp_clusters$LM_Code%in%c(602))
neighbors[[which(shp_clusters$LM_Code==602)]] <- c(neighbors[[which(shp_clusters$LM_Code==602)]], which(shp_clusters$LM_Code%in%c(604)))

# use the moran.mc function
temp <- moran.mc(x = shp_clusters %>% pull(abs_dif), 
                 # listw = list of neighbors, obtained from the nb adjacency matrix by using the nb2listw spdep function
                 listw = nb2listw(neighbors, style = "B"), 
                 # last, number of permutations for the permutation-based moran test
                 nsim = 9999)
# extract moran statistic and p value
globalmoran <- data.frame(statistic=temp$statistic[1], pval=temp$p.value) %>% 
  mutate(out=paste0(round(statistic, digits=3), " (", ifelse(pval<0.001, "<0.001", round(pval, digits=3)), ")")) 
globalmoran
## DIFFERENCE LE (FROM 1990-1994 TO 2015-2019) FOR MEN: 0.313 (<0.001)


## WOMEN ----

### baseline LE ----
shp_clusters <- right_join(shp_cz, results_cbsa_cz %>% 
                             filter(type%in%"cz",
                                    year_type%in%"year5", 
                                    year%in%"1990-1994", 
                                    gender%in%"Women"), 
                           by=c("LM_Code"="id")) %>% arrange(gender, LM_Code)

neighbors <- poly2nb(shp_clusters %>% group_by(LM_Code))
# fixing San Juan County, WA, CZ 604 (island, but manually input King County CZ 602 as neighbor)
neighbors[[which(shp_clusters$LM_Code==604)]] <- which(shp_clusters$LM_Code%in%c(602))
neighbors[[which(shp_clusters$LM_Code==602)]] <- c(neighbors[[which(shp_clusters$LM_Code==602)]], which(shp_clusters$LM_Code%in%c(604)))

# use the moran.mc function
temp <- moran.mc(x = shp_clusters %>% pull(le), 
                 # listw = list of neighbors, obtained from the nb adjacency matrix by using the nb2listw spdep function
                 listw = nb2listw(neighbors, style = "B"), 
                 # last, number of permutations for the permutation-based moran test
                 nsim = 9999)
# extract moran statistic and p value
globalmoran <- data.frame(statistic=temp$statistic[1], pval=temp$p.value) %>% 
  mutate(out=paste0(round(statistic, digits=3), " (", ifelse(pval<0.001, "<0.001", round(pval, digits=3)), ")")) 
globalmoran
## BASELINE LE FOR WOMEN: 0.614 (<0.001)



### change in LE ----
shp_clusters <- right_join(shp_cz, results_cbsa_cz %>% get_bivarite(., gender_mw="Women", yr_35="year5", cbsa_cz="cz"), by=c("LM_Code"="id"))

neighbors <- poly2nb(shp_clusters %>% group_by(LM_Code))
# fixing San Juan County, WA, CZ 604 (island, but manually input King County CZ 602 as neighbor)
neighbors[[which(shp_clusters$LM_Code==604)]] <- which(shp_clusters$LM_Code%in%c(602))
neighbors[[which(shp_clusters$LM_Code==602)]] <- c(neighbors[[which(shp_clusters$LM_Code==602)]], which(shp_clusters$LM_Code%in%c(604)))

# use the moran.mc function
temp <- moran.mc(x = shp_clusters %>% pull(abs_dif), 
                 # listw = list of neighbors, obtained from the nb adjacency matrix by using the nb2listw spdep function
                 listw = nb2listw(neighbors, style = "B"), 
                 # last, number of permutations for the permutation-based moran test
                 nsim = 9999)
# extract moran statistic and p value
globalmoran <- data.frame(statistic=temp$statistic[1], pval=temp$p.value) %>% 
  mutate(out=paste0(round(statistic, digits=3), " (", ifelse(pval<0.001, "<0.001", round(pval, digits=3)), ")")) 
globalmoran
## DIFFERENCE LE (FROM 1990-1994 TO 2015-2019) FOR WOMEN: 0.355 (<0.001)

