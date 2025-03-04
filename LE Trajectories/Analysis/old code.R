


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

