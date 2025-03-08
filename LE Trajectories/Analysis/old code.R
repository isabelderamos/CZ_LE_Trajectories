


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





## testing resampled getis ord

getis_ord_fun_iter<-function(shp, data, gender_select, significance, correct="none"){
  #gender_select<-"Men";significance<-0.05; correct="none"
  #shp<-shp_cz; data<-results_czyr5_iter %>% filter(iter==1)
  shp_clusters <- right_join(shp, 
                             data %>% 
                               get_bivarite(., gender_mw=gender_select, 
                                            yr_35="year5", cbsa_cz="cz"), by=c("LM_Code"="id")) 
  shp_clusters <- shp_clusters %>% filter(!LM_Code%in%"587")
  queen_w <- queen_weights(shp_clusters)
  gstar_baseline <- local_gstar(queen_w,shp_clusters %>% select(le), significance_cutoff = 1) 
  #significance_cutoff = significance, permutations = 99999)
  gstar_change <- local_gstar(queen_w,shp_clusters %>% select(abs_dif), significance_cutoff = 1)
  #significance_cutoff = significance, permutations = 99999)
  shp_clusters$diffLE_gstar_cluster_pval<-lisa_pvalues(gstar_change)
  shp_clusters$baseline_gstar_cluster_pval<-lisa_pvalues(gstar_baseline)
  if (correct=="none"){
    shp_clusters$baseline_gstar_cluster <- lisa_clusters(gstar_baseline)
    shp_clusters$diffLE_gstar_cluster <- lisa_clusters(gstar_change)
  } else if (correct=="fdr"){
    fdr_baseline<-lisa_fdr(gstar_baseline, current_p = significance)
    fdr_change<-lisa_fdr(gstar_change, current_p = significance)
    shp_clusters$baseline_gstar_cluster <- lisa_clusters(gstar_baseline, cutoff = fdr_baseline)
    shp_clusters$diffLE_gstar_cluster <- lisa_clusters(gstar_change, cutoff=fdr_change)
  } else if (correct=="bonf"){
    bonf_baseline<-lisa_bo(gstar_baseline, current_p = significance)
    bonf_change<-lisa_fdr(gstar_change, current_p = significance)
    shp_clusters$baseline_gstar_cluster <- lisa_clusters(gstar_baseline, cutoff = bonf_baseline)
    shp_clusters$diffLE_gstar_cluster <- lisa_clusters(gstar_change, cutoff=bonf_change)
  }
  shp_clusters <- shp_clusters %>%
    mutate(baseline_gstar_cluster=factor(baseline_gstar_cluster, levels=0:4, labels=lisa_labels(gstar_baseline)),
           diffLE_gstar_cluster=factor(diffLE_gstar_cluster, levels=0:4, labels=lisa_labels(gstar_change)),
           gender=gender_select,
           significance=significance)
  print(table(shp_clusters$baseline_gstar_cluster, shp_clusters$diffLE_gstar_cluster))
  
  # combining both df's and droplevel unused factor levels (Undefined, Isolated)
  df <- shp_clusters %>% st_drop_geometry %>% 
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
  print(df %>% count(type))
  
  # now that everyhting is checked, merge back into the original one
  shp_clusters<-full_join(shp_clusters, df %>% select(LM_Code, type))
  
  return(shp_clusters)
}

test<-results_czyr5_iter %>% 
  group_by(iter) %>% 
  group_modify(~{
    getis_ord_fun_iter(shp = shp_cz, data = .x, gender_select="Men", significance=0.05, correct = "none")
  })
test<-test %>% 
  filter(baseline_gstar_cluster!="Isolated") %>% 
  mutate(c_high=ifelse(baseline_gstar_cluster=="High-High", 1, 0),
         c_low=ifelse(baseline_gstar_cluster=="Low-Low", 1, 0)) %>% 
  group_by(LM_Code) %>% 
  summarise(prob_high=mean(c_high),
            prob_low=mean(c_low))



all_men<-getis_ord_fun(gender_select="Men", significance=0.05, correct = "none")
all_women<-getis_ord_fun(gender_select="Women", significance=0.05, correct = "none")
