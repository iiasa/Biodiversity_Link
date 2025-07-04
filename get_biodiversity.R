args <- commandArgs(trailing = TRUE)


library(tidyverse)
library(fs)

# Define_functions ------

#' Auxiliary function to recode BII coefficients
#' according to forest and non-forest ecosystems
recode_if <- function(x, condition, ...) {
  ifelse(condition, recode(x, ...), x)
}

#' Auxiliary function to calculate the
#' amount of restoration/abandonment per period
calc_lag <- function(x,n){
  k <- length(x)
  y <- rep(NA,k)
  y[1:n] <- 0
  for (i in (n+1):k){
    y[i] <- x[i] - x[i-n]
  }
  return(y)
}

#' Auxiliary function to rescale the amount
#' of abandoned land per period, according
#' to the total availability of OthNatLnd
rebalance <- function(x){

  n <- length(x)

  tot_abn <- sum(x[2:(n-1)])

  if(tot_abn > x[n]){

    shares <- x[2:(n-1)] / tot_abn
    delta <- x[n] * shares
    x[2:(n-1)] <- x[2:(n-1)] - delta

  } else {
    dif <- x[n] - tot_abn
    x[2:(n-1)] <- 0
    #x[1] <- x[1] - dif
  }

  return(x)
}

#' Auxiliary function to compute the
#' age classes of abandoned land
get_age <- function(x,lnd){

  x <- x %>% spread(times,value) %>% ungroup() %>% dplyr::select(-lu.to)
  if (dim(x)[1] > 0) x <- x %>% mutate('2000'=0,'2010'=0) %>% relocate(c('2000','2010'),.before = '2020')

  if (dim(x)[1] > 0){
    # Get restoration from 0 to 10 years
    ini <- 3
    fin <- dim(x)[2]
    restored_10 <- x
    restored_10[,ini:fin] <- apply(x[,ini:fin],1,function(x) calc_lag(x,1)) %>% t()

    # Get remaining restoration classes
    restored_aux <- x
    restored_aux[,ini:fin] <- 0
    restored_20 <- restored_30 <- restored_40 <- restored_50 <-
      restored_60 <- restored_70 <- restored_80 <- restored_aux

    # Move afforestation along the time periods to the next age class
    if ((ini+1) <= fin) restored_20[,(ini+1):fin] <- restored_10[,ini:(fin-1)]
    if ((ini+2) <= fin) restored_30[,(ini+2):fin] <- restored_10[,ini:(fin-2)]
    if ((ini+3) <= fin) restored_40[,(ini+3):fin] <- restored_10[,ini:(fin-3)]
    if ((ini+4) <= fin) restored_50[,(ini+4):fin] <- restored_10[,ini:(fin-4)]
    if ((ini+5) <= fin) restored_60[,(ini+5):fin] <- restored_10[,ini:(fin-5)]
    if ((ini+6) <= fin) restored_70[,(ini+6):fin] <- restored_10[,ini:(fin-6)]
    if ((ini+7) <= fin) restored_80[,ini+7] <- restored_10[,ini]
    if ((ini+8) <= fin) restored_80[,ini+8] <- restored_10[,ini+1] + restored_10[,ini]
    if ((ini+9) <= fin) restored_80[,ini+9] <- restored_10[,ini+2] + restored_10[,ini+1] + restored_10[,ini]
    if ((ini+10) <= fin) restored_80[,ini+10] <-   restored_10[,ini+3] + restored_10[,ini+2] + restored_10[,ini+1] + restored_10[,ini]

    # Combined the different age classes
    restoration <- restored_10 %>% gather(times,value,-c(REGION,ns)) %>% rename(restored_10=value) %>%
      left_join(restored_20 %>% gather(times,value,-c(REGION,ns)) %>% rename(restored_20=value)) %>%
      left_join(restored_30 %>% gather(times,value,-c(REGION,ns)) %>% rename(restored_30=value)) %>%
      left_join(restored_40 %>% gather(times,value,-c(REGION,ns)) %>% rename(restored_40=value)) %>%
      left_join(restored_50 %>% gather(times,value,-c(REGION,ns)) %>% rename(restored_50=value)) %>%
      left_join(restored_60 %>% gather(times,value,-c(REGION,ns)) %>% rename(restored_60=value)) %>%
      left_join(restored_70 %>% gather(times,value,-c(REGION,ns)) %>% rename(restored_70=value)) %>%
      left_join(restored_80 %>% gather(times,value,-c(REGION,ns)) %>% rename(restored_80=value)) %>%
      mutate(times=as.integer(times))

    colnames(restoration)[4:11] <- str_c(lnd,"_",colnames(restoration)[4:11])

    # Aggregate and clean
    restoration <- restoration %>% gather(lu.to,value,-c( REGION,times,ns))

  } else {
    restoration <- tibble(REGION=NA, ns=NA, times=NA, lu.to=NA, value = NA)
  }

  return(restoration)

}

#' Function to add BII coefficients and derive impacts
process_bii <- function(x){

  out <- x %>%
    left_join(lc_map) %>% rename(LUclass=LC_BTC) %>% gather(pnv,f_nf_share,-c(REGION,times,ns,LC_TYPE,LUclass,value,pot_npp)) %>%
    mutate(LUclass = recode_if(LUclass, pnv == "forest_share", "cropland_other" = "cropland_other_f"),
           LUclass = recode_if(LUclass, pnv == "nonforest_share", "cropland_other" = "cropland_other_nf"),
           LUclass = recode_if(LUclass, pnv == "forest_share", "built.up" = "built.up_f"),
           LUclass = recode_if(LUclass, pnv == "nonforest_share", "built.up" = "built.up_nf"),
           LUclass = recode_if(LUclass, pnv == "forest_share", "cropland_2Gbioen" = "cropland_2Gbioen_f"),
           LUclass = recode_if(LUclass, pnv == "nonforest_share", "cropland_2Gbioen" = "cropland_2Gbioen_nf"),
           LUclass = recode_if(LUclass, pnv == "forest_share", "grassland" = "grassland_f"),
           LUclass = recode_if(LUclass, pnv == "nonforest_share", "grassland" = "grassland_nf"),
           LUclass = recode_if(LUclass, pnv == "forest_share", "other" = "other_f"),
           LUclass = recode_if(LUclass, pnv == "nonforest_share", "other" = "other_nf"),
           LUclass = recode_if(LUclass, pnv == "forest_share", "AbnCrpLnd_restored_10" = "AbnCrpLnd_restored_10_f"),
           LUclass = recode_if(LUclass, pnv == "nonforest_share", "AbnCrpLnd_restored_10" = "AbnCrpLnd_restored_10_nf"),
           LUclass = recode_if(LUclass, pnv == "forest_share", "AbnCrpLnd_restored_20" = "AbnCrpLnd_restored_20_f"),
           LUclass = recode_if(LUclass, pnv == "nonforest_share", "AbnCrpLnd_restored_20" = "AbnCrpLnd_restored_20_nf"),
           LUclass = recode_if(LUclass, pnv == "forest_share", "AbnCrpLnd_restored_30" = "AbnCrpLnd_restored_30_f"),
           LUclass = recode_if(LUclass, pnv == "nonforest_share", "AbnCrpLnd_restored_30" = "AbnCrpLnd_restored_30_nf"),
           LUclass = recode_if(LUclass, pnv == "forest_share", "AbnCrpLnd_restored_40" = "AbnCrpLnd_restored_40_f"),
           LUclass = recode_if(LUclass, pnv == "nonforest_share", "AbnCrpLnd_restored_40" = "AbnCrpLnd_restored_40_nf"),
           LUclass = recode_if(LUclass, pnv == "forest_share", "AbnCrpLnd_restored_50" = "AbnCrpLnd_restored_50_f"),
           LUclass = recode_if(LUclass, pnv == "nonforest_share", "AbnCrpLnd_restored_50" = "AbnCrpLnd_restored_50_nf"),
           LUclass = recode_if(LUclass, pnv == "forest_share", "AbnCrpLnd_restored_60" = "AbnCrpLnd_restored_60_f"),
           LUclass = recode_if(LUclass, pnv == "nonforest_share", "AbnCrpLnd_restored_60" = "AbnCrpLnd_restored_60_nf"),
           LUclass = recode_if(LUclass, pnv == "forest_share", "AbnCrpLnd_restored_70" = "AbnCrpLnd_restored_70_f"),
           LUclass = recode_if(LUclass, pnv == "nonforest_share", "AbnCrpLnd_restored_70" = "AbnCrpLnd_restored_70_nf"),
           LUclass = recode_if(LUclass, pnv == "forest_share", "AbnCrpLnd_restored_80" = "AbnCrpLnd_restored_80_f"),
           LUclass = recode_if(LUclass, pnv == "nonforest_share", "AbnCrpLnd_restored_80" = "AbnCrpLnd_restored_80_nf"),
           LUclass = recode_if(LUclass, pnv == "forest_share", "AbnGrsLnd_restored_10" = "AbnGrsLnd_restored_10_f"),
           LUclass = recode_if(LUclass, pnv == "nonforest_share", "AbnGrsLnd_restored_10" = "AbnGrsLnd_restored_10_nf"),
           LUclass = recode_if(LUclass, pnv == "forest_share", "AbnGrsLnd_restored_20" = "AbnGrsLnd_restored_20_f"),
           LUclass = recode_if(LUclass, pnv == "nonforest_share", "AbnGrsLnd_restored_20" = "AbnGrsLnd_restored_20_nf"),
           LUclass = recode_if(LUclass, pnv == "forest_share", "AbnGrsLnd_restored_30" = "AbnGrsLnd_restored_30_f"),
           LUclass = recode_if(LUclass, pnv == "nonforest_share", "AbnGrsLnd_restored_30" = "AbnGrsLnd_restored_30_nf"),
           LUclass = recode_if(LUclass, pnv == "forest_share", "AbnGrsLnd_restored_40" = "AbnGrsLnd_restored_40_f"),
           LUclass = recode_if(LUclass, pnv == "nonforest_share", "AbnGrsLnd_restored_40" = "AbnGrsLnd_restored_40_nf"),
           LUclass = recode_if(LUclass, pnv == "forest_share", "AbnGrsLnd_restored_50" = "AbnGrsLnd_restored_50_f"),
           LUclass = recode_if(LUclass, pnv == "nonforest_share", "AbnGrsLnd_restored_50" = "AbnGrsLnd_restored_50_nf"),
           LUclass = recode_if(LUclass, pnv == "forest_share", "AbnGrsLnd_restored_60" = "AbnGrsLnd_restored_60_f"),
           LUclass = recode_if(LUclass, pnv == "nonforest_share", "AbnGrsLnd_restored_60" = "AbnGrsLnd_restored_60_nf"),
           LUclass = recode_if(LUclass, pnv == "forest_share", "AbnGrsLnd_restored_70" = "AbnGrsLnd_restored_70_f"),
           LUclass = recode_if(LUclass, pnv == "nonforest_share", "AbnGrsLnd_restored_70" = "AbnGrsLnd_restored_70_nf"),
           LUclass = recode_if(LUclass, pnv == "forest_share", "AbnGrsLnd_restored_80" = "AbnGrsLnd_restored_80_f"),
           LUclass = recode_if(LUclass, pnv == "nonforest_share", "AbnGrsLnd_restored_80" = "AbnGrsLnd_restored_80_nf"),
           LUclass = recode_if(LUclass, pnv == "forest_share", "AbnPltFor_restored_10" = "AbnPltFor_restored_10_f"),
           LUclass = recode_if(LUclass, pnv == "nonforest_share", "AbnPltFor_restored_10" = "AbnPltFor_restored_10_nf"),
           LUclass = recode_if(LUclass, pnv == "forest_share", "AbnPltFor_restored_20" = "AbnPltFor_restored_20_f"),
           LUclass = recode_if(LUclass, pnv == "nonforest_share", "AbnPltFor_restored_20" = "AbnPltFor_restored_20_nf"),
           LUclass = recode_if(LUclass, pnv == "forest_share", "AbnPltFor_restored_30" = "AbnPltFor_restored_30_f"),
           LUclass = recode_if(LUclass, pnv == "nonforest_share", "AbnPltFor_restored_30" = "AbnPltFor_restored_30_nf"),
           LUclass = recode_if(LUclass, pnv == "forest_share", "AbnPltFor_restored_40" = "AbnPltFor_restored_40_f"),
           LUclass = recode_if(LUclass, pnv == "nonforest_share", "AbnPltFor_restored_40" = "AbnPltFor_restored_40_nf"),
           LUclass = recode_if(LUclass, pnv == "forest_share", "AbnPltFor_restored_50" = "AbnPltFor_restored_50_f"),
           LUclass = recode_if(LUclass, pnv == "nonforest_share", "AbnPltFor_restored_50" = "AbnPltFor_restored_50_nf"),
           LUclass = recode_if(LUclass, pnv == "forest_share", "AbnPltFor_restored_60" = "AbnPltFor_restored_60_f"),
           LUclass = recode_if(LUclass, pnv == "nonforest_share", "AbnPltFor_restored_60" = "AbnPltFor_restored_60_nf"),
           LUclass = recode_if(LUclass, pnv == "forest_share", "AbnPltFor_restored_70" = "AbnPltFor_restored_70_f"),
           LUclass = recode_if(LUclass, pnv == "nonforest_share", "AbnPltFor_restored_70" = "AbnPltFor_restored_70_nf"),
           LUclass = recode_if(LUclass, pnv == "forest_share", "AbnPltFor_restored_80" = "AbnPltFor_restored_80_f"),
           LUclass = recode_if(LUclass, pnv == "nonforest_share", "AbnPltFor_restored_80" = "AbnPltFor_restored_80_nf")) %>%
    filter(! LC_TYPE %in% c( "forest_new_ha","forest_old_ha")) %>% left_join(bii_coefs) %>%
    mutate(value=ifelse(is.na(value),0,value)) %>% filter(value > 0, f_nf_share > 0)

  return(out)
}

# Function to add cSAR coefficients  and derive impacts
process_csar <- function(x){
  out <- x %>%
    left_join(lc_map) %>% left_join(csar_lc_map) %>%
    filter(! LC_TYPE %in% c( "forest_new_ha","forest_old_ha")) %>% left_join(eco_map) %>%
    left_join(csar_coefs) %>%
    mutate(value=ifelse(is.na(value),0,value)) %>% filter(value > 0) %>%
    mutate(CF=CF/10^3) %>%
    mutate(extinctions=value*CF*simu_share_in_eco) %>% drop_na()

}

# Get scenario------
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied", call.=FALSE)
} else {
  # default output file
  scen <- as.integer(args[1])
}


# Load G4M output conversion function------
source("G4M_DS_to_simu_link_final_limpopo.R")

# Get cluster number of downscaling post-processing------
config <- readRDS(path("Input","config.RData"))

project <- config[[1]]
lab <- config[[2]]
scen_map <- config[[3]]
cluster_nr <- config[[4]]
get_bii <- config[[5]]
get_csar <- config[[6]]


# Define loading function
rfunc <- function(x) readRDS(x)[[3]]$out.res

# Define scenarios
current_scen <- scen_map %>% filter(ScenNr == scen)

# Get scenario configuration

scen1 <- current_scen$SCEN1
scen2 <- current_scen$SCEN2
scen3 <- current_scen$SCEN3

# Load in downscalr output
downscalr_out0 <- rfunc(path(str_glue( "output_", cluster_nr, ".",
                                      sprintf("%06d",current_scen$ScenNr),".RData")))
#  group_by(times) %>% group_split()

# Process RstLnd in res---
cat("==> Process RstLnd in res","\n")

# mapping: ns to g4m 0.5degree grids
mapping <- readRDS(file='G4m_mapping.RData')[[1]]
mapping <- data.frame(apply(mapping, 2, as.numeric))

# mapping_g4mid_xy_new <- readRDS(file="mapping_for_G4MDSlink.RData")[[1]]
maplayer_isforest_new <- readRDS(file="mapping_for_G4MDSlink.RData")[[2]]

res0 <- downscalr_out0
mapping$SimUID <- as.character(mapping$SimUID)
mapping$g4m_05_id <- as.character(mapping$g4m_05_id)

# 1) mapping dat to g4mid
res1 <- res0 %>% left_join(mapping, by=c("ns"="SimUID"))

# Assign RstLnd to afforestable or non-afforestable
res2 <- res1 %>% left_join(maplayer_isforest_new) %>%
  mutate(IsForest=ifelse(is.na(IsForest),0,IsForest)) %>%
  mutate(lu.from=ifelse(IsForest==1,
                        recode(lu.from,"RstLnd"="RstLnd_YesAffor"),
                        recode(lu.from,"RstLnd"="RstLnd_NoAffor")))%>%
  mutate(lu.to=ifelse(IsForest==1,
                      recode(lu.to,"RstLnd"="RstLnd_YesAffor"),
                      recode(lu.to,"RstLnd"="RstLnd_NoAffor")))

# Reallocate Rstlnd_YesAffor to OtherNatLnd
res3 <- res2 %>%
  mutate(lu.from=recode(lu.from,"RstLnd_YesAffor"="OthNatLnd")) %>%
  mutate(lu.to=recode(lu.to,"RstLnd_YesAffor"="OthNatLnd"))

res4 <- res3 %>%
  group_by(REGION,times,ns,lu.to,lu.from,g4m_05_id) %>% summarise(value=sum(value)) %>%
  select(-c(g4m_05_id)) %>% ungroup()

downscalr_out <- res4
# END: Process RstLnd in res---

# Get output of merged Downscale & G4M spatial layer at ns resolution ------
results <- g4mid_to_simuid(downscalr_out,lab,project,scen1,scen2,scen3)

#results <- downscalr_out %>% map_df(~g4mid_to_simuid(.,lab,project,scen1,scen2,scen3)) %>% rbind
# ------------------------------------------------------------------------------

mapping_LC_BiodLink <- readRDS(file="mapping_LC_BiodLink.RData")[[1]]
results0 <- results
results <- results %>%
  left_join(mapping_LC_BiodLink,by=c("lu.from"="lu.linkoutput")) %>%
  select(-c(lu.from)) %>%
  rename(lu.from=lu.new) %>%
  left_join(mapping_LC_BiodLink,by=c("lu.to"="lu.linkoutput")) %>%
  select(-c(lu.to)) %>%
  rename(lu.to=lu.new) %>%
  select(REGION,times,ns,lu.to,value,lu.from)%>%
  group_by(REGION,times,ns,lu.to,lu.from) %>%
  summarise(value=sum(value)) %>%
  ungroup()

## YW: separate restored_nf age class-----
## The current method is to merge it directly into "OthNatLnd", so the RstLnd_nf will be identified as abandoned land as well, in the step "Get abandoned land"; further these restored_nf land will be mapped to "AbnCrpLnd_restored_nf" (mostly, because of f_nf_share for cells with RstLnd_nf are mostly 1), and this is desired.
results$lu.to[results$lu.to=="RstLnd_NoAffor"]="OthNatLnd"
results$lu.from[results$lu.from=="RstLnd_NoAffor"]="OthNatLnd"

results <- results %>%
  group_by(REGION,times,ns,lu.to,lu.from) %>%
  summarise(value=sum(value)) %>%
  ungroup()

## END YW separate restored_nf age class


# ------------------------------------------------------------------------------#
# start biodiversity computation ------
print("start biodiversity computation")

# Get potential npp
pot_npp <- readRDS(path("Input","pot_npp.RData")) %>%
  rename(ns=SimUID) %>% mutate(ns=as.factor(ns)) %>% dplyr::select(ns,pot_npp)


# Get built up areas
built <- readRDS(path("Input","built_area.RData")) %>% rename(ns=SimUID) %>%
  mutate(ns=as.factor(ns))


# Get utilization
g4m_utilization <- read.csv(str_glue("area_harvest_map_{project}_{lab}_{scen1}_{scen3}_{scen2}.csv")) %>%
  dplyr::select(g4m_id,used,year) %>% mutate(g4m_id=as.factor(g4m_id)) %>% rename(times=year)

g4m_mapping <- readRDS("g4m_mapping.RData")[[1]] %>% rename(g4m_id=g4m_05_id) %>% expand_grid(seq(2000,2100,10)) %>%
  rename(times=`seq(2000, 2100, 10)`)

managed_forests <- g4m_mapping %>% left_join(g4m_utilization) %>% drop_na() %>% dplyr::select(-g4m_id) %>%
  rename(ns=SimUID)

# Values for 2020
initial <- results %>% filter(times==2010)


#' Get abandoned land - !!! for now all classes are bundled into abandoned land------
#' since in GLOBIOM Trunk they are not differentiated. The BII and CFs are then
#' given by the average values of abandoned cropland and grassland

natlnd <- results %>% filter(lu.to == "OthNatLnd") %>%
  tidyr::complete(REGION,times,ns,lu.from,lu.to) %>% replace_na(list(value=0)) %>%
  group_by(REGION,times,ns) %>% summarise(value=sum(value))

natlnd_converted <- results %>% filter(lu.from == "RES", lu.to != "RES") %>%
  tidyr::complete(times,ns,lu.from,lu.to) %>% replace_na(list(value=0)) %>%
  group_by(times,ns) %>% summarise(value=sum(value)) %>%
  group_by(times,ns) %>% summarise(value=cumsum(value))

# Abandoned cropland
abn_crplnd <- results %>% filter(lu.from == "CrpLnd", lu.to == "OthNatLnd") %>%
  tidyr::complete(REGION,times,ns,lu.from,lu.to) %>% replace_na(list(value=0)) %>%
  group_by(REGION,ns,lu.from) %>% summarise(times=times,value=cumsum(value)) %>%
  mutate(lu.to="AbnCrpLnd") %>% relocate(times, .before=ns)

# Abandoned grassland
abn_grslnd <- results %>% filter(lu.from == "Grass", lu.to == "OthNatLnd") %>%
  tidyr::complete(REGION,times,ns,lu.from,lu.to) %>% replace_na(list(value=0)) %>%
  group_by(REGION,ns,lu.from) %>% summarise(times=times,value=cumsum(value)) %>%
  mutate(lu.to="AbnGrsLnd") %>% relocate(times, .before=ns)

# Abandoned forest plantations
abn_pltfor <- results %>% filter(lu.from == "PltFor", lu.to == "OthNatLnd")  %>%
  tidyr::complete(REGION,times,ns,lu.from,lu.to) %>% replace_na(list(value=0)) %>%
  group_by(REGION,ns,lu.from) %>% summarise(times=times,value=cumsum(value)) %>%
  mutate(lu.to="AbnPltFor") %>% relocate(times, .before=ns)

# Abandoned managed forests
abn_mngfor <- results %>% filter(lu.from %in% c("forest_old_ha","forest_new_ha"), lu.to == "OthNatLnd")  %>%
  tidyr::complete(REGION,times,ns,lu.from,lu.to) %>% replace_na(list(value=0)) %>%
  group_by(REGION,ns,lu.from) %>% summarise(times=times,value=cumsum(value)) %>%
  mutate(lu.from="forest",lu.to="AbnMngFor") %>% relocate(times, .before=ns) %>%
  group_by(REGION,times,ns,lu.from,lu.to) %>% summarise(value=sum(value))

# Merge land covers and correct OthNatLnd
natlnd_aux <- natlnd %>% rename(OthNatLnd=value) %>%
  left_join(abn_crplnd %>% rename(AbnCrpLnd=value) %>% ungroup %>% dplyr::select(-c(lu.to,lu.from))) %>%
  left_join(abn_grslnd %>% rename(AbnGrsLnd=value) %>% ungroup %>% dplyr::select(-c(lu.to,lu.from))) %>%
  left_join(abn_pltfor %>% rename(AbnPltFor=value) %>% ungroup %>% dplyr::select(-c(lu.to,lu.from))) %>%
  left_join(abn_mngfor %>% rename(AbnMngFor=value) %>% ungroup %>% dplyr::select(-c(lu.to,lu.from))) %>%
  left_join(natlnd_converted %>% rename(converted=value) %>% ungroup) %>%
  replace_na(list(AbnCrpLnd=0,AbnGrsLnd=0,AbnPltFor=0,AbnMngFor=0,converted=0)) #%>%
  #mutate(balance=OthNatLnd-(AbnCrpLnd+AbnGrsLnd+AbnPltFor+AbnMngFor))

# Rebalance cells with more abandoned land than natlnd
ini_col <- 4
ncols <- 9

natlnd_aux[,ini_col:ncols] <- apply(natlnd_aux[,ini_col:ncols],1,function(x) rebalance(x)) %>% t()

# Combined natlnd and abandoned areas
natlnd_cons <- natlnd_aux %>% dplyr::select(-converted) %>%
  mutate(OthNatLnd=OthNatLnd-(AbnCrpLnd+AbnGrsLnd+AbnPltFor+AbnMngFor)) %>%
  gather(lu.to,value,-c(REGION,times,ns)) #%>% replace(. < 0, 0)

# Compute age classes of abandoned areas
r_crplnd <- get_age(natlnd_cons %>% filter(lu.to=="AbnCrpLnd"),"AbnCrpLnd")
r_grslnd <- get_age(natlnd_cons %>% filter(lu.to=="AbnGrsLnd"),"AbnGrsLnd")
r_pltfor <- get_age(natlnd_cons %>% filter(lu.to=="AbnPltFor"),"AbnPltFor")
r_mngfor <- get_age(natlnd_cons %>% filter(lu.to=="AbnMngFor"),"AbnMngFor")

# Merge results
natlnd_cons <- natlnd_cons %>% filter(lu.to=="OthNatLnd") %>%
  rbind(r_crplnd,r_grslnd,r_pltfor,r_mngfor) %>% drop_na()


#' Get forest restored land  - for now afforestation areas that are not managed along ----
#' the simulation period are considered restoration.

# Get forest and non-forest shares
f_nf_share <- readRDS(path("input","forest_share.RData")) %>% rename(ns=SimUID) %>% mutate(ns=as.factor(ns))
f_nf_share$forest_share[which(f_nf_share$forest_share > 1)] <- 1
f_nf_share$nonforest_share[which(f_nf_share$nonforest_share < 0)] <- 0

# Get afforestation areas where the potential natural vegetation is forest and
# areas are not managed along the simulation period

restored <- results %>% filter(lu.to == "forest_new_ha") %>% group_by(REGION,times,ns) %>%
  summarise(value=sum(value)) %>% left_join(f_nf_share) %>% left_join(managed_forests) %>%
  mutate(value=value*forest_share) %>% filter(value > 0, used==0) %>%
  dplyr::select(-c(forest_share,nonforest_share,used)) %>% arrange(ns,times) %>%
  spread(times,value) %>% drop_na()

# if (dim(restored)[1] > 0) restored <- restored %>% mutate(`2000`=0) %>% relocate(`2000`,.before = `2010`) #sometimes 2010 transition may not exist
if (dim(restored)[1] > 0) restored <- restored %>% mutate(`2000`=0) %>% relocate(`2000`,.after = `ns`)


# Define simulation length
ncols <- dim(restored)[2]
ini_col <- 3

if (dim(restored)[1] > 0){
  # Get restoration from 0 to 10 years
  restored_10 <- restored
  restored_10[,ini_col:ncols] <- apply(restored[,ini_col:ncols],1,function(x) calc_lag(x,1)) %>% t()

  # Get remaining restoration classes
  restored_aux <- restored
  restored_aux[,ini_col:ncols] <- 0
  restored_20 <- restored_30 <- restored_40 <- restored_50 <-
  restored_60 <- restored_70 <- restored_80 <- restored_aux

  # Move afforestation along the time periods to the next age class
  if ((ini_col+1) <= ncols) restored_20[,(ini_col+1):ncols] <- restored_10[,ini_col:(ncols-1)]
  if ((ini_col+2) <= ncols) restored_30[,(ini_col+2):ncols] <- restored_10[,ini_col:(ncols-2)]
  if ((ini_col+3) <= ncols) restored_40[,(ini_col+3):ncols] <- restored_10[,ini_col:(ncols-3)]
  if ((ini_col+4) <= ncols) restored_50[,(ini_col+4):ncols] <- restored_10[,ini_col:(ncols-4)]
  if ((ini_col+5) <= ncols) restored_60[,(ini_col+5):ncols] <- restored_10[,ini_col:(ncols-5)]
  if ((ini_col+6) <= ncols) restored_70[,(ini_col+6):ncols] <- restored_10[,ini_col:(ncols-6)]
  if ((ini_col+7) <= ncols) restored_80[,ini_col+7] <- restored_10[,ini_col]
  if ((ini_col+8) <= ncols) restored_80[,ini_col+8] <- restored_10[,ini_col+1] + restored_10[,ini_col]
  if ((ini_col+9) <= ncols) restored_80[,ini_col+9] <- restored_10[,ini_col+2] + restored_10[,ini_col+1] + restored_10[,ini_col]
  if ((ini_col+10) <= ncols) restored_80[,ini_col+10] <- restored_10[,ini_col+3] + restored_10[,ini_col+2] + restored_10[,ini_col+1] + restored_10[,ini_col]

  # Combined the different age classes
  restoration <- restored_10 %>% gather(times,value,-c(REGION,ns)) %>% rename(restored_10=value) %>%
    left_join(restored_20 %>% gather(times,value,-c(REGION,ns)) %>% rename(restored_20=value)) %>%
    left_join(restored_30 %>% gather(times,value,-c(REGION,ns)) %>% rename(restored_30=value)) %>%
    left_join(restored_40 %>% gather(times,value,-c(REGION,ns)) %>% rename(restored_40=value)) %>%
    left_join(restored_50 %>% gather(times,value,-c(REGION,ns)) %>% rename(restored_50=value)) %>%
    left_join(restored_60 %>% gather(times,value,-c(REGION,ns)) %>% rename(restored_60=value)) %>%
    left_join(restored_70 %>% gather(times,value,-c(REGION,ns)) %>% rename(restored_70=value)) %>%
    left_join(restored_80 %>% gather(times,value,-c(REGION,ns)) %>% rename(restored_80=value)) %>%
    mutate(times=as.integer(times))

  # Afforestation area
  afflnd <- results %>% filter(lu.to == "forest_new_ha") %>%
    tidyr::complete(REGION,times,ns,lu.from,lu.to) %>% replace_na(list(value=0)) %>%
    group_by(REGION,times,ns) %>% summarise(value=sum(value))

  afflnd_aux <- afflnd %>% rename(forest_new_ha=value) %>%
    left_join(restoration) %>% replace(is.na(.), 0) %>%
    mutate(forest_new_ha=forest_new_ha-rowSums(across(starts_with("restored")), na.rm = T))

  # Aggregate and clean
  afflnd_cons <- afflnd_aux[,1:min(ncols,dim(afflnd_aux)[2])] %>%
    gather(lu.to,value,-c( REGION,times,ns)) %>% replace(. < 0, 0)
} else {
  afflnd_cons <- results %>% slice(1) %>% dplyr::select(-lu.from) %>%
    mutate(lu.to="forest_new_ha",value=0)
}



# Get BII indicator ------
if (get_bii) {
  print("bii calculation")

  # Get BII coeficients and add restored abandoned land
  bii_coefs <- readRDS(path("Input","BII_coefs.RData"))


  #' Compute time dynamics of abandoned areas. Here the same assumption
  #' employed in Leclere et al. (2020) for the cSAR17 is made: restored
  #' areas are assume to linearly increase in quality  with the final
  #' state 80% as good as the undisturbed vegetation when the age is > 70 years.
  #' For land covers where the disturbed ecosystem is better than the
  #' undisturbed * 0.8 (MngFor and PltFor), it
  #' is assumed that they converge to undisturbed ecosystem values at
  #' the end of the recovery period

  age <- seq(10,80,10)
  ecosystem <- c("f","nf")

  # Land cover classes
  lorg <- c("cropland_other","grassland","cropland_2Gbioen","forest_managed")
  lurest <- c("AbnCrpLnd","AbnGrsLnd","AbnPltFor","AbnMngFor")
  lu <- tibble(LUclass=lorg,LC_TYPE=lurest)

  # Mapping between land and BII land
  land_map <- tibble(LUclass=lorg,LC_TYPE=lurest)

  # Undisturbed BII coefficients
  org_map <- tibble(Eco=ecosystem,BII_ORG=c(bii_coefs$BII[which(bii_coefs$LUclass=="forest_unmanaged")],
                                             bii_coefs$BII[which(bii_coefs$LUclass=="other_nf")]))

  # Compute coefficients according to age of abandonment
  restored_coef <- expand.grid(lurest,age,ecosystem) %>% rename(LC_TYPE = Var1, Age = Var2, Eco = Var3) %>%
    mutate(class=str_c(LC_TYPE,"_restored_",Age,"_",Eco)) %>%
    filter(! (LC_TYPE=="AbnMngFor" & Eco=="nf")) %>%
    left_join(land_map) %>% mutate(LUclass=str_c(LUclass,"_",Eco)) %>%
    mutate(LUclass = recode(LUclass, "forest_managed_f" = "forest_managed")) %>%
    left_join(bii_coefs) %>% left_join(org_map) %>% mutate(BII_ORG=ifelse(BII<BII_ORG*.8,BII_ORG*.8,BII_ORG)) %>%
    mutate(BII_FINAL= ifelse( Age >= 20 & Age <= 70 , BII + (Age - 10) /(70 - 5) * (BII_ORG - BII),
                              ifelse( Age > 70, BII_ORG, BII))) %>%
    dplyr::select(class,BII_FINAL,weightSum) %>%
    rename(LUclass=class,BII=BII_FINAL)

  # Compute coefficients according to age of restoration
  aff_restored <- expand.grid("restored_f",age) %>% rename(LUclass = Var1, Age = Var2) %>%
    mutate(class=str_c(LUclass,"_",Age)) %>% left_join(bii_coefs) %>%
    mutate(BII_INI=bii_coefs$BII[which(bii_coefs$LUclass=="forest_managed")]) %>%
    mutate(BII_FINAL= ifelse( Age >= 20 & Age <= 70 , BII_INI + (Age - 10) /(70 - 5) * (BII - BII_INI),
                             ifelse( Age > 70, BII, BII_INI))) %>%
    dplyr::select(class,BII_FINAL,weightSum) %>%
    rename(LUclass=class,BII=BII_FINAL)

  ## Merge bii_coefs coeficients-----
  bii_coefs <- bii_coefs %>% rbind(restored_coef,aff_restored)

  ## Get mapping between GLOBIOM land cover types and BII land cover types
  lc_map <- readRDS(path("Input","BTC_LC_MAP.RData"))

  ## Define lc_map mapping -----
  lc_type <- c("Grass","OthNatLnd","CrpLnd","forest_managed","forest_unmanaged","PltFor","built_area",
               "AbnCrpLnd","AbnGrsLnd","AbnPltFor","AbnMngFor",
               "restored_10","restored_20","restored_30","restored_40","restored_50","restored_60","restored_70","restored_80",
               "AbnCrpLnd_restored_10","AbnCrpLnd_restored_20","AbnCrpLnd_restored_30","AbnCrpLnd_restored_40","AbnCrpLnd_restored_50","AbnCrpLnd_restored_60","AbnCrpLnd_restored_70","AbnCrpLnd_restored_80",
               "AbnGrsLnd_restored_10","AbnGrsLnd_restored_20","AbnGrsLnd_restored_30","AbnGrsLnd_restored_40","AbnGrsLnd_restored_50","AbnGrsLnd_restored_60","AbnGrsLnd_restored_70","AbnGrsLnd_restored_80",
               "AbnPltFor_restored_10","AbnPltFor_restored_20","AbnPltFor_restored_30","AbnPltFor_restored_40","AbnPltFor_restored_50","AbnPltFor_restored_60","AbnPltFor_restored_70","AbnPltFor_restored_80",
               "AbnMngFor_restored_10","AbnMngFor_restored_20","AbnMngFor_restored_30","AbnMngFor_restored_40","AbnMngFor_restored_50","AbnMngFor_restored_60","AbnMngFor_restored_70","AbnMngFor_restored_80","RstLnd","protected_priforest","protected_other","RstLnd_NoAffor")
  # add the mapping for RstLnd_NoAffor

  bii_type <- c("grassland","other","cropland_other","forest_managed","forest_unmanaged","cropland_2Gbioen","built.up",
                "abn_cropland_other","abn_grassland","abn_cropland_2Gbioen","abn_forest_managed",
                "restored_f_10","restored_f_20","restored_f_30","restored_f_40","restored_f_50","restored_f_60","restored_f_70","restored_f_80",
                "AbnCrpLnd_restored_10","AbnCrpLnd_restored_20","AbnCrpLnd_restored_30","AbnCrpLnd_restored_40","AbnCrpLnd_restored_50","AbnCrpLnd_restored_60","AbnCrpLnd_restored_70","AbnCrpLnd_restored_80",
                "AbnGrsLnd_restored_10","AbnGrsLnd_restored_20","AbnGrsLnd_restored_30","AbnGrsLnd_restored_40","AbnGrsLnd_restored_50","AbnGrsLnd_restored_60","AbnGrsLnd_restored_70","AbnGrsLnd_restored_80",
                "AbnPltFor_restored_10","AbnPltFor_restored_20","AbnPltFor_restored_30","AbnPltFor_restored_40","AbnPltFor_restored_50","AbnPltFor_restored_60","AbnPltFor_restored_70","AbnPltFor_restored_80",
                "AbnMngFor_restored_10_f","AbnMngFor_restored_20_f","AbnMngFor_restored_30_f","AbnMngFor_restored_40_f","AbnMngFor_restored_50_f","AbnMngFor_restored_60_f","AbnMngFor_restored_70_f","AbnMngFor_restored_80_f","restored_nf","forest_unmanaged","other","restored_nf")

  lc_map <- tibble(LC_TYPE = lc_type,LC_BTC=bii_type)

  #results_bii <- results %>% filter(! lu.to %in% c("OthNatLnd","forest_new_ha")) %>%
  results_bii <- results %>% filter(! lu.to %in% c("OthNatLnd","forest_new_ha","RstLnd_NoAffor")) %>%
    group_by(REGION,times,ns,lu.to) %>% summarise(value=sum(value)) %>%
    rbind(natlnd_cons,afflnd_cons) %>%
    left_join(f_nf_share) %>%
    left_join(managed_forests) %>% left_join(built) %>% spread(lu.to,value) %>%
    mutate(forest_managed=used*(forest_new_ha+forest_old_ha),
           forest_unmanaged=(1-used)*(forest_new_ha+forest_old_ha)) %>%
    dplyr::select(-c(used,forest_new_ha,forest_old_ha)) %>%
    gather(LC_TYPE,value,-c(REGION,times,ns,forest_share,nonforest_share))


  initial_bii <- initial %>% group_by(REGION,times,ns,lu.from) %>% summarise(value=sum(value)) %>%
    left_join(f_nf_share) %>%
    left_join(managed_forests) %>% left_join(built) %>% spread(lu.from,value) %>%
    mutate(forest_managed=used*(forest_old_ha),forest_unmanaged=(1-used)*(forest_old_ha)) %>%
    dplyr::select(-c(used,forest_old_ha)) %>%
    gather(LC_TYPE,value,-c(REGION,times,ns,forest_share,nonforest_share)) %>%
    mutate(times=2000)


  results_bii <- results_bii %>% rbind(initial_bii) %>%
    left_join(pot_npp) %>% replace(is.na(.),0)

  global_npp <- results_bii %>% filter(times==2010) %>% mutate(pot_npp=pot_npp*value) %>%
    group_by() %>% summarise(pot_npp=sum(pot_npp)) %>% pull()

  # Merge BII coefficients - for now we split the tibble to speed up processing
  results_bii <- results_bii %>% group_by(times) %>% group_split()

  # Merge results
  results_bii <- results_bii %>% map_df(~process_bii(.)) %>% rbind

  # Calculate BII at the Simu and aggregate level
  bii_simu <- results_bii %>% mutate(value=value*f_nf_share) %>% mutate(pot_npp=10^6*pot_npp*value/global_npp) %>%
    mutate(score=value*BII,score_prod=pot_npp*BII) %>% group_by(REGION,times,ns) %>%
    summarise(score=sum(score,na.rm = TRUE)/sum(value,na.rm = TRUE),score_prod=sum(score_prod,na.rm = TRUE)/sum(pot_npp,na.rm = TRUE),value=sum(value,na.rm = TRUE)) %>% mutate(ScenNr=scen)

  bii <- results_bii %>% mutate(value=value*f_nf_share) %>% mutate(pot_npp=10^6*pot_npp*value/global_npp) %>%
    mutate(score=value*BII,score_prod=pot_npp*BII) %>% group_by(times) %>%
    summarise(score=sum(score,na.rm = TRUE)/sum(value,na.rm = TRUE),value=sum(value,na.rm = TRUE),score_prod=sum(score_prod,na.rm = TRUE)/sum(pot_npp,na.rm = TRUE)) %>%
    mutate(ScenNr=scen)

} else {
  bii_simu <- bii <- NULL
}

# Get cSAR indicator------
if (get_csar){
  print("csar calculation")

  ## Get cSAR coeficients-----
  csar_coefs <- readRDS(path("Input","PSLglobal_extended.RData")) %>%
    mutate(forest_managed=0.5*(extensive_forestry+intensive_forestry)) %>%#,
    gather(LUclass,CF,-c(ecoregion))


  #' Compute time dynamics of abandoned areas. Here the same assumption
  #' employed in Leclere et al. (2020) for the cSAR17 is made: restored
  #' areas are assume to linearly increase in quality  with the final
  #' state 80% as good as the undisturbed vegetation when the age is > 70 years.
  #' For land covers where the disturbed ecosystem is better than the
  #' undisturbed * 0.8 (MngFor and PltFor in non-forest ecosystem), it
  #' was assumed that they converge to undsturbed ecosystem values at
  #' the end of the recovery period

  age <- seq(10,80,10)
  ecoreg <- csar_coefs$ecoregion %>% unique()

  ## Land cover classes
  lorg <- c("annual_crops","pasture","permanent_crops","forest_managed")
  lurest <- c("AbnCrpLnd","AbnGrsLnd","AbnPltFor","AbnMngFor")
  lu <- tibble(LUclass=lorg,LC_TYPE=lurest)

  ## Mapping between land and BII land-----
  land_map <- tibble(LUclass=lorg,LC_TYPE=lurest)

  ## Compute coefficients for abandoned areas-----
  restored_coef <- expand.grid(lurest,age,ecoreg) %>% rename(LC_TYPE = Var1, Age = Var2, ecoregion = Var3) %>%
    mutate(class=str_c(LC_TYPE,"_restored_",Age)) %>%
    left_join(land_map) %>% left_join(csar_coefs) %>%
    mutate(CF_ORG=CF*0.2) %>% mutate(CF_FINAL=
                                       ifelse( Age >= 20 & Age <= 70 , CF - ((Age - 10)/(70-5))*(CF - 0.2*CF),
                                               ifelse(Age > 70 , CF_ORG, CF))) %>%
    dplyr::select(ecoregion,class,CF_FINAL) %>% rename(LUclass=class,CF=CF_FINAL)

  csar_coefs <- csar_coefs %>% rbind(restored_coef)

  # Get simu to ecoregion map
  eco_map <- readRDS(path("Input","ecoregions_share.RData")) %>% rename(ns=SimUID) %>%
    mutate(ns=as.factor(ns))

  # Excluded ecoregions
  my_excluded_ers <- c('AA0101','AA0102','AA0103','AA0108','AA0109','AA0110','AA0114','AA0118','AA0119','AA0125','AA0126','AA0202','AA0401','AA0407','AA1101','AA1401',
                       'AN1101','AN1102','AN1103','AN1104',
                       'AT0105','AT0113','AT0120','AT0127','AT0201','AT0703','AT0720','AT0802','AT0803','AT0902','AT1301','AT1308','AT1318','AT1405',
                       'IM0101','IM0110','IM0114','IM0125','IM0127','IM0133','IM0143','IM0148','IM0156','IM0170','IM1401','IM1403',
                       'NA0301','NA0525','NA0615','NA1102','NA1109','NA1112','NA1113','NA1117',
                       'NT0102','NT0110','NT0116','NT0123','NT0134','NT0172','NT0179','NT0216','NT0218','NT0220','NT0226','NT0301','NT0401','NT0403','NT0705','NT1301','NT1305','NT1306','NT1307','NT1311','NT1314','NT1318','NT1402','NT1404','NT1406',
                       'OC0101','OC0102','OC0103','OC0104','OC0105','OC0106','OC0107','OC0108','OC0109','OC0110','OC0111','OC0112','OC0113','OC0114','OC0115','OC0116','OC0117','OC0201','OC0202','OC0203','OC0204','OC0701','OC0702','OC0703',
                       'PA0403','PA0425','PA0438','PA0807','PA1101','PA1109','PA1113','PA1203'
  )

  # Get mapping between GLOBIOM land cover types and cSAR land cover types
  lc_type <- c("Grass","OthNatLnd","CrpLnd","forest_managed","forest_unmanaged","PltFor","built_area",
               "restored_10","restored_20","restored_30","restored_40","restored_50","restored_60","restored_70","restored_80",
               "AbnCrpLnd_restored_10","AbnCrpLnd_restored_20","AbnCrpLnd_restored_30","AbnCrpLnd_restored_40","AbnCrpLnd_restored_50","AbnCrpLnd_restored_60","AbnCrpLnd_restored_70","AbnCrpLnd_restored_80",
               "AbnGrsLnd_restored_10","AbnGrsLnd_restored_20","AbnGrsLnd_restored_30","AbnGrsLnd_restored_40","AbnGrsLnd_restored_50","AbnGrsLnd_restored_60","AbnGrsLnd_restored_70","AbnGrsLnd_restored_80",
               "AbnPltFor_restored_10","AbnPltFor_restored_20","AbnPltFor_restored_30","AbnPltFor_restored_40","AbnPltFor_restored_50","AbnPltFor_restored_60","AbnPltFor_restored_70","AbnPltFor_restored_80",
               "AbnMngFor_restored_10","AbnMngFor_restored_20","AbnMngFor_restored_30","AbnMngFor_restored_40","AbnMngFor_restored_50","AbnMngFor_restored_60","AbnMngFor_restored_70","AbnMngFor_restored_80","RstLnd_NoAffor")

  btc_type <- c("grassland","other","cropland_other","forest_managed","forest_unmanaged","cropland_2Gbioen","built_area",
                "restored_10","restored_20","restored_30","restored_40","restored_50","restored_60","restored_70","restored_80",
                "AbnCrpLnd_restored_10","AbnCrpLnd_restored_20","AbnCrpLnd_restored_30","AbnCrpLnd_restored_40","AbnCrpLnd_restored_50","AbnCrpLnd_restored_60","AbnCrpLnd_restored_70","AbnCrpLnd_restored_80",
                "AbnGrsLnd_restored_10","AbnGrsLnd_restored_20","AbnGrsLnd_restored_30","AbnGrsLnd_restored_40","AbnGrsLnd_restored_50","AbnGrsLnd_restored_60","AbnGrsLnd_restored_70","AbnGrsLnd_restored_80",
                "AbnPltFor_restored_10","AbnPltFor_restored_20","AbnPltFor_restored_30","AbnPltFor_restored_40","AbnPltFor_restored_50","AbnPltFor_restored_60","AbnPltFor_restored_70","AbnPltFor_restored_80",
                "AbnMngFor_restored_10","AbnMngFor_restored_20","AbnMngFor_restored_30","AbnMngFor_restored_40","AbnMngFor_restored_50","AbnMngFor_restored_60","AbnMngFor_restored_70","AbnMngFor_restored_80","restored_nf")

  # add the mapping for RstLnd_NoAffor

  csar_type <- c("pasture",NA,"annual_crops","forest_managed","primary","permanent_crops","urban",
                 "restored0to10Years","restored10to20Years","restored20to30Years","restored30to40Years",
                 "restored40to50Years","restored50to60Years","restored60to70Years","restored70YearsOrOlder",
                 "AbnCrpLnd_restored_10","AbnCrpLnd_restored_20","AbnCrpLnd_restored_30","AbnCrpLnd_restored_40","AbnCrpLnd_restored_50","AbnCrpLnd_restored_60","AbnCrpLnd_restored_70","AbnCrpLnd_restored_80",
                 "AbnGrsLnd_restored_10","AbnGrsLnd_restored_20","AbnGrsLnd_restored_30","AbnGrsLnd_restored_40","AbnGrsLnd_restored_50","AbnGrsLnd_restored_60","AbnGrsLnd_restored_70","AbnGrsLnd_restored_80",
                 "AbnPltFor_restored_10","AbnPltFor_restored_20","AbnPltFor_restored_30","AbnPltFor_restored_40","AbnPltFor_restored_50","AbnPltFor_restored_60","AbnPltFor_restored_70","AbnPltFor_restored_80",
                 "AbnMngFor_restored_10","AbnMngFor_restored_20","AbnMngFor_restored_30","AbnMngFor_restored_40","AbnMngFor_restored_50","AbnMngFor_restored_60","AbnMngFor_restored_70","AbnMngFor_restored_80","AbnCrpLnd_restored_30")
  #YW: currently cannot find a good mapper for RstLnd_NoAffor; using "AbnCrpLnd_restored_30" or maybe "restored0to10Years"


  lc_map <- tibble(LC_TYPE = lc_type,LC_BTC=btc_type)
  csar_lc_map <- tibble(LC_BTC=btc_type,LUclass=csar_type)

  #results_csar <- results %>% filter(! lu.to %in% c("OthNatLnd","forest_new_ha")) %>%
 results_csar <- results %>% filter(! lu.to %in% c("OthNatLnd","forest_new_ha","RstLnd_NoAffor")) %>%
    group_by(REGION,times,ns,lu.to) %>% summarise(value=sum(value)) %>%
    rbind(natlnd_cons,afflnd_cons) %>%
    left_join(managed_forests) %>% left_join(built) %>% spread(lu.to,value) %>%
    mutate(forest_managed=used*(forest_new_ha+forest_old_ha),forest_unmanaged=(1-used)*(forest_new_ha+forest_old_ha)) %>%
    dplyr::select(-used) %>%
    gather(LC_TYPE,value,-c(REGION,times,ns))

  initial_csar <- initial %>% group_by(REGION,times,ns,lu.from) %>% summarise(value=sum(value)) %>%
    left_join(managed_forests) %>% left_join(built) %>% spread(lu.from,value) %>%
    mutate(forest_managed=used*(forest_old_ha),forest_unmanaged=(1-used)*(forest_old_ha)) %>%
    dplyr::select(-used) %>%
    gather(LC_TYPE,value,-c(REGION,times,ns)) %>%
    mutate(times=2000)

  results_csar <- results_csar %>% rbind(initial_csar)

  # Merge coefficients
  results_csar <- results_csar %>% group_by(times) %>% group_split()

  # Add coefficients and compute impact
  results_csar <- results_csar %>% map_df(~process_csar(.)) %>% rbind

  # Calculate extinctions at ecoregion level
  csar_simu <- results_csar %>% filter(! ecoregion %in% my_excluded_ers) %>%
    group_by(REGION,times,ns) %>%
    summarise(extinctions=sum(extinctions)) %>%
    mutate(score=1-extinctions) %>% mutate(score=ifelse(score <0,0,score)) %>% mutate(ScenNr=scen)

  # Calculate total extinctions
  csar <- results_csar %>% filter(! ecoregion %in% my_excluded_ers) %>%
    group_by(REGION,times,ecoregion) %>%
    summarise(extinctions=sum(extinctions)) %>%
    mutate(score=extinctions) %>% mutate(score=ifelse(score <0,0,score)) %>%
    group_by(REGION,times) %>% summarise(score=1-sum(score)) %>% mutate(ScenNr=scen)

} else {
  csar_simu <- csar <- NULL
}

# Save ns and REGION level biodiv results ------
#saveRDS(list(bii_simu,bii,csar_simu,csar,downscalr_out0,downscalr_out,results0,results),paste0("Output/biodiversity_",project,"_",lab,".RData")) ## for diagnosis: unload also intermediate LUC results
saveRDS(list(bii_simu,bii,csar_simu,csar,results),paste0("Output/biodiversity_",project,"_",lab,".RData")) ## for normal pipeline runs: only unload important results

