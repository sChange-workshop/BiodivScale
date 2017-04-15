##======================================================================
##	BBS case study for scale paper...
##======================================================================
rm(list=ls())
##	
library(dggridR)
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(scales)
library(vegan)
library(iNEXT)
library(nlme)
library(purrr)

##	local path (add your own as appropriate...) to the dropbox data
path <- '~/Dropbox/BioTIMELatest/'
##	get the data and metadata
##	Get the raw data locally 

bt <- read_csv(paste0(path, 'BioTIMEMarch.csv'))
#bt <- read_csv('/home/sarah/Dropbox/BioTIMELatest/bioTIMEFeb.csv')
#bt <- read_csv('/Users/HotStuff/Dropbox/BioTIMELatest/bioTIMEFeb.csv')

##	Get the meta data locally
meta <- read_csv(paste0(path, 'bioTIMEmetadataFeb.csv'))
#meta <- read.csv('/home/sarah/Dropbox/BioTIMELatest/bioTIMEmetadataFeb.csv')
#meta <- read.csv('/Users/HotStuff/Dropbox/BioTIMELatest/bioTIMEmetadataFeb.csv')

##	join abundance records with the metadata
bt <- inner_join(meta, bt, by='STUDY_ID')

# Add columns for Species species name 
bt <- bt %>% unite(col=Species, GENUS, SPECIES, remove=FALSE) %>%
	select(-GENUS_SPECIES)

# get the BBS data
BBS <- bt %>% filter(STUDY_ID==195)

##======================================================================
##	want to put these data into cells of increasing size 
##	NB: nested hexagonal cells overlap?! Is this a problem for inferences regarding
##	the temporal trends at different scales...

##	res=16 ~ 1.2 km2
dgg16 <- dgconstruct(res=16)
##	res=14 ~ 10.7 km2
dgg14 <- dgconstruct(res=14)
##	res=12 ~ 96 km2
dgg12 <- dgconstruct(res=12)
##	res=10 ~ 863.8 km2
dgg10 <- dgconstruct(res=10)
##	res=8 ~ 7774.2 km2
dgg8 <- dgconstruct(res=8)
##	res=6 ~ 70 000 km2
dgg6 <- dgconstruct(res=6)
##	res=4 = 629,710.64410 km2
dgg4 <- dgconstruct(res=4)
##	res=2 = 5,667,395.79693 km2
dgg2 <- dgconstruct(res=2)
##	res=0 = 51,006,562.17241 km2
dgg0 <- dgconstruct(res=0)


##	get the corresponding grid cells for all observations (and for each scale)
BBS_0 <- BBS %>% mutate(
	resolution = 0,
	cell = dgtransform(dgg0, lat=LATITUDE, lon=LONGITUDE))

BBS_2 <- BBS %>% mutate(
	resolution = 2,
	cell = dgtransform(dgg2, lat=LATITUDE, lon=LONGITUDE))

BBS_4 <- BBS %>% mutate(
	resolution = 4,
	cell = dgtransform(dgg4, lat=LATITUDE, lon=LONGITUDE))

##	res=6
BBS_6 <- BBS %>% mutate(
	resolution = 6,
	cell = dgtransform(dgg6, lat=LATITUDE, lon=LONGITUDE))

##	res=8
BBS_8 <- BBS %>% mutate(
	resolution = 8,
	cell = dgtransform(dgg8, lat=LATITUDE, lon=LONGITUDE))

##	res=10
BBS_10 <- BBS %>% mutate(
	resolution = 10,
	cell = dgtransform(dgg10, lat=LATITUDE, lon=LONGITUDE))

##	res=12
BBS_12 <- BBS %>% mutate(
	resolution = 12,
	cell = dgtransform(dgg12, lat=LATITUDE, lon=LONGITUDE))

##	res=14
BBS_14 <- BBS %>% mutate(
	resolution = 14,
	cell = dgtransform(dgg14, lat=LATITUDE, lon=LONGITUDE))

##	res=16
BBS_16 <- BBS %>% mutate(
	resolution = 16,
	cell = dgtransform(dgg16, lat=LATITUDE, lon=LONGITUDE))
	
##	put 'em back together
BBS_multi <- bind_rows(BBS_0, BBS_2, BBS_4, BBS_6, BBS_8, BBS_10, BBS_12, BBS_14, BBS_16)	

##	what are the cell counts?	
check_cells <- BBS_multi %>% group_by(resolution, YEAR) %>% summarise(N_cell = n_distinct(cell))	
##	looks like there is very little gained by going smaller than res=10...perhaps 12?
ggplot(check_cells) +
	facet_wrap(~resolution, nrow=3) + 
	geom_point(aes(x=YEAR, y= N_cell), alpha=0.5, colour='grey') +
	stat_smooth(method='lm', aes(x=YEAR, y= N_cell), lwd=0.1, alpha=0.25, se=F) +
#	scale_y_log10() +
	theme_bw()

##======================================================================
##	rename abundance column 
BBS_multi <- BBS_multi %>%
	mutate(Abundance = sum.allrawdata.ABUNDANCE) %>%
	select(-sum.allrawdata.ABUNDANCE, -sum.allrawdata.BIOMASS)

##	collate Species within the new cells 
BBS_multi <- BBS_multi %>%
	group_by(REALM, CLIMATE, TAXA, resolution, YEAR, cell, Species) %>%
	summarise(
		Abundance = sum(Abundance)) %>%
	ungroup()

##	time series of species richness within cells
checkS <- BBS_multi %>% group_by(resolution, YEAR, cell) %>% summarise(N_spp = n_distinct(Species))

ggplot(filter(checkS, resolution < 11 )) +
	facet_wrap(~resolution, nrow=3, scales='free') + 
	geom_point(aes(x=YEAR, y=N_spp), alpha=0.5, colour='grey') +
	stat_smooth(method='lm', aes(x=YEAR, y=N_spp, group=cell), lwd=0.1, alpha=0.25, se=F) +
#	scale_y_log10() +
	theme_bw()
	
##	calculate coverage of the new BBS time series
BBS_coverage <- BBS_multi %>% 
  group_by(resolution, YEAR, cell) %>%
  summarise(
    # how many singletons
    singletons = sum(Abundance==1),
    # how many doubletons
    doubletons = sum(Abundance==2),
    # how many individuals in total sample
    N = sum(Abundance),
    # eqn 4a in Chao & Jost 2012 Ecology (==eqn 12 in Chao et al 2014 Ecol Monogr)
    Chat = 1 - (singletons/N) * (((N-1)*singletons)/((N-1)*singletons + 2*doubletons)),
    # with fix from Chao et al 2014 Subroutine for iNEXT code (appendix)
    # correction for communities with no doubletons (this prevents NaN results for either singletons==0 or doubletons==0)
    Chat_corrected = ifelse(doubletons>0, 
                            1 - (singletons/N) * (((N-1)*singletons)/((N-1)*singletons + 2*doubletons)), 
                            1 - (singletons/N) * (((N-1)*singletons)/((N-1)*singletons + 2))),
    # the iNEXT coverage calculation has some extras in it that I don't understand
    # e.g., it corrects using f0.hat, the number of unobserved species (Chao1 estimator)? Why? Exactly how? Ref?
    coverage = DataInfo(Abundance)$SC) %>%
  ungroup()
  
 ggplot(BBS_coverage) +
 	facet_wrap(~resolution) +
 	geom_point(aes(x=YEAR, y=Chat_corrected), alpha=0.5, colour='grey') +
	stat_smooth(method='lm', aes(x=YEAR, y=Chat_corrected, group=cell), lwd=0.1, alpha=0.25, se=F) +
	theme_bw() 
	
BBS_coverage %>% filter(Chat_corrected < 0.85)	# minimum coverage is 0.82, proceed with all data...
##======================================================================
##	I think coverage- or individual-based rarefaction would be best here. 
##	But for consistency with other analyses Patrick has done, a similar sample-based
##	approach is OK.

##======================================================================
##	models comparable with those fit by Patrick for scale-analyses of BioTIME:
##	teomporal autocorrelation, though no grouping variables here; each scale modellled
##	independently

##	calculate metrics (this would be done in the rarefaction step to be included later...)
BBS_biodiversity <- BBS_multi %>%
	# calculate metrics for resolution/year/cell combinations
	group_by(resolution, YEAR, cell) %>%
	summarise(
		N = sum(Abundance),
		S = n_distinct(Species),
		hill1 = renyi(Abundance, scales=1, hill=TRUE)[1],
		PIE = diversity(Abundance, index='simpson'),
		ENSPIE = diversity(Abundance, index='invsimpson')) %>%
	ungroup()

##	centre year on 1992
BBS_biodiversity %>% distinct(YEAR) %>% summarise(min(YEAR), max(YEAR), median(YEAR))
BBS_biodiversity$cYEAR <- BBS_biodiversity$YEAR - 1992

##	first nested df
BBS_nest <- BBS_biodiversity %>%
	filter(resolution < 11) %>%
	group_by(resolution, cell) %>%                                          
	nest()

##	models to be fit to the data: I've modelled the richness related metrics here on
## 	the natural scale because I think that is what Patrick did??
model_list <- list(
	S_gls = function(x) gls(S ~ cYEAR, correlation = corAR1(form = ~cYEAR), data=x),
	logN_gls = function(x) gls(log(N) ~ cYEAR, correlation = corAR1(form = ~cYEAR), data=x),
	ENSPIE_gls = function(x) gls(ENSPIE ~ cYEAR, correlation = corAR1(form = ~cYEAR), data=x),
	hill1_gls = function(x) gls(hill1 ~ cYEAR, correlation = corAR1(form = ~cYEAR), data=x))


##	check with one time series
tmp_data <- BBS_nest %>% filter(resolution==0 & cell==2)
S_gls(unnest(tmp_data))	# ok

fn_model <- function(.model, df){
  # possibly replaces errors with a default value (here set to be NULL)
  df$model <- purrr::map(df$data, possibly(.model, otherwise=NULL))
  df
}

##	For each element of the list of model functions,
##	run the inner-loop function, and row-bind the results into a df 
gls_model_nested <- model_list %>%
	map_df(fn_model, BBS_nest, .id='model_id') %>%
	# add check for fits that failed (and returned NULL)
	mutate(is_null = map_lgl(model, is.null))

## how many failed?	Good!
sum(gls_model_nested$is_null)	

##	no tidier for gls! Quick hack...need to fix your own path to the github repository here
path2 <- '~/Desktop/current/BioTime/Scale_MS/BiodivScale/'
source(paste0(path2, 'tidy_gls.R'))

gls_model_params <- gls_model_nested %>%
	dplyr::select(model_id, resolution, cell, model) %>%
	# tidy function 
	mutate(coef = purrr::map(model, tidy.gls)) %>%
	dplyr::select(-model) %>%
	# and unpack
	unnest()

##	we want the slopes (and possibly the se's for variance-weighted meta-analytic models...)
gls_slopes <- gls_model_params %>%
	dplyr::select(model_id, resolution, cell, term, estimate, std.error, statistic, p.value) %>%
	spread(key=term, value=estimate) %>%
	dplyr::filter(!is.na(cYEAR)) %>%
	# rename YEAR as slope
	mutate(slope = cYEAR) %>%
	# drop the intercept term
	dplyr::select(model_id, resolution, cell, slope, std.error, statistic, p.value) 


##	add a better scale covariate for plotting
gls_slopes %>% distinct(resolution)
gls_slopes <- gls_slopes %>%
	mutate(
		# values for each resolution came from https://cran.r-project.org/web/packages/dggridR/vignettes/dggridR.html 
		scale = ifelse(resolution==0, round(51006562.17241),
			ifelse(resolution==2, round(5667395.79693),
			ifelse(resolution==4, round(629710.64410),
			ifelse(resolution==6, round(69967.84934),
			ifelse(resolution==8, round(7774.20548),
			ifelse(resolution==10, round(863.80061), NA)))))))

##	first look...
ggplot(gls_slopes, aes(x=scale, y=slope)) +
	# allow scales to vary as N was modelled on log-scale
	facet_wrap(~model_id, scales='free') +
	geom_point() +
	scale_x_log10() +
	geom_hline(yintercept=0, lty=2) +
	# gam to see if we are likely to recover hump-shaped prediction?
	stat_smooth(method='gam', se=F, formula = y ~ s(x, bs='cr', k=4)) +
	# or 2nd-order polynomial on quartiles (as per one of Patrick's other plots)... 
	stat_quantile(formula = y ~ poly(x, 2), lty=2, quantiles=c(0.25, 0.75)) +
	xlab('Scale (km2)') +
	theme_bw()
	

ggplot(gls_slopes, aes(x=scale, y=std.error)) +
	# allow scales to vary as N was modelled on log-scale
	facet_wrap(~model_id, scales='free') +
	geom_point() +
	scale_x_log10() +
	geom_hline(yintercept=0, lty=2) +
	# gam to see if we are likely to recover hump-shaped prediction?
	stat_smooth(method='gam', se=F, formula = y ~ s(x, bs='cr', k=4)) +
	# or 2nd-order polynomial on quartiles (as per one of Patrick's other plots)... 
	stat_quantile(formula = y ~ poly(x, 2), lty=2, quantiles=c(0.25, 0.75)) +
	xlab('Scale (km2)') +
	theme_bw()
