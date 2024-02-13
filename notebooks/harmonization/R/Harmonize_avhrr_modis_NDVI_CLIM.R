
# This script uses the MGCV Generalised Additive Modelling package to calibrate
# monthly AVHRR NDVI (CDR version) to approximate MODIS NDVI. Predictions are run in
# parallel so computation is fast (~15 mins with 11 cores), but the script requires a
# 'hugemem' queue on the NCI (12 cores, 366 GiB works well), using the R-Studio Rocker Image.
# The input datasets have been split into regions to increase the accuracy of the final predictions - 
# 'trees' and 'nontrees' defined using woody cover fraction.

# This version of the script uses climate features (rainfall,srad) as inputs to the GAM.

# Written by Sami Rifai and Chad Burton Oct. 2023

library(stars); library(tidyverse); library(data.table); library(lubridate)
library(dtplyr, warn.conflicts = FALSE); library(mgcv); library(RcppArmadillo)
library(metrica); library(scico); library(mgcViz); library(foreach); library(doParallel)

##########################################################################
# Analysis Parameters ----------------------------------------------------
# (Adjust parameters to suit)
##########################################################################

n_cores <- 12
base_path <- "/g/data/os22/chad_tmp/AusENDVI/data/NDVI_harmonization/GAM/"
regions = list('trees', 'nontrees')

#loop through the regions and create models/predictions
for (i in regions) {
    
    setDTthreads(threads=0)
    
    ##########################################################################
    # Open datasets and wrangle -----------------------------------------------
    ##########################################################################
    
    # These AVHRR datasets have already been filtered/cleaned in a python script,
    # so the usual data-table filtering has been removed in this script.
    avhrr_path <- paste(base_path, i, '_AVHRR_NDVI_5km_monthly_1982_2013_GAMinput.nc', sep='')
    modis_path <- paste(base_path, i, '_MODIS_NDVI_5km_monthly_200003_202212.nc', sep='')
    
    ## Read AVHRR data 'bands' separately
    tmp_median <- stars::read_ncdf(avhrr_path, var="NDVI_avhrr", make_time = T, proxy=F)
    tmp_median <- tmp_median %>% set_names(c("ndvi_cdr"))
    
    tmp_mod_median <- stars::read_ncdf(avhrr_path, var="NDVI_modis_median", make_time = T, proxy=F)
    tmp_mod_median <- tmp_mod_median %>% set_names(c("ndvi_modis_median"))
    
    tmp_mod_min <- stars::read_ncdf(avhrr_path, var="NDVI_modis_min", make_time = T, proxy=F)
    tmp_mod_min <- tmp_mod_min %>% set_names(c("ndvi_modis_min"))
    
    # tmp_median_1b <- stars::read_ncdf(avhrr_path, var="NDVI_avhrr_1b", make_time = T, proxy=F)
    # tmp_median_1b <- tmp_median_1b %>% set_names(c("ndvi_cdr_1b"))
    # 
    # tmp_median_1f <- stars::read_ncdf(avhrr_path, var="NDVI_avhrr_1f", make_time = T, proxy=F)
    # tmp_median_1f <- tmp_median_1f %>% set_names(c("ndvi_cdr_1f"))
    
    tmp_mod_max <- stars::read_ncdf(avhrr_path, var="NDVI_modis_max", make_time = T, proxy=F)
    tmp_mod_max <- tmp_mod_max %>% set_names(c("ndvi_modis_max"))
    
    tmp_rain_cml3 <- stars::read_ncdf(avhrr_path, var="rain_cml3", make_time = T, proxy=F)
    tmp_rain_cml3 <- tmp_rain_cml3 %>% set_names(c("rain_cml3"))
    
    tmp_srad <- stars::read_ncdf(avhrr_path, var="srad", make_time = T, proxy=F)
    tmp_srad <- tmp_srad %>% set_names(c("srad"))
    
    tmp_sza <- stars::read_ncdf(avhrr_path, var='SZEN_median', make_time = T, proxy=F)
    tmp_sza <- tmp_sza %>% set_names(c('sza'))
    
    tmp_tod <- stars::read_ncdf(avhrr_path, var='TIMEOFDAY_median', make_time = T, proxy=F)
    tmp_tod <- tmp_tod %>% set_names(c('tod'))
    
    ## Convert AVHRR data into data tables, add 'month' var.
    tmp <- c(tmp_median, tmp_sza, tmp_tod, tmp_mod_median,# tmp_median_1b, tmp_median_1f,
                 tmp_mod_min, tmp_mod_max, tmp_rain_cml3, tmp_srad)
    
    d_cdr <- tmp %>% units::drop_units() %>% as.data.table()
    d_cdr <- d_cdr %>% mutate(month=month(time)) %>% as.data.table()
    
    #clear up some memory
    rm(tmp_median, tmp_sza, tmp_tod, tmp_mod_median, # tmp_median_1b, tmp_median_1f,
                 tmp_mod_min, tmp_mod_max, tmp_rain_cml3, tmp_srad)
    
    gc(reset = T, full=T)
    
    # Import MODIS NDVI--------------------------------------------------------
    tmp_nm <- stars::read_ncdf(modis_path, var="NDVI_median", make_time = T, proxy=F)
    tmp_nm <- tmp_nm %>% set_names(c('ndvi_mcd'))
    
    d_mcd <- tmp_nm %>% units::drop_units() %>% as.data.table()
    d_mcd <- d_mcd[between(time,ymd("2000-03-01"),ymd("2013-12-31"))==T]
    rm(tmp_nm)
    gc(reset = T, full=T)
    
    ##########################################################################
    # Calibrate AVHRR CDR to approximate MODIS NDVI --------------------------
    ##########################################################################
    
    # pre-setting the keys hypothetically makes the merge faster
    setkeyv(d_cdr,cols=c("longitude","latitude","time"))
    setkeyv(d_mcd,cols=c("longitude","latitude","time"))
    
    # merge the modis and avhrr datasets # ~ 2 minutes
    dc2 <- merge(d_mcd[,.(longitude,latitude,time,ndvi_mcd)],
                 d_cdr[,.(longitude,latitude, time, ndvi_cdr, ndvi_modis_min, #ndvi_cdr_1b, ndvi_cdr_1f,
                        ndvi_modis_max, sza, ndvi_modis_median, rain_cml3, srad, tod)], 
               all=TRUE,
               by=c("longitude","latitude","time"))
    dc2 <- dc2[,`:=`(month=month(time))]
    
    # training and testing samples
    set.seed(3)
    dc2_train <- dc2[is.na(ndvi_mcd)==F][is.na(ndvi_cdr)==F][
      between(time,ymd("2000-03-01"),ymd("2013-12-31"))==T][sample(.N, 1e6)]
    
    # Create a BAM-GAM with climate terms    
    mc11 <- bam(ndvi_mcd ~ 
            s(ndvi_cdr, bs='cs')+
            s(rain_cml3, bs='ts')+
            s(srad, bs='ts')+
            s(ndvi_modis_median, bs='ts')+
            s(ndvi_modis_min, bs='ts')+
            s(ndvi_modis_max, bs='ts')+
            ti(month, sza, bs=c('cc','ts'))+
            ti(ndvi_cdr, tod, bs='ts')+
            te(longitude, latitude, by=ndvi_cdr),
          data=dc2_train[ndvi_mcd>0.05][ndvi_modis_median > 0],
          discrete=T, 
          select=T)
  
    ##########################################################################
    # Apply Calibration prediction in parallel -------------------------------
    ##########################################################################
    cl <- makeCluster(n_cores, outfile="")
    
    dc2 <- mutate(dc2, year=data.table::year(time),
                month=data.table::month(time)) %>% 
                as.data.table()
    
    # na's in the predictor cols can break stuff
    tmp <- dc2 %>% select(time,month,year,
      longitude,latitude,ndvi_modis_median,ndvi_cdr,ndvi_modis_max,# ndvi_cdr_1b, ndvi_cdr_1f,
       ndvi_modis_min, month, sza, tod, srad, rain_cml3) %>% 
      drop_na() %>% as.data.table()
    
    tmp[,proc_id := 1:nrow(tmp)]
    
    # increase number of chunks
    vec_ids <- tmp$proc_id
    n_chunks <- n_cores*10
    l_ids <- split(vec_ids, cut(seq_along(vec_ids), n_chunks, labels = FALSE))
    
    out <- foreach(i = 1:length(l_ids), 
                 .packages = c("mgcv","data.table","tidyverse"),
                 .combine=rbind) %dopar% {
                   out <- tmp[proc_id%in%l_ids[[i]]] %>% 
                     mutate(ndvi_mcd_pred = predict(mc11, 
                                                    newdata=., 
                                                    newdata.guaranteed = T,
                                                    discrete = TRUE)) %>% 
                     as.data.table() # force computation
                   gc(full=T, reset=T)
                   out
                 } 
    ## ---prepare for export-------
    setkeyv(out, c("longitude","latitude",'time'))
    
    d_export <- merge(dc2,
                    out[, .(longitude,latitude, time, 
                            ndvi_mcd_pred)],
                    by=c("longitude","latitude",'time'),
                    all=TRUE)
    
    tmp3 <- st_as_stars(d_export, dims = c("longitude","latitude","time"))
    
    ## requires stars 0.6-1 or greater
    stars::write_mdim(tmp3,
                      filename=
                        paste(base_path,
                            'NDVI_', i, '_CLIM_GAM_harmonize_5km_monthly_1982_2013.nc',
                            sep=''),
                      layer = c("ndvi_mcd", "ndvi_cdr", 
                                "ndvi_mcd_pred", 
                                "month", "year"
                      ))
    
    ##########################################################################
    ## Plots, check residuals through time, space, etc -----------------------
    ##########################################################################
    
    # dc2_test <- dc2[is.na(ndvi_mcd)==F][is.na(ndvi_cdr)==F][
    #   between(time,ymd("2000-03-01"),ymd("2013-12-31"))==T][sample(.N, 1e6)]
    
    # summarize and test results
    # summary(mc11) # 94.5%
    # getViz(mc11) %>% 
    #   plot(allTerms=T) %>% 
    #   print()
    
    # dc2_test %>% 
    #   mutate(pred = predict(mc11,newdata=.,type='response')) %>% 
    #   select(ndvi_mcd,pred) %>% 
    #   drop_na() %>% as.data.table() %>% 
    #   .[,.(rmse = metrica::RMSE(
    #     obs=.$ndvi_mcd,
    #     pred=.$pred), 
    #     R2 = metrica::R2(
    #       obs=.$ndvi_mcd,
    #       pred=.$pred))] 
    # 
    # #plot 1:1 line
    # dc2_test %>% 
    #   mutate(pred = predict(mc11,newdata=.,type='response')) %>% 
    #   select(ndvi_mcd,pred) %>% 
    #   drop_na() %>% as.data.table() %>%
    #   .[sample(.N,1000)] %>% 
    #   ggplot(aes(pred,ndvi_mcd))+
    #   geom_point()+
    #   geom_abline(col='red') +
    #   coord_equal()
    # 
    # #plot difference between modis/avhrr
    # dc2_test[sample(.N,10000)] %>% 
    #   mutate(pred = predict(mc11,newdata=.,type='response')) %>% 
    #   select(longitude, latitude, ndvi_mcd,pred) %>% 
    #   drop_na() %>% as.data.table() %>%
    #   ggplot(aes(longitude, latitude, color=ndvi_mcd - pred))+
    #   geom_point()+
    #   scico::scale_color_scico(palette='roma',midpoint=0, 
    #                            limits=c(-0.2,0.2),
    #                            oob=scales::squish)+
    #   coord_sf()
    # 
    
    # ss1 <- d_export[is.na(ndvi_mcd)==F][,.(
    #   res = mean(ndvi_mcd - ndvi_mcd_pred,na.rm=T),
    #   res_sd =- sd(ndvi_mcd - ndvi_mcd_pred,na.rm=T)
    # ),by=time]
    # 
    # ss1 %>%
    #   ggplot(aes(time, res))+
    #   geom_hline(yintercept = 0, lwd=1, col='grey30')+
    #   geom_ribbon(aes(time,ymin=res+res_sd,ymax=res-res_sd),
    #               alpha=0.3,col='transparent')+
    #   geom_line(col='blue') +
    #   labs(y = expression(paste(NDVI[MCD] - NDVI[pred])))
    
    # tmp[,.(SZA = median(sza,na.rm=T)),by=time] %>%
    #   ggplot(aes(time, SZA))+
    #   geom_line()
    # tmp[,.(TOD = median(tod,na.rm=T)),by=time] %>%
    #   ggplot(aes(time, TOD))+
    #   geom_line()
    # 
    # ## check residuals through time
    # ss2 <- d_export[is.na(ndvi_mcd)==F][,.(
    #   res = mean(ndvi_mcd - ndvi_mcd_pred,na.rm=T),
    #   ndvi_u = mean(ndvi_mcd,na.rm=T),
    #   res_sd =- sd(ndvi_mcd - ndvi_mcd_pred,na.rm=T)
    # ),by=.(longitude,latitude)]
    # 
    # ss2 %>%
    #   ggplot(aes(longitude, latitude, fill=100 * res/ndvi_u))+
    #   geom_raster()+
    #   coord_sf() +
    #   scico::scale_fill_scico(palette='roma',midpoint=0,
    #                           limits=c(-25, 25),
    #                           oob=scales::squish) +
    #   labs(fill = "% bias")

}


