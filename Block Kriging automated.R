library(raster)
library(gstat)
library(rgeos)
library(sf)
library(tidyverse)
library(sp)
library(rgdal)

##############################################################################################################################################
Region_grid_empty = raster(xmn= , ymn= , xmx = ,ymx = , resolution = ,
                          crs = "UTM")
Region_grid_empty[] = 0 # set to zero
###############################################################################################################
ls = list.files()
for(i in ls){
  tryCatch({ # This function is used to skip errors if happened and continue to loop through folders
    ############################## Import Raw yield data and its boundaries boundaries #########################################################
    read_csv=list.files(paste0(getwd(),"/",i),pattern=".csv$",full.names=T)
    yield_df = lapply(read_csv,read.csv)
    yield_df = data.frame(yield_df)
    ###############################  Procees yield data to remove outliers  using normalisation function #######################################
    yield_df = subset(yield_df, yield_df$Yield.Raw.tonne.ha. >=0)
    yield_df = subset(yield_df, yield_df$Yield.Raw.tonne.ha. <=10)
    m = mean(yield_df$Yield.Raw.tonne.ha.)
    sd =sd(yield_df$Yield.Raw.tonne.ha.)
    out.u<-m+2*sd # Change to 3 if you want another level of cleaning
    out.l<-m-2*sd # Change to 3 if you want another level of cleaning
    Cleaned = subset(yield_df,yield_df$Yield.Raw.tonne.ha. < out.u & yield_df$Yield.Raw.tonne.ha. > out.l)
    ################################################# Extracting yield from the dataframe #####################################################
    # Cleaned = Cleaned[1:100,]# Use this line to test your code on small portion of data
    Cleaned_2 = dplyr::select ( Cleaned,Longitude,Latitude,Yield.Raw.tonne.ha.)# Use this code to generate new df as some csv file columns are not consistent
    colnames(Cleaned_2)[3] = "Yield" # Here you can change the name of your column name
    ################################################# Convert dataframe spatialDataFrame ######################################################
    coordinates(Cleaned_2) = c("Longitude", "Latitude") # Use this if your Spatial dataframe is long lat
    proj4string(Cleaned_2) = CRS("+init=epsg:4326") # Assign the WGS 84 first
    CRS.new = CRS("")# Obtain the UTM from other data # You can use EPSG: Spatial code of the area of interest
    Cleaned_2 = spTransform(Cleaned_2, CRS.new)# Transform the new SpatialDataFrame
    buffer_polygon = gBuffer(Cleaned_2,width = 5,quadsegs = 10) # Create buffer around the spatial dataframe
    ################################################ Make sure boundaries has same coordinate system  ##########################################
    empty_RS_for_Krig = crop(Region_grid_empty,buffer_polygon)# Crop big raster we created for the entire study area to the paddock extent
    empty_RS_for_Krig = mask(empty_RS_for_Krig,buffer_polygon)# mask the pixels we don't want to krige on them
    SpAsGrid = rasterToPoints(empty_RS_for_Krig, spatial = TRUE)# Convert paddock empty raster to points
    gridded(SpAsGrid) = TRUE # make spatialDataFrame to convert it to spatial pixel
    SpAsGrid = as(SpAsGrid, "SpatialPixels")# And now convert it to SpatialPixel which is used for kriging
    ##########################################################################################################################################
    # Now Data is ready to fit the variogam and perform Block Kriging
    ##########################################################################################################################################
    variog_Yield = variogram(Cleaned_2$Yield~1, Cleaned_2)# fit the model (variogram)
    # plot(variog_Yield)# plot variogram
    variogram_plot = as.data.frame(variog_Yield)
    options(warn = -1) # don't print warnings
    Yield.vgm = fit.variogram(variog_Yield, vgm(c("Exp", "Mat", "Sph", "Gau")), fit.kappa = TRUE)# This will choose the best variogram
    model_fitting = as.data.frame(Yield.vgm)# This csv will be exported to check model later
    Block_Krig = krige(Yield ~1,Cleaned_2,SpAsGrid, maxdist = 100,nmax = 100, Yield.vgm, block = c(20, 20))# Perform block kriging
    ##################################### Convert results to Data frame #######################################################################
    Block_Krig_Df = as.data.frame(Block_Krig, xy =T)# Convert the results to dataframe
    Block_Krig_Df_Pred = Block_Krig_Df[1:3]# Extract predicted yield
    Block_Krig_Df_Standard_Error = Block_Krig_Df[c(1,2,4)]# Extract predicted yield
    # ############################################################################################################################################
    # # Now make raster of the predictions and Variance
    # # 1- Predictions Raster
    coordinates(Block_Krig_Df_Pred) = ~x+y
    gridded(Block_Krig_Df_Pred) = TRUE
    # plot(Block_Krig_Df_Pred)
    Ras_Results = rasterFromXYZ(Block_Krig_Df_Pred)
    names(Ras_Results)[1] = "Yield"
    crs(Ras_Results) = "+init=epsg:XXXX" #epsg will depend on your area of interest
    Ras_Results = clamp(Ras_Results,0,10)
    # 2- Standard Error Raster
    coordinates(Block_Krig_Df_Standard_Error) = ~x+y
    gridded(Block_Krig_Df_Standard_Error) = TRUE
    Ras_Results_Variance = rasterFromXYZ(Block_Krig_Df_Standard_Error)
    names(Ras_Results_Variance)[1] = "Variance"
    crs(Ras_Results_Variance) = "+init=epsg:XXXX" #epsg will depend on your area of interest
    ##############################################################################################################################################
    # # Final step is to export the results
    AOIname = paste0(i) # Get the name of your file to extract the Yield raster name, SE name and the polygon name.
    Aggname = substr(i, 0, 19) # This will depend on how many characters your file name has
    Aggname = gsub(" ", "", Aggname, fixed = TRUE) # Reomve spaces in case the file name has them.
    Farmname1 = strsplit(i, "-", "[", 1)[[1]][2]
    Farmname2 = strsplit(i, "-", "[", 1)[[1]][3]
    Farmname1 = gsub(" ", "", Farmname1, fixed = TRUE)
    Farmname2 = gsub(" ", "", Farmname2, fixed = TRUE)
    YR = str_extract(i, "\\d{4}")
    CropType = substr(AOIname, 1, regexpr("\\[", AOIname)-1)
    CropType = gsub(" ", "", CropType, fixed = TRUE)
    CropType = gsub("-", "_", CropType, fixed = TRUE)
    Crop = sub('.*\\_', '', CropType)
    finalName = str_c(Aggname, Farmname1,Farmname2,YR,Crop, sep = "_")
    
    ## Now export the shapefile of the field boundaries
    poly =  as(buffer_polygon, "SpatialPolygonsDataFrame")
    writeOGR(poly,dsn= paste0(getwd(),"/",finalName, ""),finalName, driver="ESRI Shapefile")
    writeRaster(Ras_Results,paste0(getwd(),"/",finalName, "","/","Yield.tif"))
    writeRaster(Ras_Results_Variance,paste0(getwd(),"/",finalName, "","/","Variance.tif"))
    write.csv(model_fitting,paste0(getwd(),"/",finalName, "","/","_variogram_model.csv"))
    write.csv(variogram_plot,paste0(getwd(),"/",finalName, "","/","_variogram_plot_as_df.csv"))
    cat(i,"\n")
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}



