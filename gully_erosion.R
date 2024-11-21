# Title: 

# Author details: 

# Script and data info: This script performs the  (XXX et al.,) 

# Copyright statement: This script is the product of the work of 

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
# OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

# Packages used

install.packages("whitebox")
whitebox::install_whitebox()
install.packages("terra")
install.packages("data.table")
install.packages("splitstackshape")
install.packages("car")
install.packages("groupdata2")
install.packages("splitstackshape")
install.packages("earth")
install.packages("caret")
install.packages("pROC")

library(whitebox)
library(terra)
library(data.table)
library(splitstackshape)
library(car)
library(groupdata2) 
library(splitstackshape) 
library(earth) 
library(caret) 
library(pROC) 



###################################################################################################################
# START OF THE SCRIPT
###################################################################################################################

# Set working directory to your location
wd <- ("D:/GRAZIA/PHD/TURKEY/R")
setwd(wd)

#### INPUT DATA PREPARATION ####
##### Upload DEM #####
# resampling 
wbt_resample(
  input='./lidar_USDA_clip.tif',       # Path to the input raster file
  output='./dem_4m.tif',               # Path to save the output resampled raster
  cell_size = 4,                       # Desired cell size for the output raster 
  method = "bilinear",                 # Resampling method: 'nn' (nearest neighbor), 'bilinear', or 'cc' (cubic convolution)
  wd = NULL,                           # Change the working directory (optional)
  verbose_mode = NULL,                 # Enable verbose output for debugging (optional)
  compress_rasters = NULL,             # Option to compress the output raster file (optional)
  command_only = FALSE                 # Set to TRUE to return the command as a string without executing it
)


##### Upload the vector of the Study Area  #####
# convert to raster
wbt_vector_polygons_to_raster(
  input="./tiles.shp",                 # Path to the input shapefile containing vector polygons
  output="./tiles.tif",                # Path to save the output raster file
  field = "OBJECTID",                  # Attribute field from the shapefile to be rasterized
  nodata = TRUE,                       # Enable NoData values for areas without polygons; default is 0.0 if not set
  cell_size = NULL,                    # Define the resolution of the raster; uses the base raster if not specified
  base = "./dem_4m.tif",               # Path to a base raster file for spatial reference and resolution
  wd = NULL,                           # Change the working directory (optional)
  verbose_mode = NULL,                 # Enable verbose output for debugging (optional)
  compress_rasters = NULL,             # Option to compress the output raster file (optional)
  command_only = FALSE                 # Set to TRUE to return the command as a string without running it
)


##### Upload the vector of the Gullies #####
# Convert the vector lines of gullies into a raster format
wbt_vector_lines_to_raster(
  input="./G0_linear.shp",            # Path to the input shapefile containing gully lines
  output="./G0.tif",                  # Path to save the output raster file
  field = "IDN",                      # Attribute field to use for rasterization
  nodata = TRUE,                      # Enable NoData values for cells with no corresponding vector data
  cell_size = NULL,                   # Specify the resolution of the raster; uses the base raster if not specified
  base = "./dem_4m.tif",              # Path to a base raster file for spatial reference and resolution
  wd = NULL,                          # Change the working directory (optional)
  verbose_mode = NULL,                # Enable verbose output for debugging (optional)
  compress_rasters = NULL,            # Option to compress the output raster file (optional)
  command_only = FALSE                # Set to TRUE to return the command as a string without running it
)

# Convert the gully raster into vector points
wbt_raster_to_vector_points(
  './G0.tif',                         # Path to the input raster file (gullies)
  './G0.shp',                         # Path to save the output vector point shapefile
  wd = NULL,                          # Change the working directory (optional)
  verbose_mode = NULL,                # Enable verbose output for debugging (optional)
  compress_rasters = NULL,            # Option to compress the output shapefile (optional)
  command_only = FALSE                # Set to TRUE to return the command as a string without executing it
)


#### DEM PRE-PROCESSING, CHANNELS EXTRACTIONS AND GULLY SNAPPING ####
##### Fill depression - Wang & Liu #####
wbt_fill_depressions_wang_and_liu(
  dem='./dem4m.tif',                  # Path to the input DEM raster file
  output='/fill_wang.tif',            # Path to save the output DEM with depressions filled
  fix_flats = TRUE,                   # Optional flag indicating whether flat areas should have a small gradient applied 
  flat_increment = NULL,              # Specify the elevation increment for flat areas (optional, uses default if not set)
  wd = NULL,                          # Change the working directory (optional)
  verbose_mode = NULL,                # Enable verbose output for debugging (optional)
  compress_rasters = NULL,            # Option to compress the output raster file (optional)
  command_only = FALSE                # Set to TRUE to return the command as a string without executing it
)


##### Extracting D8 Contributing Area #####
# Contributing Area (CA m²) represents the total upstream area contributing flow to a point. 
# Calculate the D8 contributing area (CA) from the input DEM
wbt_d8_flow_accumulation(
  input = '/fill_wang.tif',           # Path to the input filled DEM raster file
  output = './CA.tif',                # Path to save the output raster showing the total contributing area
  out_type = 'area'                   # Output type set to total contributing area, measuring the total upstream flow area
)

##### Extracting D8 Specific Contributing Area #####
# Specific Contributing Area (SCA m²/m)) normalizes CA by flow width, accounting for slope.
# Calculate the D8 specific contributing area (SCA) from the input DEM
wbt_d8_flow_accumulation(
  input = '/fill_wang.tif',               # Path to the input filled DEM raster file
  output = './SCA.tif',                   # Path to save the output raster showing the specific contributing area
  out_type = 'Specific Contributing Area' # Output type set to specific contributing area, measuring flow per unit width
)


##### Extracting Channels with Different Thresholds #####

# Set stream extraction thresholds based on contributing area (CA)
T1 <- 1000                             # Threshold 1: Minimum contributing area (in m²) to define a stream
T2 <- 2000                             # Threshold 2: A higher threshold to define streams with larger contributing areas
T3 <- 5000                             # Threshold 3: An even larger threshold for major streams

# Extract streams using Threshold 1 (T1)
wbt_extract_streams(
  './CA.tif',                          # Path to the input Contributing Area (CA) raster file
  './CH_T1.tif',                       # Path to save the output raster with streams extracted using T1
  T1,                                  # Threshold value for contributing area to define streams
  zero_background = FALSE,             # Background values will not be set to zero (optional)
  wd = NULL,                           # Specify the working directory (optional)
  verbose_mode = NULL,                 # Enable verbose output for debugging (optional)
  compress_rasters = NULL,             # Option to compress the output raster file (optional)
  command_only = FALSE                 # Set to TRUE to return the command as a string without executing it
)

# Extract streams using Threshold 2 (T2)
wbt_extract_streams(
  './CA.tif',                          # Path to the input Contributing Area (CA) raster file
  './CH_T2.tif',                       # Path to save the output raster with streams extracted using T2
  T2,                                  # Threshold value for contributing area to define streams
  zero_background = FALSE,             # Background values will not be set to zero (optional)
  wd = NULL,                           # Specify the working directory (optional)
  verbose_mode = NULL,                 # Enable verbose output for debugging (optional)
  compress_rasters = NULL,             # Option to compress the output raster file (optional)
  command_only = FALSE                 # Set to TRUE to return the command as a string without executing it
)

# Extract streams using Threshold 3 (T3)
wbt_extract_streams(
  './CA.tif',                          # Path to the input Contributing Area (CA) raster file
  './CH_T3.tif',                       # Path to save the output raster with streams extracted using T3
  T3,                                  # Threshold value for contributing area to define streams
  zero_background = FALSE,             # Background values will not be set to zero (optional)
  wd = NULL,                           # Specify the working directory (optional)
  verbose_mode = NULL,                 # Enable verbose output for debugging (optional)
  compress_rasters = NULL,             # Option to compress the output raster file (optional)
  command_only = FALSE                 # Set to TRUE to return the command as a string without executing it
)

##### Snapping Gully Points Using Different Distances and Channels #####

# Set snapping distances (SD) for aligning gully points to the nearest streams
SD1 <- 4                              # Snapping Distance 1: Snap points within a distance of 4 map units
SD2 <- 8                              # Snapping Distance 2: Snap points within a distance of 8  map units
SD3 <- 16                             # Snapping Distance 3: Snap points within a distance of 16  map units

# Snapping gully points to streams using SD1 (distance = 4)
wbt_jenson_snap_pour_points(
  pour_pts = "./G0.shp",              # Path to the input shapefile of gully points
  streams = './CH_T1.tif',            # Path to the input stream raster (Threshold T1)
  output = "./G1.shp",                # Path to save the output snapped points for T1
  snap_dist = SD1,                    # Snapping distance set to SD1 (4 units)
  wd = NULL,                          # Specify the working directory (optional)
  verbose_mode = NULL,                # Enable verbose output for debugging (optional)
  compress_rasters = NULL,            # Option to compress the output file (optional)
  command_only = FALSE                # Set to TRUE to return the command as a string without executing it
)

wbt_jenson_snap_pour_points(
  pour_pts = "./G0.shp",
  streams = './CH_T2.tif',            # Use stream raster for Threshold T2
  output = "./G2.shp",                # Output snapped points for T2
  snap_dist = SD1,
  wd = NULL,
  verbose_mode = NULL,
  compress_rasters = NULL,
  command_only = FALSE
)

wbt_jenson_snap_pour_points(
  pour_pts = "./G0.shp",
  streams = './CH_T3.tif',            # Use stream raster for Threshold T3
  output = "./G3.shp",                # Output snapped points for T3
  snap_dist = SD1,
  wd = NULL,
  verbose_mode = NULL,
  compress_rasters = NULL,
  command_only = FALSE
)

# Snapping gully points to streams using SD2 (distance = 8)
wbt_jenson_snap_pour_points(
  pour_pts = "./G0.shp",
  streams = './CH_T1.tif',
  output = "./G4.shp",                # Output snapped points for T1 using SD2
  snap_dist = SD2,                    # Snapping distance set to SD2 (8 units)
  wd = NULL,
  verbose_mode = NULL,
  compress_rasters = NULL,
  command_only = FALSE
)

wbt_jenson_snap_pour_points(
  pour_pts = "./G0.shp",
  streams = './CH_T2.tif',
  output = "./G5.shp",                # Output snapped points for T2 using SD2
  snap_dist = SD2,
  wd = NULL,
  verbose_mode = NULL,
  compress_rasters = NULL,
  command_only = FALSE
)

wbt_jenson_snap_pour_points(
  pour_pts = "./G0.shp",
  streams = './CH_T3.tif',
  output = "./G6.shp",                # Output snapped points for T3 using SD2
  snap_dist = SD2,
  wd = NULL,
  verbose_mode = NULL,
  compress_rasters = NULL,
  command_only = FALSE
)

# Snapping gully points to streams using SD3 (distance = 16)
wbt_jenson_snap_pour_points(
  pour_pts = "./G0.shp",
  streams = './CH_T1.tif',
  output = "./G7.shp",                # Output snapped points for T1 using SD3
  snap_dist = SD3,                    # Snapping distance set to SD3 (16 units)
  wd = NULL,
  verbose_mode = NULL,
  compress_rasters = NULL,
  command_only = FALSE
)

wbt_jenson_snap_pour_points(
  pour_pts = "./G0.shp",
  streams = './CH_T2.tif',
  output = "./G8.shp",                # Output snapped points for T2 using SD3
  snap_dist = SD3,
  wd = NULL,
  verbose_mode = NULL,
  compress_rasters = NULL,
  command_only = FALSE
)

wbt_jenson_snap_pour_points(
  pour_pts = "./G0.shp",
  streams = './CH_T3.tif',
  output = "./G9.shp",                # Output snapped points for T3 using SD3
  snap_dist = SD3,
  wd = NULL,
  verbose_mode = NULL,
  compress_rasters = NULL,
  command_only = FALSE
)

##### Create predictor Variables #####

##### CA Independent Variables #####

# Extracting slope angle from the filled DEM (using Wang & Liu method)
wbt_slope(
  '/fill_wang.tif',                    # Input filled DEM 
  './slope_deg.tif',                   # Output raster file for slope angle in degrees
  zfactor = NULL,                      # Optional: Vertical scaling factor (leave NULL to use default)
  units = "degrees",                   # Specify units for slope angle (degrees)
  wd = NULL,                           # Optional: Specify working directory
  verbose_mode = NULL,                 # Optional: Enable verbose output for debugging
  compress_rasters = NULL,             # Optional: Option to compress the output raster
  command_only = FALSE                 # Set to TRUE to return the command as a string without execution
)

# Extracting plan curvature from the filled DEM
wbt_plan_curvature( 
  '/fill_wang.tif',                    # Input filled DEM 
  './PLANC.tif',                       # Output raster file for plan curvature
  log = FALSE,                         # Optional: Logarithmic transformation (FALSE to skip)
  zfactor = NULL,                      # Optional: Vertical scaling factor (leave NULL for default)
  wd = NULL,                           # Optional: Working directory
  verbose_mode = NULL,                 # Optional: Enable verbose output
  compress_rasters = NULL,             # Optional: Compress the output raster
  command_only = FALSE                 # Set to TRUE to generate command without execution
)

# Extracting deviation from mean elevation (a measure of local topographic variation)
wbt_dev_from_mean_elev(
  '/fill_wang.tif',                    # Input filled DEM 
  './dev_mean_ele.tif',                # Output raster file for deviation from mean elevation
  filterx = 11,                        # Filter size in x direction (size of the neighborhood for computation)
  filtery = 11,                        # Filter size in y direction
  wd = NULL,                           # Optional: Working directory
  verbose_mode = NULL,                 # Optional: Enable verbose output for debugging
  compress_rasters = NULL,             # Optional: Compress output raster file
  command_only = FALSE                 # Set to TRUE to generate the command as a string without execution
)

# Extracting elevation percentile (a statistical measure of local elevation)
wbt_elev_percentile(
  '/fill_wang.tif',                    # Input filled DEM 
  './elev_perc.tif',                   # Output raster file for elevation percentile
  filterx = 11,                        # Filter size in x direction
  filtery = 11,                        # Filter size in y direction
  sig_digits = 2,                      # Optional: Number of significant digits for output
  wd = NULL,                           # Optional: Working directory
  verbose_mode = NULL,                 # Optional: Enable verbose output
  compress_rasters = NULL,             # Optional: Compress the output raster
  command_only = FALSE                 # Set to TRUE to return the command as a string without execution
)

# Extracting geomorphons (topographic features or landforms) from the filled DEM
wbt_geomorphons( 
  '/fill_wang.tif',                    # Input filled DEM 
  './geomorphons.tif',                 # Output raster file for geomorphons
  search = 50,                         # Search radius for geomorphon calculation
  threshold = 0,                       # Threshold for form classification
  fdist = 0,                           # Optional: Distance factor (0 for no influence)
  skip = 0,                            # Optional: Number of classes to skip
  forms = TRUE,                        # Set to TRUE to output geomorphon forms (shapes)
  residuals = FALSE,                   # Set to TRUE to output residuals (optional)
  wd = NULL,                           # Optional: Working directory
  verbose_mode = NULL,                 # Optional: Enable verbose output for debugging
  compress_rasters = NULL,             # Optional: Compress the output raster file
  command_only = FALSE                 # Set to TRUE to return the command as a string without execution
)

##### CA Dependent Variables #####

# Extracting GORD (D8 pointer) from the filled DEM to calculate flow directions
wbt_d8_pointer( 
    '/fill_wang.tif',                   # Input filled DEM 
    './pntr.tif',                       # Output raster file for D8 pointer (flow direction map)
    esri_pntr = FALSE,                  # Optional: Set to TRUE for ESRI-style pointers, FALSE for default
    wd = NULL,                          # Optional: Working directory
    verbose_mode = NULL,                # Optional: Enable verbose output for debugging
    compress_rasters = NULL,            # Optional: Compress the output raster
    command_only = FALSE                # Set to TRUE to generate the command as a string without execution
)

# Extracting streams using a threshold of 0 from the contributing area raster (CA)
wbt_extract_streams( 
    './CA.tif',                         # Input contributing area (CA) raster
    './CH0m.tif',                       # Output raster for extracted streams with threshold 0
    threshold = 0,                      # Set threshold for stream extraction (0 means all cells are included)
    zero_background = FALSE,            # Set to TRUE to set no data values to zero, FALSE to retain original no data values
    wd = NULL,                          # Optional: Working directory
    verbose_mode = NULL,                # Optional: Enable verbose output for debugging
    compress_rasters = NULL,            # Optional: Compress the output raster
    command_only = FALSE                # Set to TRUE to return the command only, without execution
)

# Extracting the Strahler stream order from the flow direction pointers and streams
wbt_strahler_stream_order( 
    './pntr.tif',                       # Input D8 flow direction pointer raster
    './CH0m.tif',                       # Input stream raster (from the previous step)
    './strahler.tif',                   # Output raster for Strahler stream order
    esri_pntr = FALSE,                  # Optional: Set to TRUE for ESRI-style pointers, FALSE for default
    zero_background = FALSE,            # Set to TRUE to assign zero to no data values in the output
    wd = NULL,                          # Optional: Working directory
    verbose_mode = NULL,                # Optional: Enable verbose output for debugging
    compress_rasters = NULL,            # Optional: Compress the output raster
    command_only = FALSE                # Set to TRUE to return the command only, without execution
)

# Extracting Stream Power Index (SPI) from the specific contributing area (SCA) and slope rasters
wbt_stream_power_index(
    './SCA.tif',                        # Input raster for Specific Contributing Area (SCA)
    './slope_deg.tif',                  # Input raster for slope angle (in degrees)
    './SPI.tif',                        # Output raster for Stream Power Index (SPI)
    exponent = 1,                       # Exponent for SCA (default value is 1)
    wd = NULL,                          # Optional: Working directory
    verbose_mode = NULL,                # Optional: Enable verbose output for debugging
    compress_rasters = NULL,            # Optional: Compress the output raster
    command_only = FALSE                # Set to TRUE to return the command only, without execution
)

# Extracting Length Slope Factor (LSF) as part of the Sediment Transport Index (STI)
wbt_sediment_transport_index(
    './SCA.tif',                        # Input raster for Specific Contributing Area (SCA)
    './slope_deg.tif',                  # Input raster for slope angle (in degrees)
    './LS.tif',                         # Output raster for Length Slope Factor (LS)
    sca_exponent = 0.4,                 # Exponent for SCA (default value is 0.4)
    slope_exponent = 1.3,               # Exponent for slope (default value is 1.3)
    wd = NULL,                          # Optional: Working directory
    verbose_mode = NULL,                # Optional: Enable verbose output for debugging
    compress_rasters = NULL,            # Optional: Compress the output raster
    command_only = FALSE                # Set to TRUE to return the command only, without execution
)

# Extracting Topographic Wetness Index (TWI) based on SCA and slope angle rasters
wbt_wetness_index( 
    './SCA.tif',                        # Input raster for Specific Contributing Area (SCA)
    './slope_deg.tif',                  # Input raster for slope angle (in degrees)
    './TWI.tif',                        # Output raster for Topographic Wetness Index (TWI)
    wd = NULL,                          # Optional: Working directory
    verbose_mode = NULL,                # Optional: Enable verbose output for debugging
    compress_rasters = NULL,            # Optional: Compress the output raster
    command_only = FALSE                # Set to TRUE to return the command only, without execution
)


##### Dataframe preparation #####

# loading gully grids
G0<-rast("./G0.tif")

G1<-vect("./G1.shp")
G1<-rasterize(G1, field="VALUE",G0)
G2<-vect("./G2.shp")
G2<-rasterize(G2, field="VALUE",G0)
G3<-vect("./G3.shp")
G3<-rasterize(G3, field="VALUE",G0)
G4<-vect("./G4.shp")
G4<-rasterize(G4, field="VALUE",G0)
G5<-vect("./G5.shp")
G5<-rasterize(G5, field="VALUE",G0)
G6<-vect("./G6.shp")
G6<-rasterize(G6, field="VALUE",G0)
G7<-vect("./G7.shp")
G7<-rasterize(G7, field="VALUE",G0)
G8<-vect("./G8.shp")
G8<-rasterize(G8, field="VALUE",G0)
G9<-vect("./G9.shp")
G9<-rasterize(G9, field="VALUE",G0)

#set the reference system for all the raster
crs(G0)<-"EPSG:26914"
crs(G2)<-"EPSG:26914"
crs(G1)<-"EPSG:26914"
crs(G3)<-"EPSG:26914"
crs(G4)<-"EPSG:26914"
crs(G5)<-"EPSG:26914"
crs(G6)<-"EPSG:26914"
crs(G7)<-"EPSG:26914"
crs(G8)<-"EPSG:26914"
crs(G9)<-"EPSG:26914"

# loading channels rasters
CH_T1<-rast("./CH_T1.tif")
CH_T2<-rast("./CH_T2.tif")
CH_T3<-rast("./CH_T3.tif")
crs(CH_T1)<-"EPSG:26914"
crs(CH_T2)<-"EPSG:26914"
crs(CH_T3)<-"EPSG:26914"

# loading study area raster
tiles<-rast("./tiles.tif")
crs(tiles)<-"EPSG:26914"

# loading elevation raster
ELE<-rast("./dem_4m.tif")
crs(ELE)<-"EPSG:26914"

# loading contributing area raster
CA<-rast("./CA.tif")
crs(CA)<-"EPSG:26914"

# loading predictors rasters
SLO<-rast("./slope_deg.tif")
crs(SLO)<-"EPSG:26914"
PLC<-rast("./PLANC.tif")
crs(PLC)<-"EPSG:26914"
DEV<-rast("./dev_mean_ele.tif")
crs(DEV)<-"EPSG:26914"
EP<-rast("./elev_perc.tif")
crs(EP)<-"EPSG:26914"
GEO<-rast("./geomorphons.tif")
crs(GEO)<-"EPSG:26914"
GORD<-rast("./strahler.tif")
crs(GORD)<-"EPSG:26914"
SPI<-rast("./SPI.tif")
crs(SPI)<-"EPSG:26914"
TWI<-rast("./TWI.tif")
crs(TWI)<-"EPSG:26914"
LSF<-rast("./LS.tif")
crs(LSF)<-"EPSG:26914"

# rasters stacking
rasters<-c(G0,G1,G2,G3,G4,G5,G6,G7,G8,G9,CH_T1,CH_T2,CH_T3,tiles,ELE,CA,SLO,PLC,DEV,EP,GEO,GORD,SPI,TWI,LSF)

# dataframe creation
data <- as.data.frame(rasters, xy=TRUE)
colnames(data) <- c("x","y","G0","G1","G2","G3","G4","G5","G6","G7","G8","G9","CH_T1","CH_T2","CH_T3","tiles","ELE","CA","SLO","PLC","DEV","EP","GEO","GORD","SPI","TWI","LSF")

#transform the NA value on zero value
data$G0[is.na(data$G0)] = 0
data$G1[is.na(data$G1)] = 0
data$G2[is.na(data$G2)] = 0
data$G3[is.na(data$G3)] = 0
data$G4[is.na(data$G4)] = 0
data$G5[is.na(data$G5)] = 0
data$G6[is.na(data$G6)] = 0
data$G7[is.na(data$G7)] = 0
data$G8[is.na(data$G8)] = 0
data$G9[is.na(data$G9)] = 0
data$CH_T1[is.na(data$CH_T1)] = 0
data$CH_T2[is.na(data$CH_T2)] = 0
data$CH_T3[is.na(data$CH_T3)] = 0

# Convert the values in the 'G0' column to a binary factor with levels "0" and "1"
# If the value in 'G0' is 0, it becomes "0"; otherwise, it becomes "1"
data$G0<-factor( ifelse(data$G0 == 0, "0", "1") )
# Repeat the same transformation for each of the following columns (G1 to G9)
# This will convert each column into a binary factor with levels "0" and "1"
#data$G1<-factor( ifelse(data$G1 == 0, "0", "1") )
#data$G2<-factor( ifelse(data$G2 == 0, "0", "1") )
#data$G3<-factor( ifelse(data$G3 == 0, "0", "1") )
#data$G4<-factor( ifelse(data$G4 == 0, "0", "1") )
#data$G5<-factor( ifelse(data$G5 == 0, "0", "1") )
#data$G6<-factor( ifelse(data$G6 == 0, "0", "1") )
#data$G7<-factor( ifelse(data$G7 == 0, "0", "1") )
#data$G8<-factor( ifelse(data$G8 == 0, "0", "1") )
#data$G9<-factor( ifelse(data$G9 == 0, "0", "1") )

data$CH_T1<-as.factor(data$CH_T1)
data$CH_T2<-as.factor(data$CH_T2)
data$CH_T3<-as.factor(data$CH_T3)
data$tiles<-as.factor(data$tiles)
data$GEO<-as.factor(data$GEO)

# Remove the missing values (NA)
data<-na.omit(data)
# Drop unused levels from all factor variables in the data frame 'data'
data<-droplevels(data)

# Convert your existing dataframe to a data.table
setDT(data)

# Create the new column 'CA_threshold'
# the new column stores the 'CA' values associated with the maximum 'ELE' for each group in G0, G1, ..., G9
data[, G0_CA_threshold := CA[which.max(ELE)], by = G0]
data[, G1_CA_threshold := CA[which.max(ELE)], by = G1]
data[, G2_CA_threshold := CA[which.max(ELE)], by = G2]
data[, G3_CA_threshold := CA[which.max(ELE)], by = G3]
data[, G4_CA_threshold := CA[which.max(ELE)], by = G4]
data[, G5_CA_threshold := CA[which.max(ELE)], by = G5]
data[, G6_CA_threshold := CA[which.max(ELE)], by = G6]
data[, G7_CA_threshold := CA[which.max(ELE)], by = G7]
data[, G8_CA_threshold := CA[which.max(ELE)], by = G8]
data[, G9_CA_threshold := CA[which.max(ELE)], by = G9]

# Transform each column G0, G1, G2, ..., G9 into a binary factor with levels '1' and '0'
# Set to '1' if the value of the column G0, G1, G2,..., G9 is not 0 and CA is greater than or equal to the threshold G0_CA, otherwise '0'
data$G0<-factor( ifelse(data$G0 != 0 & data$CA>=data$G0_CA_threshold, "1", "0") )
data$G1<-factor( ifelse(data$G1 != 0 & data$CA>=data$G1_CA_threshold, "1", "0") )
data$G2<-factor( ifelse(data$G2 != 0 & data$CA>=data$G2_CA_threshold, "1", "0") )
data$G3<-factor( ifelse(data$G3 != 0 & data$CA>=data$G3_CA_threshold, "1", "0") )
data$G4<-factor( ifelse(data$G4 != 0 & data$CA>=data$G4_CA_threshold, "1", "0") )
data$G5<-factor( ifelse(data$G5 != 0 & data$CA>=data$G5_CA_threshold, "1", "0") )
data$G6<-factor( ifelse(data$G6 != 0 & data$CA>=data$G6_CA_threshold, "1", "0") )
data$G7<-factor( ifelse(data$G7 != 0 & data$CA>=data$G7_CA_threshold, "1", "0") )
data$G8<-factor( ifelse(data$G8 != 0 & data$CA>=data$G8_CA_threshold, "1", "0") )
data$G9<-factor( ifelse(data$G9 != 0 & data$CA>=data$G9_CA_threshold, "1", "0") )

# For each column G1 to G9, convert the values into a binary factor (1 or 0)
# The condition is: set to '1' if the value in the column G is not 0 and the corresponding CH_Tx is equal to 1, otherwise set to '0'
data$G1<-factor( ifelse(data$G1 != 0 & data$CH_T1 == 1, "1", "0") )
data$G2<-factor( ifelse(data$G2 != 0 & data$CH_T2 == 1, "1", "0") )
data$G3<-factor( ifelse(data$G3 != 0 & data$CH_T3 == 1, "1", "0") )
data$G4<-factor( ifelse(data$G4 != 0 & data$CH_T1 == 1, "1", "0") )
data$G5<-factor( ifelse(data$G5 != 0 & data$CH_T2 == 1, "1", "0") )
data$G6<-factor( ifelse(data$G6 != 0 & data$CH_T3 == 1, "1", "0") )
data$G7<-factor( ifelse(data$G7 != 0 & data$CH_T1 == 1, "1", "0") )
data$G8<-factor( ifelse(data$G8 != 0 & data$CH_T2 == 1, "1", "0") )
data$G9<-factor( ifelse(data$G9 != 0 & data$CH_T3 == 1, "1", "0") )


#### PRELIMINARY DATA ANALYSIS  ####

#set path as working directory
setwd("D:/GRAZIA/PHD/TURKEY/R")

# load data #
load("D:/GRAZIA/PHD/TURKEY/R/data.RData")


##### Stratified Sampling for the Mann-Whitney-Wilcoxon Test and Spineplots  ######

# Initialize empty lists for storing stratified samples for each group (G0 to G9)
sample_G0<-NULL
# Loop 10 times to perform stratified sampling on the data
# Perform stratified sampling on the 'G0' variable, selecting 50% of rows where G0 == 1
# Store the result in the list 'sample_G0' at index 'i'
for (i in 1:10){
  sample_G0[[i]]<-stratified(data, "G0",.5*nrow(data[which(data$G0==1),]))
}

# Combine all 10 stratified samples into a single data frame
sample_G0<-do.call("rbind", sample_G0)

# Repeat the same process for other groups (G1 to G9)
# This process ensures that we create balanced stratified samples for each group
sample_G1<-NULL
for (i in 1:10){
  sample_G1[[i]]<-stratified(data, "G1",.5*nrow(data[which(data$G1==1),]))
}
sample_G1<-do.call("rbind", sample_G1)

sample_G2<-NULL
for (i in 1:10){
  sample_G2[[i]]<-stratified(data, "G2",.5*nrow(data[which(data$G2==1),]))
}
sample_G2<-do.call("rbind", sample_G2)

sample_G3<-NULL
for (i in 1:10){
  sample_G3[[i]]<-stratified(data, "G3",.5*nrow(data[which(data$G3==1),]))
}
sample_G3<-do.call("rbind", sample_G3)

sample_G4<-NULL
for (i in 1:10){
  sample_G4[[i]]<-stratified(data, "G4",.5*nrow(data[which(data$G4==1),]))
}
sample_G4<-do.call("rbind", sample_G4)

sample_G5<-NULL
for (i in 1:10){
  sample_G5[[i]]<-stratified(data, "G5",.5*nrow(data[which(data$G5==1),]))
}
sample_G5<-do.call("rbind", sample_G5)

sample_G6<-NULL
for (i in 1:10){
  sample_G6[[i]]<-stratified(data, "G6",.5*nrow(data[which(data$G6==1),]))
}
sample_G6<-do.call("rbind", sample_G6)

sample_G7<-NULL
for (i in 1:10){
  sample_G7[[i]]<-stratified(data, "G7",.5*nrow(data[which(data$G7==1),]))
}
sample_G7<-do.call("rbind", sample_G7)

sample_G8<-NULL
for (i in 1:10){
  sample_G8[[i]]<-stratified(data, "G8",.5*nrow(data[which(data$G8==1),]))
}
sample_G8<-do.call("rbind", sample_G8)

sample_G9<-NULL
for (i in 1:10){
  sample_G9[[i]]<-stratified(data, "G9",.5*nrow(data[which(data$G9==1),]))
}
sample_G9<-do.call("rbind", sample_G9)


##### Mann-Whitney-Wilcoxon Test #####
# Split the stratified samples into 'presence' and 'absence' groups for each variable
# 'presence' represents rows where the group variable equals 1
# 'absence' represents rows where the group variable equals 0
presence_G0<-sample_G0[which(sample_G0$G0==1),]
presence_G1<-sample_G1[which(sample_G1$G1==1),]
presence_G2<-sample_G2[which(sample_G2$G2==1),]
presence_G3<-sample_G3[which(sample_G3$G3==1),]
presence_G4<-sample_G4[which(sample_G4$G4==1),]
presence_G5<-sample_G5[which(sample_G5$G5==1),]
presence_G6<-sample_G6[which(sample_G6$G6==1),]
presence_G7<-sample_G7[which(sample_G7$G7==1),]
presence_G8<-sample_G8[which(sample_G8$G8==1),]
presence_G9<-sample_G9[which(sample_G9$G9==1),]

absence_G0<-sample_G0[which(sample_G0$G0==0),]
absence_G1<-sample_G1[which(sample_G1$G1==0),]
absence_G2<-sample_G2[which(sample_G2$G2==0),]
absence_G3<-sample_G3[which(sample_G3$G3==0),]
absence_G4<-sample_G4[which(sample_G4$G4==0),]
absence_G5<-sample_G5[which(sample_G5$G5==0),]
absence_G6<-sample_G6[which(sample_G6$G6==0),]
absence_G7<-sample_G7[which(sample_G7$G7==0),]
absence_G8<-sample_G8[which(sample_G8$G8==0),]
absence_G9<-sample_G9[which(sample_G9$G9==0),]

# Perform Mann-Whitney-Wilcoxon tests comparing the presence and absence groups for each variable
wilcox.test(presence_G0$SLO,absence_G0$SLO, paired=FALSE)
wilcox.test(presence_G0$PLC,absence_G0$PLC, paired=FALSE)
wilcox.test(presence_G0$DEV,absence_G0$DEV, paired=FALSE)
wilcox.test(presence_G0$EP,absence_G0$EP, paired=FALSE)
wilcox.test(presence_G0$SPI,absence_G0$SPI, paired=FALSE)
wilcox.test(presence_G0$TWI,absence_G0$TWI, paired=FALSE)
wilcox.test(presence_G0$LSF,absence_G0$LSF, paired=FALSE)

wilcox.test(presence_G1$SLO,absence_G1$SLO, paired=FALSE)
wilcox.test(presence_G1$PLC,absence_G1$PLC, paired=FALSE)
wilcox.test(presence_G1$DEV,absence_G1$DEV, paired=FALSE)
wilcox.test(presence_G1$EP,absence_G1$EP, paired=FALSE)
wilcox.test(presence_G1$SPI,absence_G1$SPI, paired=FALSE)
wilcox.test(presence_G1$TWI,absence_G1$TWI, paired=FALSE)
wilcox.test(presence_G1$LSF,absence_G1$LSF, paired=FALSE)

wilcox.test(presence_G2$SLO,absence_G2$SLO, paired=FALSE)
wilcox.test(presence_G2$PLC,absence_G2$PLC, paired=FALSE)
wilcox.test(presence_G2$DEV,absence_G2$DEV, paired=FALSE)
wilcox.test(presence_G2$EP,absence_G2$EP, paired=FALSE)
wilcox.test(presence_G2$SPI,absence_G2$SPI, paired=FALSE)
wilcox.test(presence_G2$TWI,absence_G2$TWI, paired=FALSE)
wilcox.test(presence_G2$LSF,absence_G2$LSF, paired=FALSE)

wilcox.test(presence_G3$SLO,absence_G3$SLO, paired=FALSE)
wilcox.test(presence_G3$PLC,absence_G3$PLC, paired=FALSE)
wilcox.test(presence_G3$DEV,absence_G3$DEV, paired=FALSE)
wilcox.test(presence_G3$EP,absence_G3$EP, paired=FALSE)
wilcox.test(presence_G3$SPI,absence_G3$SPI, paired=FALSE)
wilcox.test(presence_G3$TWI,absence_G3$TWI, paired=FALSE)
wilcox.test(presence_G3$LSF,absence_G3$LSF, paired=FALSE)

wilcox.test(presence_G4$SLO,absence_G4$SLO, paired=FALSE)
wilcox.test(presence_G4$PLC,absence_G4$PLC, paired=FALSE)
wilcox.test(presence_G4$DEV,absence_G4$DEV, paired=FALSE)
wilcox.test(presence_G4$EP,absence_G4$EP, paired=FALSE)
wilcox.test(presence_G4$SPI,absence_G4$SPI, paired=FALSE)
wilcox.test(presence_G4$TWI,absence_G4$TWI, paired=FALSE)
wilcox.test(presence_G4$LSF,absence_G4$LSF, paired=FALSE)

wilcox.test(presence_G5$SLO,absence_G5$SLO, paired=FALSE)
wilcox.test(presence_G5$PLC,absence_G5$PLC, paired=FALSE)
wilcox.test(presence_G5$DEV,absence_G5$DEV, paired=FALSE)
wilcox.test(presence_G5$EP,absence_G5$EP, paired=FALSE)
wilcox.test(presence_G5$SPI,absence_G5$SPI, paired=FALSE)
wilcox.test(presence_G5$TWI,absence_G5$TWI, paired=FALSE)
wilcox.test(presence_G5$LSF,absence_G5$LSF, paired=FALSE)

wilcox.test(presence_G6$SLO,absence_G6$SLO, paired=FALSE)
wilcox.test(presence_G6$PLC,absence_G6$PLC, paired=FALSE)
wilcox.test(presence_G6$DEV,absence_G6$DEV, paired=FALSE)
wilcox.test(presence_G6$EP,absence_G6$EP, paired=FALSE)
wilcox.test(presence_G6$SPI,absence_G6$SPI, paired=FALSE)
wilcox.test(presence_G6$TWI,absence_G6$TWI, paired=FALSE)
wilcox.test(presence_G6$LSF,absence_G6$LSF, paired=FALSE)

wilcox.test(presence_G7$SLO,absence_G7$SLO, paired=FALSE)
wilcox.test(presence_G7$PLC,absence_G7$PLC, paired=FALSE)
wilcox.test(presence_G7$DEV,absence_G7$DEV, paired=FALSE)
wilcox.test(presence_G7$EP,absence_G7$EP, paired=FALSE)
wilcox.test(presence_G7$SPI,absence_G7$SPI, paired=FALSE)
wilcox.test(presence_G7$TWI,absence_G7$TWI, paired=FALSE)
wilcox.test(presence_G7$LSF,absence_G7$LSF, paired=FALSE)

wilcox.test(presence_G8$SLO,absence_G8$SLO, paired=FALSE)
wilcox.test(presence_G8$PLC,absence_G8$PLC, paired=FALSE)
wilcox.test(presence_G8$DEV,absence_G8$DEV, paired=FALSE)
wilcox.test(presence_G8$EP,absence_G8$EP, paired=FALSE)
wilcox.test(presence_G8$SPI,absence_G8$SPI, paired=FALSE)
wilcox.test(presence_G8$TWI,absence_G8$TWI, paired=FALSE)
wilcox.test(presence_G8$LSF,absence_G8$LSF, paired=FALSE)

wilcox.test(presence_G9$SLO,absence_G9$SLO, paired=FALSE)
wilcox.test(presence_G9$PLC,absence_G9$PLC, paired=FALSE)
wilcox.test(presence_G9$DEV,absence_G9$DEV, paired=FALSE)
wilcox.test(presence_G9$EP,absence_G9$EP, paired=FALSE)
wilcox.test(presence_G9$SPI,absence_G9$SPI, paired=FALSE)
wilcox.test(presence_G9$TWI,absence_G9$TWI, paired=FALSE)
wilcox.test(presence_G9$LSF,absence_G9$LSF, paired=FALSE)

##### Spineplots #####
# Set the working directory for saving the figures
setwd("D:/GRAZIA/PHD/TURKEY/R/figures")

# Save spineplots as a PDF
pdf("spineplots.pdf",width=14, height=17)

# Adjust graphical parameters for multi-panel layout
par(mfrow=c(9,4))
par(mar = c(3.5, 2.5, 0.8, 2.5))
par(mgp=c(2,0.75,0))

# Generate spineplots for each variable and group
spineplot(G0 ~ SLO, sample_G0, col=c("#2c7fb8","#edf8b1"),xlab="", ylab="",xaxlabels =c("0%","25%","50%","75%","100%"),breaks = quantile(sample_G0$SLO))
text(0.06,0.9,substitute(paste(bold("SLO"))),cex=1)
spineplot(G1 ~ SLO, sample_G1, col=c("#2c7fb8","#edf8b1"),xlab="", ylab="",xaxlabels =c("0%","25%","50%","75%","100%"),breaks = quantile(sample_G1$SLO))
text(0.06,0.9,substitute(paste(bold("SLO"))),cex=1)
spineplot(G5 ~ SLO, sample_G5, col=c("#2c7fb8","#edf8b1"),xlab="", ylab="",xaxlabels =c("0%","25%","50%","75%","100%"),breaks = quantile(sample_G5$SLO))
text(0.06,0.9,substitute(paste(bold("SLO"))),cex=1)
spineplot(G9 ~ SLO, sample_G9, col=c("#2c7fb8","#edf8b1"),xlab="", ylab="",xaxlabels =c("0%","25%","50%","75%","100%"),breaks = quantile(sample_G9$SLO))
text(0.06,0.9,substitute(paste(bold("SLO"))),cex=1)


spineplot(G0 ~ PLC, sample_G0, col=c("#2c7fb8","#edf8b1"),xlab="", ylab="",xaxlabels =c("0%","25%","50%","75%","100%"),breaks = quantile(sample_G0$PLC))
text(0.06,0.9,substitute(paste(bold("PLC"))),cex=1)
spineplot(G1 ~ PLC, sample_G1, col=c("#2c7fb8","#edf8b1"),xlab="", ylab="",xaxlabels =c("0%","25%","50%","75%","100%"),breaks = quantile(sample_G1$PLC))
text(0.06,0.9,substitute(paste(bold("PLC"))),cex=1)
spineplot(G5 ~ PLC, sample_G5, col=c("#2c7fb8","#edf8b1"),xlab="", ylab="",xaxlabels =c("0%","25%","50%","75%","100%"),breaks = quantile(sample_G5$PLC))
text(0.06,0.9,substitute(paste(bold("PLC"))),cex=1)
spineplot(G9 ~ PLC, sample_G9, col=c("#2c7fb8","#edf8b1"),xlab="", ylab="",xaxlabels =c("0%","25%","50%","75%","100%"),breaks = quantile(sample_G9$PLC))
text(0.06,0.9,substitute(paste(bold("PLC"))),cex=1)


spineplot(G0 ~ DEV, sample_G0, col=c("#2c7fb8","#edf8b1"),xlab="", ylab="",xaxlabels =c("0%","25%","50%","75%","100%"),breaks = quantile(sample_G0$DEV))
text(0.06,0.9,substitute(paste(bold("DEV"))),cex=1)
spineplot(G1 ~ DEV, sample_G1, col=c("#2c7fb8","#edf8b1"),xlab="", ylab="",xaxlabels =c("0%","25%","50%","75%","100%"),breaks = quantile(sample_G1$DEV))
text(0.06,0.9,substitute(paste(bold("DEV"))),cex=1)
spineplot(G5 ~ DEV, sample_G5, col=c("#2c7fb8","#edf8b1"),xlab="", ylab="",xaxlabels =c("0%","25%","50%","75%","100%"),breaks = quantile(sample_G5$DEV))
text(0.06,0.9,substitute(paste(bold("DEV"))),cex=1)
spineplot(G9 ~ DEV, sample_G9, col=c("#2c7fb8","#edf8b1"),xlab="", ylab="",xaxlabels =c("0%","25%","50%","75%","100%"),breaks = quantile(sample_G9$DEV))
text(0.06,0.9,substitute(paste(bold("DEV"))),cex=1)


spineplot(G0 ~ EP, sample_G0, col=c("#2c7fb8","#edf8b1"),xlab="", ylab="",xaxlabels =c("0%","25%","50%","75%","100%"),breaks = quantile(sample_G0$EP))
text(0.06,0.95,substitute(paste(bold("EP"))),cex=1)
spineplot(G1 ~ EP, sample_G1, col=c("#2c7fb8","#edf8b1"),xlab="", ylab="",xaxlabels =c("0%","25%","50%","75%","100%"),breaks = quantile(sample_G1$EP))
text(0.06,0.95,substitute(paste(bold("EP"))),cex=1)
spineplot(G5 ~ EP, sample_G5, col=c("#2c7fb8","#edf8b1"),xlab="", ylab="",xaxlabels =c("0%","25%","50%","75%","100%"),breaks = quantile(sample_G5$EP))
text(0.06,0.95,substitute(paste(bold("EP"))),cex=1)
spineplot(G9 ~ EP, sample_G9, col=c("#2c7fb8","#edf8b1"),xlab="", ylab="",xaxlabels =c("0%","25%","50%","75%","100%"),breaks = quantile(sample_G9$EP))
text(0.06,0.95,substitute(paste(bold("EP"))),cex=1)


spineplot(G0 ~ GEO, sample_G0, col=c("#2c7fb8","#edf8b1"),xlab="", ylab="")
text(0.175,0.90,substitute(paste(bold("GEO"))),cex=1)
spineplot(G1 ~ GEO, sample_G1, col=c("#2c7fb8","#edf8b1"),xlab="", ylab="")
text(0.175,0.90,substitute(paste(bold("GEO"))),cex=1)
spineplot(G5 ~ GEO, sample_G5, col=c("#2c7fb8","#edf8b1"),xlab="", ylab="")
text(0.175,0.90,substitute(paste(bold("GEO"))),cex=1)
spineplot(G9 ~ GEO, sample_G9, col=c("#2c7fb8","#edf8b1"),xlab="", ylab="")
text(0.175,0.90,substitute(paste(bold("GEO"))),cex=1)


spineplot(G0 ~ GORD, sample_G0, col=c("#2c7fb8","#edf8b1"),xlab="", ylab="",
          breaks = c(1,2,3,4,5,6,10))
text(0.1,0.90,substitute(paste(bold("GORD"))),cex=1)
spineplot(G1 ~ GORD, sample_G1, col=c("#2c7fb8","#edf8b1"),xlab="", ylab="",
          breaks = c(1,2,3,4,5,6,10))
text(0.1,0.90,substitute(paste(bold("GORD"))),cex=1)
spineplot(G5 ~ GORD, sample_G5, col=c("#2c7fb8","#edf8b1"),xlab="", ylab="",
          breaks = c(1,2,3,4,5,6,10))
text(0.1,0.90,substitute(paste(bold("GORD"))),cex=1)
spineplot(G9 ~ GORD, sample_G9, col=c("#2c7fb8","#edf8b1"),xlab="", ylab="",
          breaks = c(1,2,3,4,5,6,10))
text(0.1,0.90,substitute(paste(bold("GORD"))),cex=1)


spineplot(G0 ~ SPI, sample_G0, col=c("#2c7fb8","#edf8b1"),xlab="", ylab="",xaxlabels =c("0%","25%","50%","75%","100%"),breaks = quantile(sample_G0$SPI))
text(0.06,0.90,substitute(paste(bold("SPI"))),cex=1)
spineplot(G1 ~ SPI, sample_G1, col=c("#2c7fb8","#edf8b1"),xlab="", ylab="",xaxlabels =c("0%","25%","50%","75%","100%"),breaks = quantile(sample_G1$SPI))
text(0.06,0.90,substitute(paste(bold("SPI"))),cex=1)
spineplot(G5 ~ SPI, sample_G5, col=c("#2c7fb8","#edf8b1"),xlab="", ylab="",xaxlabels =c("0%","25%","50%","75%","100%"),breaks = quantile(sample_G5$SPI))
text(0.06,0.90,substitute(paste(bold("SPI"))),cex=1)
spineplot(G9 ~ SPI, sample_G9, col=c("#2c7fb8","#edf8b1"),xlab="", ylab="",xaxlabels =c("0%","25%","50%","75%","100%"),breaks = quantile(sample_G9$SPI))
text(0.06,0.90,substitute(paste(bold("SPI"))),cex=1)


spineplot(G0 ~ TWI, sample_G0, col=c("#2c7fb8","#edf8b1"),xlab="", ylab="",xaxlabels =c("0%","25%","50%","75%","100%"),breaks = quantile(sample_G0$TWI))
text(0.06,0.90,substitute(paste(bold("TWI"))),cex=1)
spineplot(G1 ~ TWI, sample_G1, col=c("#2c7fb8","#edf8b1"),xlab="", ylab="",xaxlabels =c("0%","25%","50%","75%","100%"),breaks = quantile(sample_G1$TWI))
text(0.06,0.90,substitute(paste(bold("TWI"))),cex=1)
spineplot(G5 ~ TWI, sample_G5, col=c("#2c7fb8","#edf8b1"),xlab="", ylab="",xaxlabels =c("0%","25%","50%","75%","100%"),breaks = quantile(sample_G5$TWI))
text(0.06,0.90,substitute(paste(bold("TWI"))),cex=1)
spineplot(G9 ~ TWI, sample_G9, col=c("#2c7fb8","#edf8b1"),xlab="", ylab="",xaxlabels =c("0%","25%","50%","75%","100%"),breaks = quantile(sample_G9$TWI))
text(0.06,0.90,substitute(paste(bold("TWI"))),cex=1)

spineplot(G0 ~ LSF, sample_G0, col=c("#2c7fb8","#edf8b1"),xlab="", ylab="",xaxlabels =c("0%","25%","50%","75%","100%"),breaks = quantile(sample_G0$LSF))
text(0.06,0.9,substitute(paste(bold("LSF"))),cex=1)
spineplot(G1 ~ LSF, sample_G1, col=c("#2c7fb8","#edf8b1"),xlab="", ylab="",xaxlabels =c("0%","25%","50%","75%","100%"),breaks = quantile(sample_G1$LSF))
text(0.06,0.9,substitute(paste(bold("LSF"))),cex=1)
spineplot(G5 ~ LSF, sample_G5, col=c("#2c7fb8","#edf8b1"),xlab="", ylab="",xaxlabels =c("0%","25%","50%","75%","100%"),breaks = quantile(sample_G5$LSF))
text(0.06,0.9,substitute(paste(bold("LSF"))),cex=1)
spineplot(G9 ~ LSF, sample_G9, col=c("#2c7fb8","#edf8b1"),xlab="", ylab="",xaxlabels =c("0%","25%","50%","75%","100%"),breaks = quantile(sample_G9$LSF))
text(0.06,0.9,substitute(paste(bold("LSF"))),cex=1)

dev.off()

##### VIF calculation #####

data_vif <- data[,c(19,20,21,22,25:27)]  # Select columns 19, 20, 21, 22, and 25 to 27 as predictor variables.

target <- as.factor(data$G0)             # Convert the 'G0' column into a factor to create the target variable.

model <- glm(target ~ .,                 # Fit a logistic regression model with 'target' as the dependent variable
             data = data_vif,            # and selected columns ('data_vif') as predictors.
             family = binomial)          # Specify logistic regression using 'family = binomial'.

vif_values <- vif(model)                 # Calculate the Variance Inflation Factor (VIF) for multicollinearity assessment.
                                         # High VIF values (e.g., > 5 or > 10) indicate strong multicollinearity and may warrant further investigation.

print(vif_values)                        # Print the VIF values to identify highly collinear predictors.




#### MARS MODELLING ####
# Define the formulas for the models
formula_MAG0 <- G0 ~ ELE+SLO+GEO+PLC+DEV+EP
formula_MAG1 <- G1 ~ ELE+SLO+GEO+PLC+DEV+EP
formula_MAG2 <- G2 ~ ELE+SLO+GEO+PLC+DEV+EP
formula_MAG3 <- G3 ~ ELE+SLO+GEO+PLC+DEV+EP
formula_MAG4 <- G4 ~ ELE+SLO+GEO+PLC+DEV+EP
formula_MAG5 <- G5 ~ ELE+SLO+GEO+PLC+DEV+EP
formula_MAG6 <- G6 ~ ELE+SLO+GEO+PLC+DEV+EP
formula_MAG7 <- G7 ~ ELE+SLO+GEO+PLC+DEV+EP
formula_MAG8 <- G8 ~ ELE+SLO+GEO+PLC+DEV+EP
formula_MAG9 <- G9 ~ ELE+SLO+GEO+PLC+DEV+EP

formula_MBG0 <- G0 ~ ELE+SLO+GEO+PLC+DEV+EP+SPI+TWI+LSF+GORD
formula_MBG1 <- G1 ~ ELE+SLO+GEO+PLC+DEV+EP+SPI+TWI+LSF+GORD
formula_MBG2 <- G2 ~ ELE+SLO+GEO+PLC+DEV+EP+SPI+TWI+LSF+GORD
formula_MBG3 <- G3 ~ ELE+SLO+GEO+PLC+DEV+EP+SPI+TWI+LSF+GORD
formula_MBG4 <- G4 ~ ELE+SLO+GEO+PLC+DEV+EP+SPI+TWI+LSF+GORD
formula_MBG5 <- G5 ~ ELE+SLO+GEO+PLC+DEV+EP+SPI+TWI+LSF+GORD
formula_MBG6 <- G6 ~ ELE+SLO+GEO+PLC+DEV+EP+SPI+TWI+LSF+GORD
formula_MBG7 <- G7 ~ ELE+SLO+GEO+PLC+DEV+EP+SPI+TWI+LSF+GORD
formula_MBG8 <- G8 ~ ELE+SLO+GEO+PLC+DEV+EP+SPI+TWI+LSF+GORD
formula_MBG9 <- G9 ~ ELE+SLO+GEO+PLC+DEV+EP+SPI+TWI+LSF+GORD

##### Initialize Variables for the MAG0 Model #####

MAG0 <- NULL           # Initialize the MAG0 model variable as NULL.
training_set <- NULL   # Initialize the training dataset variable as NULL.
testing_set <- NULL    # Initialize the testing dataset variable as NULL.
dataf <- NULL          # Initialize the data frame variable as NULL.
auc_MAG0 <- NULL       # Initialize the variable to store AUC (Area Under the Curve) for MAG0 as NULL.
kappa_MAG0 <- NULL     # Initialize the variable to store Cohen's Kappa for MAG0 as NULL.

set.seed(123)          # Set the seed for random number generation to ensure reproducibility of results.


## Start of the 5-Fold Cross-Validation Loop ##

# Loop to perform cross-validation across 5 iterations.
for (i in 1:5) {
  
  # Generate 5 folds for cross-validation.
  # The `fold` function divides the dataset into 5 folds based on the 'tiles' column using a distance-based method.
  df_folded <- fold(
    data = data,
    k = 5,
    id_col = "tiles",
    method = "n_dist",
  )
  
  # Split data into individual folds.
  # Each fold is extracted based on the `.folds` column generated by the `fold` function.
  data_f1 <- df_folded[which(df_folded$.folds == "1"), ]
  data_f2 <- df_folded[which(df_folded$.folds == "2"), ]
  data_f3 <- df_folded[which(df_folded$.folds == "3"), ]
  data_f4 <- df_folded[which(df_folded$.folds == "4"), ]
  data_f5 <- df_folded[which(df_folded$.folds == "5"), ]

  # Stratify the data in each fold based on the target variable 'G0'.
  # Ensures balanced representation of class "1" samples in the dataset.
  sample_f1 <- stratified(data_f1, "G0",
                          nrow(data_f1[which(data_f1$G0 == "1"), ]))
  sample_f2 <- stratified(data_f2, "G0",
                          nrow(data_f2[which(data_f2$G0 == "1"), ]))
  sample_f3 <- stratified(data_f3, "G0",
                          nrow(data_f3[which(data_f3$G0 == "1"), ]))
  sample_f4 <- stratified(data_f4, "G0",
                          nrow(data_f4[which(data_f4$G0 == "1"), ]))
  sample_f5 <- stratified(data_f5, "G0",
                          nrow(data_f5[which(data_f5$G0 == "1"), ]))
  
  # Combine the stratified samples from all 5 folds into a single dataset for this iteration.
  dataf[[i]] <- rbind(sample_f1, sample_f2, sample_f3, sample_f4, sample_f5)

  # Initialize AUC and Kappa metrics matrices for the current iteration.
  auc_MAG0[[i]] <- matrix(nrow = 1, ncol = 5)
  kappa_MAG0[[i]] <- matrix(nrow = 1, ncol = 5)
  
  ## Loop Over Each Fold ##

  for (fold in 1:5) {
    # Split the dataset into training and testing sets.
    # Training data includes all folds except the current fold.
    training_set <- dataf[[i]][dataf[[i]]$.folds != fold, ]
    # Testing data is the current fold.
    testing_set <- dataf[[i]][dataf[[i]]$.folds == fold, ]

    # Train the MAG0 Model using the Earth Algorithm (MARS).
    # The Earth algorithm is a regression-based method. It is configured for logistic regression (binomial family).
    MAG0[[fold]] <- earth(formula_MAG0, data = training_set, degree = 1, trace = 1, glm = list(family = binomial))
    
    # Predict Scores on the Testing Set.
    # The trained model generates probabilities for the target class '1'.
    testing_set$score_MAG0 <- predict(MAG0[[fold]], testing_set, type = "response")
    
    # Calculate the Area Under the Curve (AUC).
    # AUC is a performance metric that evaluates the model's ability to discriminate between classes.
    auc_MAG0[[i]][, fold] <- as.vector(auc(testing_set$G0, testing_set$score_MAG0))
    
    # Convert Probabilities to Binary Class Predictions.
    # A threshold of 0.5 is used to classify samples as "0" or "1".
    testing_set$pred_MAG0 <- factor(ifelse(testing_set$score_MAG0 < 0.5, "0", "1"))
    
    # Generate the Confusion Matrix.
    # The confusion matrix is used to evaluate the agreement between predictions and true labels.
    conf_MAG0 <- confusionMatrix(testing_set$pred_MAG0, testing_set$G0, positive = "1")
    
    # Store the Kappa Statistic.
    # Kappa measures the agreement between predictions and true labels, adjusted for chance agreement.
    kappa_MAG0[[i]][, fold] <- conf_MAG0[["overall"]][["Kappa"]]
  }
}

##### Final Results for the MAG0 Model #####
# Combine AUC values from all iterations into a single vector for analysis.
auc_MAG0 <- unlist(auc_MAG0)
# Combine Kappa values from all iterations into a single vector for analysis.
kappa_MAG0 <- unlist(kappa_MAG0)


##### Initialize Variables for the MAG1 Model #####

MAG1<-NULL             
training_set<-NULL
testing_set<-NULL
dataf<-NULL
auc_MAG1<-NULL
kappa_MAG1<-NULL

set.seed(123)


for (i in 1:5){
  df_folded <- fold(
    data = data,
    k = 5,
    id_col = "tiles",
    method = "n_dist",
  )
  
  data_f1<-df_folded[which(df_folded$.folds=="1"),]
  data_f2<-df_folded[which(df_folded$.folds=="2"),]
  data_f3<-df_folded[which(df_folded$.folds=="3"),]
  data_f4<-df_folded[which(df_folded$.folds=="4"),]
  data_f5<-df_folded[which(df_folded$.folds=="5"),]
  
  sample_f1<-stratified(data_f1, "G1",
                        nrow(data_f1[which(data_f1$G1=="1"),]))
  sample_f2<-stratified(data_f2, "G1",
                        nrow(data_f2[which(data_f2$G1=="1"),]))
  sample_f3<-stratified(data_f3, "G1",
                        nrow(data_f3[which(data_f3$G1=="1"),]))
  sample_f4<-stratified(data_f4, "G1",
                        nrow(data_f4[which(data_f4$G1=="1"),]))
  sample_f5<-stratified(data_f5, "G1",
                        nrow(data_f5[which(data_f5$G1=="1"),]))
  
  dataf[[i]]<- rbind(sample_f1,sample_f2,sample_f3,sample_f4,sample_f5)
  
  auc_MAG1[[i]]<-matrix(nrow=1, ncol=5)
  kappa_MAG1[[i]]<-matrix(nrow=1, ncol=5)
  
  for (fold in 1:5){
    training_set <- dataf[[i]][dataf[[i]]$.folds != fold,]
    testing_set <- dataf[[i]][dataf[[i]]$.folds == fold,]
    
    MAG1[[fold]]<- earth(formula_MAG1, data=training_set, degree=1, trace=1, glm=list(family=binomial))
    
    testing_set$score_MAG1<-predict(MAG1[[fold]], testing_set, type="response")
    auc_MAG1[[i]][,fold] <- as.vector(auc(testing_set$G1,testing_set$score_MAG1))
    testing_set$pred_MAG1<-factor( ifelse(testing_set$score_MAG1 < 0.5, "0", "1") )
    conf_MAG1 <- confusionMatrix(testing_set$pred_MAG1,testing_set$G1, positive="1")
    kappa_MAG1[[i]][,fold]<-conf_MAG1[["overall"]][["Kappa"]]
  }
}

auc_MAG1<-unlist(auc_MAG1)
kappa_MAG1<-unlist(kappa_MAG1)

##### Initialize Variables for the MAG2 Model #####

MAG2<-NULL

training_set<-NULL
testing_set<-NULL
dataf<-NULL
auc_MAG2<-NULL
kappa_MAG2<-NULL

set.seed(123)


for (i in 1:5){
  df_folded <- fold(
    data = data,
    k = 5,
    id_col = "tiles",
    method = "n_dist",
  )
  
  data_f1<-df_folded[which(df_folded$.folds=="1"),]
  data_f2<-df_folded[which(df_folded$.folds=="2"),]
  data_f3<-df_folded[which(df_folded$.folds=="3"),]
  data_f4<-df_folded[which(df_folded$.folds=="4"),]
  data_f5<-df_folded[which(df_folded$.folds=="5"),]
  
  sample_f1<-stratified(data_f1, "G2",
                        nrow(data_f1[which(data_f1$G2=="1"),]))
  sample_f2<-stratified(data_f2, "G2",
                        nrow(data_f2[which(data_f2$G2=="1"),]))
  sample_f3<-stratified(data_f3, "G2",
                        nrow(data_f3[which(data_f3$G2=="1"),]))
  sample_f4<-stratified(data_f4, "G2",
                        nrow(data_f4[which(data_f4$G2=="1"),]))
  sample_f5<-stratified(data_f5, "G2",
                        nrow(data_f5[which(data_f5$G2=="1"),]))
  
  dataf[[i]]<- rbind(sample_f1,sample_f2,sample_f3,sample_f4,sample_f5)
  
  auc_MAG2[[i]]<-matrix(nrow=1, ncol=5)
  kappa_MAG2[[i]]<-matrix(nrow=1, ncol=5)
  
  for (fold in 1:5){
    training_set <- dataf[[i]][dataf[[i]]$.folds != fold,]
    testing_set <- dataf[[i]][dataf[[i]]$.folds == fold,]
    
    MAG2[[fold]]<- earth(formula_MAG2, data=training_set, degree=1, trace=1, glm=list(family=binomial))
    
    testing_set$score_MAG2<-predict(MAG2[[fold]], testing_set, type="response")
    auc_MAG2[[i]][,fold] <- as.vector(auc(testing_set$G2,testing_set$score_MAG2))
    testing_set$pred_MAG2<-factor( ifelse(testing_set$score_MAG2 < 0.5, "0", "1") )
    conf_MAG2 <- confusionMatrix(testing_set$pred_MAG2,testing_set$G2, positive="1")
    kappa_MAG2[[i]][,fold]<-conf_MAG2[["overall"]][["Kappa"]]
  }
}

auc_MAG2<-unlist(auc_MAG2)
kappa_MAG2<-unlist(kappa_MAG2)

##### Initialize Variables for the MAG3 Model #####

MAG3<-NULL

training_set<-NULL
testing_set<-NULL
dataf<-NULL
auc_MAG3<-NULL
kappa_MAG3<-NULL

set.seed(123)



for (i in 1:5){
  df_folded <- fold(
    data = data,
    k = 5,
    id_col = "tiles",
    method = "n_dist",
  )
  
  data_f1<-df_folded[which(df_folded$.folds=="1"),]
  data_f2<-df_folded[which(df_folded$.folds=="2"),]
  data_f3<-df_folded[which(df_folded$.folds=="3"),]
  data_f4<-df_folded[which(df_folded$.folds=="4"),]
  data_f5<-df_folded[which(df_folded$.folds=="5"),]
  
  sample_f1<-stratified(data_f1, "G3",
                        nrow(data_f1[which(data_f1$G3=="1"),]))
  sample_f2<-stratified(data_f2, "G3",
                        nrow(data_f2[which(data_f2$G3=="1"),]))
  sample_f3<-stratified(data_f3, "G3",
                        nrow(data_f3[which(data_f3$G3=="1"),]))
  sample_f4<-stratified(data_f4, "G3",
                        nrow(data_f4[which(data_f4$G3=="1"),]))
  sample_f5<-stratified(data_f5, "G3",
                        nrow(data_f5[which(data_f5$G3=="1"),]))
  
  dataf[[i]]<- rbind(sample_f1,sample_f2,sample_f3,sample_f4,sample_f5)
  
  auc_MAG3[[i]]<-matrix(nrow=1, ncol=5)
  kappa_MAG3[[i]]<-matrix(nrow=1, ncol=5)
  
  for (fold in 1:5){
    training_set <- dataf[[i]][dataf[[i]]$.folds != fold,]
    testing_set <- dataf[[i]][dataf[[i]]$.folds == fold,]
    
    MAG3[[fold]]<- earth(formula_MAG3, data=training_set, degree=1, trace=1, glm=list(family=binomial))
    
    testing_set$score_MAG3<-predict(MAG3[[fold]], testing_set, type="response")
    auc_MAG3[[i]][,fold] <- as.vector(auc(testing_set$G3,testing_set$score_MAG3))
    testing_set$pred_MAG3<-factor( ifelse(testing_set$score_MAG3 < 0.5, "0", "1") )
    conf_MAG3 <- confusionMatrix(testing_set$pred_MAG3,testing_set$G3, positive="1")
    kappa_MAG3[[i]][,fold]<-conf_MAG3[["overall"]][["Kappa"]]
  }
}

auc_MAG3<-unlist(auc_MAG3)
kappa_MAG3<-unlist(kappa_MAG3)

##### Initialize Variables for the MAG4 Model #####

MAG4<-NULL

training_set<-NULL
testing_set<-NULL
dataf<-NULL
auc_MAG4<-NULL
kappa_MAG4<-NULL

set.seed(123)



for (i in 1:5){
  df_folded <- fold(
    data = data,
    k = 5,
    id_col = "tiles",
    method = "n_dist",
  )
  
  data_f1<-df_folded[which(df_folded$.folds=="1"),]
  data_f2<-df_folded[which(df_folded$.folds=="2"),]
  data_f3<-df_folded[which(df_folded$.folds=="3"),]
  data_f4<-df_folded[which(df_folded$.folds=="4"),]
  data_f5<-df_folded[which(df_folded$.folds=="5"),]
  
  sample_f1<-stratified(data_f1, "G4",
                        nrow(data_f1[which(data_f1$G4=="1"),]))
  sample_f2<-stratified(data_f2, "G4",
                        nrow(data_f2[which(data_f2$G4=="1"),]))
  sample_f3<-stratified(data_f3, "G4",
                        nrow(data_f3[which(data_f3$G4=="1"),]))
  sample_f4<-stratified(data_f4, "G4",
                        nrow(data_f4[which(data_f4$G4=="1"),]))
  sample_f5<-stratified(data_f5, "G4",
                        nrow(data_f5[which(data_f5$G4=="1"),]))
  
  dataf[[i]]<- rbind(sample_f1,sample_f2,sample_f3,sample_f4,sample_f5)
  
  auc_MAG4[[i]]<-matrix(nrow=1, ncol=5)
  kappa_MAG4[[i]]<-matrix(nrow=1, ncol=5)
  
  for (fold in 1:5){
    training_set <- dataf[[i]][dataf[[i]]$.folds != fold,]
    testing_set <- dataf[[i]][dataf[[i]]$.folds == fold,]
    
    MAG4[[fold]]<- earth(formula_MAG4, data=training_set, degree=1, trace=1, glm=list(family=binomial))
    
    testing_set$score_MAG4<-predict(MAG4[[fold]], testing_set, type="response")
    auc_MAG4[[i]][,fold] <- as.vector(auc(testing_set$G4,testing_set$score_MAG4))
    testing_set$pred_MAG4<-factor( ifelse(testing_set$score_MAG4 < 0.5, "0", "1") )
    conf_MAG4 <- confusionMatrix(testing_set$pred_MAG4,testing_set$G4, positive="1")
    kappa_MAG4[[i]][,fold]<-conf_MAG4[["overall"]][["Kappa"]]
  }
}

auc_MAG4<-unlist(auc_MAG4)
kappa_MAG4<-unlist(kappa_MAG4)

##### Initialize Variables for the MAG5 Model #####

MAG5<-NULL

training_set<-NULL
testing_set<-NULL
dataf<-NULL
auc_MAG5<-NULL
kappa_MAG5<-NULL

set.seed(123)


for (i in 1:5){
  df_folded <- fold(
    data = data,
    k = 5,
    id_col = "tiles",
    method = "n_dist",
  )
  
  data_f1<-df_folded[which(df_folded$.folds=="1"),]
  data_f2<-df_folded[which(df_folded$.folds=="2"),]
  data_f3<-df_folded[which(df_folded$.folds=="3"),]
  data_f4<-df_folded[which(df_folded$.folds=="4"),]
  data_f5<-df_folded[which(df_folded$.folds=="5"),]
  
  sample_f1<-stratified(data_f1, "G5",
                        nrow(data_f1[which(data_f1$G5=="1"),]))
  sample_f2<-stratified(data_f2, "G5",
                        nrow(data_f2[which(data_f2$G5=="1"),]))
  sample_f3<-stratified(data_f3, "G5",
                        nrow(data_f3[which(data_f3$G5=="1"),]))
  sample_f4<-stratified(data_f4, "G5",
                        nrow(data_f4[which(data_f4$G5=="1"),]))
  sample_f5<-stratified(data_f5, "G5",
                        nrow(data_f5[which(data_f5$G5=="1"),]))
  
  dataf[[i]]<- rbind(sample_f1,sample_f2,sample_f3,sample_f4,sample_f5)
  
  auc_MAG5[[i]]<-matrix(nrow=1, ncol=5)
  kappa_MAG5[[i]]<-matrix(nrow=1, ncol=5)
  
  for (fold in 1:5){
    training_set <- dataf[[i]][dataf[[i]]$.folds != fold,]
    testing_set <- dataf[[i]][dataf[[i]]$.folds == fold,]
    
    MAG5[[fold]]<- earth(formula_MAG5, data=training_set, degree=1, trace=1, glm=list(family=binomial))
    
    testing_set$score_MAG5<-predict(MAG5[[fold]], testing_set, type="response")
    auc_MAG5[[i]][,fold] <- as.vector(auc(testing_set$G5,testing_set$score_MAG5))
    testing_set$pred_MAG5<-factor( ifelse(testing_set$score_MAG5 < 0.5, "0", "1") )
    conf_MAG5 <- confusionMatrix(testing_set$pred_MAG5,testing_set$G5, positive="1")
    kappa_MAG5[[i]][,fold]<-conf_MAG5[["overall"]][["Kappa"]]
  }
}

auc_MAG5<-unlist(auc_MAG5)
kappa_MAG5<-unlist(kappa_MAG5)

##### Initialize Variables for the MAG6 Model #####

MAG6<-NULL

training_set<-NULL
testing_set<-NULL
dataf<-NULL
auc_MAG6<-NULL
kappa_MAG6<-NULL

set.seed(123)


for (i in 1:5){
  df_folded <- fold(
    data = data,
    k = 5,
    id_col = "tiles",
    method = "n_dist",
  )
  
  data_f1<-df_folded[which(df_folded$.folds=="1"),]
  data_f2<-df_folded[which(df_folded$.folds=="2"),]
  data_f3<-df_folded[which(df_folded$.folds=="3"),]
  data_f4<-df_folded[which(df_folded$.folds=="4"),]
  data_f5<-df_folded[which(df_folded$.folds=="5"),]
  
  sample_f1<-stratified(data_f1, "G6",
                        nrow(data_f1[which(data_f1$G6=="1"),]))
  sample_f2<-stratified(data_f2, "G6",
                        nrow(data_f2[which(data_f2$G6=="1"),]))
  sample_f3<-stratified(data_f3, "G6",
                        nrow(data_f3[which(data_f3$G6=="1"),]))
  sample_f4<-stratified(data_f4, "G6",
                        nrow(data_f4[which(data_f4$G6=="1"),]))
  sample_f5<-stratified(data_f5, "G6",
                        nrow(data_f5[which(data_f5$G6=="1"),]))
  
  dataf[[i]]<- rbind(sample_f1,sample_f2,sample_f3,sample_f4,sample_f5)
  
  auc_MAG6[[i]]<-matrix(nrow=1, ncol=5)
  kappa_MAG6[[i]]<-matrix(nrow=1, ncol=5)
  
  for (fold in 1:5){
    training_set <- dataf[[i]][dataf[[i]]$.folds != fold,]
    testing_set <- dataf[[i]][dataf[[i]]$.folds == fold,]
    
    MAG6[[fold]]<- earth(formula_MAG6, data=training_set, degree=1, trace=1, glm=list(family=binomial))
    
    testing_set$score_MAG6<-predict(MAG6[[fold]], testing_set, type="response")
    auc_MAG6[[i]][,fold] <- as.vector(auc(testing_set$G6,testing_set$score_MAG6))
    testing_set$pred_MAG6<-factor( ifelse(testing_set$score_MAG6 < 0.5, "0", "1") )
    conf_MAG6 <- confusionMatrix(testing_set$pred_MAG6,testing_set$G6, positive="1")
    kappa_MAG6[[i]][,fold]<-conf_MAG6[["overall"]][["Kappa"]]
  }
}

auc_MAG6<-unlist(auc_MAG6)
kappa_MAG6<-unlist(kappa_MAG6)

##### Initialize Variables for the MAG7 Model #####

MAG7<-NULL

training_set<-NULL
testing_set<-NULL
dataf<-NULL
auc_MAG7<-NULL
kappa_MAG7<-NULL

set.seed(123)


for (i in 1:5){
  df_folded <- fold(
    data = data,
    k = 5,
    id_col = "tiles",
    method = "n_dist",
  )
  
  data_f1<-df_folded[which(df_folded$.folds=="1"),]
  data_f2<-df_folded[which(df_folded$.folds=="2"),]
  data_f3<-df_folded[which(df_folded$.folds=="3"),]
  data_f4<-df_folded[which(df_folded$.folds=="4"),]
  data_f5<-df_folded[which(df_folded$.folds=="5"),]
  
  sample_f1<-stratified(data_f1, "G7",
                        nrow(data_f1[which(data_f1$G7=="1"),]))
  sample_f2<-stratified(data_f2, "G7",
                        nrow(data_f2[which(data_f2$G7=="1"),]))
  sample_f3<-stratified(data_f3, "G7",
                        nrow(data_f3[which(data_f3$G7=="1"),]))
  sample_f4<-stratified(data_f4, "G7",
                        nrow(data_f4[which(data_f4$G7=="1"),]))
  sample_f5<-stratified(data_f5, "G7",
                        nrow(data_f5[which(data_f5$G7=="1"),]))
  
  dataf[[i]]<- rbind(sample_f1,sample_f2,sample_f3,sample_f4,sample_f5)
  
  auc_MAG7[[i]]<-matrix(nrow=1, ncol=5)
  kappa_MAG7[[i]]<-matrix(nrow=1, ncol=5)
  
  for (fold in 1:5){
    training_set <- dataf[[i]][dataf[[i]]$.folds != fold,]
    testing_set <- dataf[[i]][dataf[[i]]$.folds == fold,]
    
    MAG7[[fold]]<- earth(formula_MAG7, data=training_set, degree=1, trace=1, glm=list(family=binomial))
    
    testing_set$score_MAG7<-predict(MAG7[[fold]], testing_set, type="response")
    auc_MAG7[[i]][,fold] <- as.vector(auc(testing_set$G7,testing_set$score_MAG7))
    testing_set$pred_MAG7<-factor( ifelse(testing_set$score_MAG7 < 0.5, "0", "1") )
    conf_MAG7 <- confusionMatrix(testing_set$pred_MAG7,testing_set$G7, positive="1")
    kappa_MAG7[[i]][,fold]<-conf_MAG7[["overall"]][["Kappa"]]
  }
}

auc_MAG7<-unlist(auc_MAG7)
kappa_MAG7<-unlist(kappa_MAG7)

##### Initialize Variables for the MAG8 Model #####

MAG8<-NULL

training_set<-NULL
testing_set<-NULL
dataf<-NULL
auc_MAG8<-NULL
kappa_MAG8<-NULL

set.seed(123)



for (i in 1:5){
  df_folded <- fold(
    data = data,
    k = 5,
    id_col = "tiles",
    method = "n_dist",
  )
  
  data_f1<-df_folded[which(df_folded$.folds=="1"),]
  data_f2<-df_folded[which(df_folded$.folds=="2"),]
  data_f3<-df_folded[which(df_folded$.folds=="3"),]
  data_f4<-df_folded[which(df_folded$.folds=="4"),]
  data_f5<-df_folded[which(df_folded$.folds=="5"),]
  
  sample_f1<-stratified(data_f1, "G8",
                        nrow(data_f1[which(data_f1$G8=="1"),]))
  sample_f2<-stratified(data_f2, "G8",
                        nrow(data_f2[which(data_f2$G8=="1"),]))
  sample_f3<-stratified(data_f3, "G8",
                        nrow(data_f3[which(data_f3$G8=="1"),]))
  sample_f4<-stratified(data_f4, "G8",
                        nrow(data_f4[which(data_f4$G8=="1"),]))
  sample_f5<-stratified(data_f5, "G8",
                        nrow(data_f5[which(data_f5$G8=="1"),]))
  
  dataf[[i]]<- rbind(sample_f1,sample_f2,sample_f3,sample_f4,sample_f5)
  
  auc_MAG8[[i]]<-matrix(nrow=1, ncol=5)
  kappa_MAG8[[i]]<-matrix(nrow=1, ncol=5)
  
  for (fold in 1:5){
    training_set <- dataf[[i]][dataf[[i]]$.folds != fold,]
    testing_set <- dataf[[i]][dataf[[i]]$.folds == fold,]
    
    MAG8[[fold]]<- earth(formula_MAG8, data=training_set, degree=1, trace=1, glm=list(family=binomial))
    
    testing_set$score_MAG8<-predict(MAG8[[fold]], testing_set, type="response")
    auc_MAG8[[i]][,fold] <- as.vector(auc(testing_set$G8,testing_set$score_MAG8))
    testing_set$pred_MAG8<-factor( ifelse(testing_set$score_MAG8 < 0.5, "0", "1") )
    conf_MAG8 <- confusionMatrix(testing_set$pred_MAG8,testing_set$G8, positive="1")
    kappa_MAG8[[i]][,fold]<-conf_MAG8[["overall"]][["Kappa"]]
  }
}

auc_MAG8<-unlist(auc_MAG8)
kappa_MAG8<-unlist(kappa_MAG8)

##### Initialize Variables for the MAG9 Model #####

MAG9<-NULL

training_set<-NULL
testing_set<-NULL
dataf<-NULL
auc_MAG9<-NULL
kappa_MAG9<-NULL

set.seed(123)


for (i in 1:5){
  df_folded <- fold(
    data = data,
    k = 5,
    id_col = "tiles",
    method = "n_dist",
  )
  
  data_f1<-df_folded[which(df_folded$.folds=="1"),]
  data_f2<-df_folded[which(df_folded$.folds=="2"),]
  data_f3<-df_folded[which(df_folded$.folds=="3"),]
  data_f4<-df_folded[which(df_folded$.folds=="4"),]
  data_f5<-df_folded[which(df_folded$.folds=="5"),]
  
  sample_f1<-stratified(data_f1, "G9",
                        nrow(data_f1[which(data_f1$G9=="1"),]))
  sample_f2<-stratified(data_f2, "G9",
                        nrow(data_f2[which(data_f2$G9=="1"),]))
  sample_f3<-stratified(data_f3, "G9",
                        nrow(data_f3[which(data_f3$G9=="1"),]))
  sample_f4<-stratified(data_f4, "G9",
                        nrow(data_f4[which(data_f4$G9=="1"),]))
  sample_f5<-stratified(data_f5, "G9",
                        nrow(data_f5[which(data_f5$G9=="1"),]))
  
  dataf[[i]]<- rbind(sample_f1,sample_f2,sample_f3,sample_f4,sample_f5)
  
  auc_MAG9[[i]]<-matrix(nrow=1, ncol=5)
  kappa_MAG9[[i]]<-matrix(nrow=1, ncol=5)
  
  for (fold in 1:5){
    training_set <- dataf[[i]][dataf[[i]]$.folds != fold,]
    testing_set <- dataf[[i]][dataf[[i]]$.folds == fold,]
    
    MAG9[[fold]]<- earth(formula_MAG9, data=training_set, degree=1, trace=1, glm=list(family=binomial))
    
    testing_set$score_MAG9<-predict(MAG9[[fold]], testing_set, type="response")
    auc_MAG9[[i]][,fold] <- as.vector(auc(testing_set$G9,testing_set$score_MAG9))
    testing_set$pred_MAG9<-factor( ifelse(testing_set$score_MAG9 < 0.5, "0", "1") )
    conf_MAG9 <- confusionMatrix(testing_set$pred_MAG9,testing_set$G9, positive="1")
    kappa_MAG9[[i]][,fold]<-conf_MAG9[["overall"]][["Kappa"]]
  }
}

auc_MAG9<-unlist(auc_MAG9)
kappa_MAG9<-unlist(kappa_MAG9)


##### Initialize Variables for the MBG0 Model #####

MBG0<-NULL

training_set<-NULL
testing_set<-NULL
dataf<-NULL
auc_MBG0<-NULL
kappa_MBG0<-NULL

set.seed(123)


for (i in 1:5){
  df_folded <- fold(
    data = data,
    k = 5,
    id_col = "tiles",
    method = "n_dist",
  )
  
  data_f1<-df_folded[which(df_folded$.folds=="1"),]
  data_f2<-df_folded[which(df_folded$.folds=="2"),]
  data_f3<-df_folded[which(df_folded$.folds=="3"),]
  data_f4<-df_folded[which(df_folded$.folds=="4"),]
  data_f5<-df_folded[which(df_folded$.folds=="5"),]
  
  sample_f1<-stratified(data_f1, "G0",
                        nrow(data_f1[which(data_f1$G0=="1"),]))
  sample_f2<-stratified(data_f2, "G0",
                        nrow(data_f2[which(data_f2$G0=="1"),]))
  sample_f3<-stratified(data_f3, "G0",
                        nrow(data_f3[which(data_f3$G0=="1"),]))
  sample_f4<-stratified(data_f4, "G0",
                        nrow(data_f4[which(data_f4$G0=="1"),]))
  sample_f5<-stratified(data_f5, "G0",
                        nrow(data_f5[which(data_f5$G0=="1"),]))
  
  dataf[[i]]<- rbind(sample_f1,sample_f2,sample_f3,sample_f4,sample_f5)
  
  auc_MBG0[[i]]<-matrix(nrow=1, ncol=5)
  kappa_MBG0[[i]]<-matrix(nrow=1, ncol=5)
  
  for (fold in 1:5){
    training_set <- dataf[[i]][dataf[[i]]$.folds != fold,]
    testing_set <- dataf[[i]][dataf[[i]]$.folds == fold,]
    
    MBG0[[fold]]<- earth(formula_MBG0, data=training_set, degree=1, trace=1, glm=list(family=binomial))
    
    testing_set$score_MBG0<-predict(MBG0[[fold]], testing_set, type="response")
    auc_MBG0[[i]][,fold] <- as.vector(auc(testing_set$G0,testing_set$score_MBG0))
    testing_set$pred_MBG0<-factor( ifelse(testing_set$score_MBG0 < 0.5, "0", "1") )
    conf_MBG0 <- confusionMatrix(testing_set$pred_MBG0,testing_set$G0, positive="1")
    kappa_MBG0[[i]][,fold]<-conf_MBG0[["overall"]][["Kappa"]]
  }
}

auc_MBG0<-unlist(auc_MBG0)
kappa_MBG0<-unlist(kappa_MBG0)

##### Initialize Variables for the MBG1 Model #####

MBG1<-NULL

training_set<-NULL
testing_set<-NULL
dataf<-NULL
auc_MBG1<-NULL
kappa_MBG1<-NULL

set.seed(123)


for (i in 1:5){
  df_folded <- fold(
    data = data,
    k = 5,
    id_col = "tiles",
    method = "n_dist",
  )
  
  data_f1<-df_folded[which(df_folded$.folds=="1"),]
  data_f2<-df_folded[which(df_folded$.folds=="2"),]
  data_f3<-df_folded[which(df_folded$.folds=="3"),]
  data_f4<-df_folded[which(df_folded$.folds=="4"),]
  data_f5<-df_folded[which(df_folded$.folds=="5"),]
  
  sample_f1<-stratified(data_f1, "G1",
                        nrow(data_f1[which(data_f1$G1=="1"),]))
  sample_f2<-stratified(data_f2, "G1",
                        nrow(data_f2[which(data_f2$G1=="1"),]))
  sample_f3<-stratified(data_f3, "G1",
                        nrow(data_f3[which(data_f3$G1=="1"),]))
  sample_f4<-stratified(data_f4, "G1",
                        nrow(data_f4[which(data_f4$G1=="1"),]))
  sample_f5<-stratified(data_f5, "G1",
                        nrow(data_f5[which(data_f5$G1=="1"),]))
  
  dataf[[i]]<- rbind(sample_f1,sample_f2,sample_f3,sample_f4,sample_f5)
  
  auc_MBG1[[i]]<-matrix(nrow=1, ncol=5)
  kappa_MBG1[[i]]<-matrix(nrow=1, ncol=5)
  
  for (fold in 1:5){
    training_set <- dataf[[i]][dataf[[i]]$.folds != fold,]
    testing_set <- dataf[[i]][dataf[[i]]$.folds == fold,]
    
    MBG1[[fold]]<- earth(formula_MBG1, data=training_set, degree=1, trace=1, glm=list(family=binomial))
    
    testing_set$score_MBG1<-predict(MBG1[[fold]], testing_set, type="response")
    auc_MBG1[[i]][,fold] <- as.vector(auc(testing_set$G1,testing_set$score_MBG1))
    testing_set$pred_MBG1<-factor( ifelse(testing_set$score_MBG1 < 0.5, "0", "1") )
    conf_MBG1 <- confusionMatrix(testing_set$pred_MBG1,testing_set$G1, positive="1")
    kappa_MBG1[[i]][,fold]<-conf_MBG1[["overall"]][["Kappa"]]
  }
}

auc_MBG1<-unlist(auc_MBG1)
kappa_MBG1<-unlist(kappa_MBG1)

##### Initialize Variables for the MBG2 Model #####

MBG2<-NULL

training_set<-NULL
testing_set<-NULL
dataf<-NULL
auc_MBG2<-NULL
kappa_MBG2<-NULL

set.seed(123)


for (i in 1:5){
  df_folded <- fold(
    data = data,
    k = 5,
    id_col = "tiles",
    method = "n_dist",
  )
  
  data_f1<-df_folded[which(df_folded$.folds=="1"),]
  data_f2<-df_folded[which(df_folded$.folds=="2"),]
  data_f3<-df_folded[which(df_folded$.folds=="3"),]
  data_f4<-df_folded[which(df_folded$.folds=="4"),]
  data_f5<-df_folded[which(df_folded$.folds=="5"),]
  
  sample_f1<-stratified(data_f1, "G2",
                        nrow(data_f1[which(data_f1$G2=="1"),]))
  sample_f2<-stratified(data_f2, "G2",
                        nrow(data_f2[which(data_f2$G2=="1"),]))
  sample_f3<-stratified(data_f3, "G2",
                        nrow(data_f3[which(data_f3$G2=="1"),]))
  sample_f4<-stratified(data_f4, "G2",
                        nrow(data_f4[which(data_f4$G2=="1"),]))
  sample_f5<-stratified(data_f5, "G2",
                        nrow(data_f5[which(data_f5$G2=="1"),]))
  
  dataf[[i]]<- rbind(sample_f1,sample_f2,sample_f3,sample_f4,sample_f5)
  
  auc_MBG2[[i]]<-matrix(nrow=1, ncol=5)
  kappa_MBG2[[i]]<-matrix(nrow=1, ncol=5)
  
  for (fold in 1:5){
    training_set <- dataf[[i]][dataf[[i]]$.folds != fold,]
    testing_set <- dataf[[i]][dataf[[i]]$.folds == fold,]
    
    MBG2[[fold]]<- earth(formula_MBG2, data=training_set, degree=1, trace=1, glm=list(family=binomial))
    
    testing_set$score_MBG2<-predict(MBG2[[fold]], testing_set, type="response")
    auc_MBG2[[i]][,fold] <- as.vector(auc(testing_set$G2,testing_set$score_MBG2))
    testing_set$pred_MBG2<-factor( ifelse(testing_set$score_MBG2 < 0.5, "0", "1") )
    conf_MBG2 <- confusionMatrix(testing_set$pred_MBG2,testing_set$G2, positive="1")
    kappa_MBG2[[i]][,fold]<-conf_MBG2[["overall"]][["Kappa"]]
  }
}

auc_MBG2<-unlist(auc_MBG2)
kappa_MBG2<-unlist(kappa_MBG2)

##### Initialize Variables for the MBG3 Model #####

MBG3<-NULL

training_set<-NULL
testing_set<-NULL
dataf<-NULL
auc_MBG3<-NULL
kappa_MBG3<-NULL

set.seed(123)


for (i in 1:5){
  df_folded <- fold(
    data = data,
    k = 5,
    id_col = "tiles",
    method = "n_dist",
  )
  
  data_f1<-df_folded[which(df_folded$.folds=="1"),]
  data_f2<-df_folded[which(df_folded$.folds=="2"),]
  data_f3<-df_folded[which(df_folded$.folds=="3"),]
  data_f4<-df_folded[which(df_folded$.folds=="4"),]
  data_f5<-df_folded[which(df_folded$.folds=="5"),]
  
  sample_f1<-stratified(data_f1, "G3",
                        nrow(data_f1[which(data_f1$G3=="1"),]))
  sample_f2<-stratified(data_f2, "G3",
                        nrow(data_f2[which(data_f2$G3=="1"),]))
  sample_f3<-stratified(data_f3, "G3",
                        nrow(data_f3[which(data_f3$G3=="1"),]))
  sample_f4<-stratified(data_f4, "G3",
                        nrow(data_f4[which(data_f4$G3=="1"),]))
  sample_f5<-stratified(data_f5, "G3",
                        nrow(data_f5[which(data_f5$G3=="1"),]))
  
  dataf[[i]]<- rbind(sample_f1,sample_f2,sample_f3,sample_f4,sample_f5)
  
  auc_MBG3[[i]]<-matrix(nrow=1, ncol=5)
  kappa_MBG3[[i]]<-matrix(nrow=1, ncol=5)
  
  for (fold in 1:5){
    training_set <- dataf[[i]][dataf[[i]]$.folds != fold,]
    testing_set <- dataf[[i]][dataf[[i]]$.folds == fold,]
    
    MBG3[[fold]]<- earth(formula_MBG3, data=training_set, degree=1, trace=1, glm=list(family=binomial))
    
    testing_set$score_MBG3<-predict(MBG3[[fold]], testing_set, type="response")
    auc_MBG3[[i]][,fold] <- as.vector(auc(testing_set$G3,testing_set$score_MBG3))
    testing_set$pred_MBG3<-factor( ifelse(testing_set$score_MBG3 < 0.5, "0", "1") )
    conf_MBG3 <- confusionMatrix(testing_set$pred_MBG3,testing_set$G3, positive="1")
    kappa_MBG3[[i]][,fold]<-conf_MBG3[["overall"]][["Kappa"]]
  }
}

auc_MBG3<-unlist(auc_MBG3)
kappa_MBG3<-unlist(kappa_MBG3)

##### Initialize Variables for the MBG4 Model #####

MBG4<-NULL

training_set<-NULL
testing_set<-NULL
dataf<-NULL
auc_MBG4<-NULL
kappa_MBG4<-NULL

set.seed(123)


for (i in 1:5){
  df_folded <- fold(
    data = data,
    k = 5,
    id_col = "tiles",
    method = "n_dist",
  )
  
  data_f1<-df_folded[which(df_folded$.folds=="1"),]
  data_f2<-df_folded[which(df_folded$.folds=="2"),]
  data_f3<-df_folded[which(df_folded$.folds=="3"),]
  data_f4<-df_folded[which(df_folded$.folds=="4"),]
  data_f5<-df_folded[which(df_folded$.folds=="5"),]
  
  sample_f1<-stratified(data_f1, "G4",
                        nrow(data_f1[which(data_f1$G4=="1"),]))
  sample_f2<-stratified(data_f2, "G4",
                        nrow(data_f2[which(data_f2$G4=="1"),]))
  sample_f3<-stratified(data_f3, "G4",
                        nrow(data_f3[which(data_f3$G4=="1"),]))
  sample_f4<-stratified(data_f4, "G4",
                        nrow(data_f4[which(data_f4$G4=="1"),]))
  sample_f5<-stratified(data_f5, "G4",
                        nrow(data_f5[which(data_f5$G4=="1"),]))
  
  dataf[[i]]<- rbind(sample_f1,sample_f2,sample_f3,sample_f4,sample_f5)
  
  auc_MBG4[[i]]<-matrix(nrow=1, ncol=5)
  kappa_MBG4[[i]]<-matrix(nrow=1, ncol=5)
  
  for (fold in 1:5){
    training_set <- dataf[[i]][dataf[[i]]$.folds != fold,]
    testing_set <- dataf[[i]][dataf[[i]]$.folds == fold,]
    
    MBG4[[fold]]<- earth(formula_MBG4, data=training_set, degree=1, trace=1, glm=list(family=binomial))
    
    testing_set$score_MBG4<-predict(MBG4[[fold]], testing_set, type="response")
    auc_MBG4[[i]][,fold] <- as.vector(auc(testing_set$G4,testing_set$score_MBG4))
    testing_set$pred_MBG4<-factor( ifelse(testing_set$score_MBG4 < 0.5, "0", "1") )
    conf_MBG4 <- confusionMatrix(testing_set$pred_MBG4,testing_set$G4, positive="1")
    kappa_MBG4[[i]][,fold]<-conf_MBG4[["overall"]][["Kappa"]]
  }
}

auc_MBG4<-unlist(auc_MBG4)
kappa_MBG4<-unlist(kappa_MBG4)

##### Initialize Variables for the MBG5 Model #####

MBG5<-NULL

training_set<-NULL
testing_set<-NULL
dataf<-NULL
auc_MBG5<-NULL
kappa_MBG5<-NULL

set.seed(123)


for (i in 1:5){
  df_folded <- fold(
    data = data,
    k = 5,
    id_col = "tiles",
    method = "n_dist",
  )
  
  data_f1<-df_folded[which(df_folded$.folds=="1"),]
  data_f2<-df_folded[which(df_folded$.folds=="2"),]
  data_f3<-df_folded[which(df_folded$.folds=="3"),]
  data_f4<-df_folded[which(df_folded$.folds=="4"),]
  data_f5<-df_folded[which(df_folded$.folds=="5"),]
  
  sample_f1<-stratified(data_f1, "G5",
                        nrow(data_f1[which(data_f1$G5=="1"),]))
  sample_f2<-stratified(data_f2, "G5",
                        nrow(data_f2[which(data_f2$G5=="1"),]))
  sample_f3<-stratified(data_f3, "G5",
                        nrow(data_f3[which(data_f3$G5=="1"),]))
  sample_f4<-stratified(data_f4, "G5",
                        nrow(data_f4[which(data_f4$G5=="1"),]))
  sample_f5<-stratified(data_f5, "G5",
                        nrow(data_f5[which(data_f5$G5=="1"),]))
  
  dataf[[i]]<- rbind(sample_f1,sample_f2,sample_f3,sample_f4,sample_f5)
  
  auc_MBG5[[i]]<-matrix(nrow=1, ncol=5)
  kappa_MBG5[[i]]<-matrix(nrow=1, ncol=5)
  
  for (fold in 1:5){
    training_set <- dataf[[i]][dataf[[i]]$.folds != fold,]
    testing_set <- dataf[[i]][dataf[[i]]$.folds == fold,]
    
    MBG5[[fold]]<- earth(formula_MBG5, data=training_set, degree=1, trace=1, glm=list(family=binomial))
    
    testing_set$score_MBG5<-predict(MBG5[[fold]], testing_set, type="response")
    auc_MBG5[[i]][,fold] <- as.vector(auc(testing_set$G5,testing_set$score_MBG5))
    testing_set$pred_MBG5<-factor( ifelse(testing_set$score_MBG5 < 0.5, "0", "1") )
    conf_MBG5 <- confusionMatrix(testing_set$pred_MBG5,testing_set$G5, positive="1")
    kappa_MBG5[[i]][,fold]<-conf_MBG5[["overall"]][["Kappa"]]
  }
}

auc_MBG5<-unlist(auc_MBG5)
kappa_MBG5<-unlist(kappa_MBG5)

##### Initialize Variables for the MBG6 Model #####

MBG6<-NULL

training_set<-NULL
testing_set<-NULL
dataf<-NULL
auc_MBG6<-NULL
kappa_MBG6<-NULL

set.seed(123)


for (i in 1:5){
  df_folded <- fold(
    data = data,
    k = 5,
    id_col = "tiles",
    method = "n_dist",
  )
  
  data_f1<-df_folded[which(df_folded$.folds=="1"),]
  data_f2<-df_folded[which(df_folded$.folds=="2"),]
  data_f3<-df_folded[which(df_folded$.folds=="3"),]
  data_f4<-df_folded[which(df_folded$.folds=="4"),]
  data_f5<-df_folded[which(df_folded$.folds=="5"),]
  
  sample_f1<-stratified(data_f1, "G6",
                        nrow(data_f1[which(data_f1$G6=="1"),]))
  sample_f2<-stratified(data_f2, "G6",
                        nrow(data_f2[which(data_f2$G6=="1"),]))
  sample_f3<-stratified(data_f3, "G6",
                        nrow(data_f3[which(data_f3$G6=="1"),]))
  sample_f4<-stratified(data_f4, "G6",
                        nrow(data_f4[which(data_f4$G6=="1"),]))
  sample_f5<-stratified(data_f5, "G6",
                        nrow(data_f5[which(data_f5$G6=="1"),]))
  
  dataf[[i]]<- rbind(sample_f1,sample_f2,sample_f3,sample_f4,sample_f5)
  
  auc_MBG6[[i]]<-matrix(nrow=1, ncol=5)
  kappa_MBG6[[i]]<-matrix(nrow=1, ncol=5)
  
  for (fold in 1:5){
    training_set <- dataf[[i]][dataf[[i]]$.folds != fold,]
    testing_set <- dataf[[i]][dataf[[i]]$.folds == fold,]
    
    MBG6[[fold]]<- earth(formula_MBG6, data=training_set, degree=1, trace=1, glm=list(family=binomial))
    
    testing_set$score_MBG6<-predict(MBG6[[fold]], testing_set, type="response")
    auc_MBG6[[i]][,fold] <- as.vector(auc(testing_set$G6,testing_set$score_MBG6))
    testing_set$pred_MBG6<-factor( ifelse(testing_set$score_MBG6 < 0.5, "0", "1") )
    conf_MBG6 <- confusionMatrix(testing_set$pred_MBG6,testing_set$G6, positive="1")
    kappa_MBG6[[i]][,fold]<-conf_MBG6[["overall"]][["Kappa"]]
  }
}

auc_MBG6<-unlist(auc_MBG6)
kappa_MBG6<-unlist(kappa_MBG6)

##### Initialize Variables for the MBG7 Model #####

MBG7<-NULL

training_set<-NULL
testing_set<-NULL
dataf<-NULL
auc_MBG7<-NULL
kappa_MBG7<-NULL

set.seed(123)


for (i in 1:5){
  df_folded <- fold(
    data = data,
    k = 5,
    id_col = "tiles",
    method = "n_dist",
  )
  
  data_f1<-df_folded[which(df_folded$.folds=="1"),]
  data_f2<-df_folded[which(df_folded$.folds=="2"),]
  data_f3<-df_folded[which(df_folded$.folds=="3"),]
  data_f4<-df_folded[which(df_folded$.folds=="4"),]
  data_f5<-df_folded[which(df_folded$.folds=="5"),]
  
  sample_f1<-stratified(data_f1, "G7",
                        nrow(data_f1[which(data_f1$G7=="1"),]))
  sample_f2<-stratified(data_f2, "G7",
                        nrow(data_f2[which(data_f2$G7=="1"),]))
  sample_f3<-stratified(data_f3, "G7",
                        nrow(data_f3[which(data_f3$G7=="1"),]))
  sample_f4<-stratified(data_f4, "G7",
                        nrow(data_f4[which(data_f4$G7=="1"),]))
  sample_f5<-stratified(data_f5, "G7",
                        nrow(data_f5[which(data_f5$G7=="1"),]))
  
  dataf[[i]]<- rbind(sample_f1,sample_f2,sample_f3,sample_f4,sample_f5)
  
  auc_MBG7[[i]]<-matrix(nrow=1, ncol=5)
  kappa_MBG7[[i]]<-matrix(nrow=1, ncol=5)
  
  for (fold in 1:5){
    training_set <- dataf[[i]][dataf[[i]]$.folds != fold,]
    testing_set <- dataf[[i]][dataf[[i]]$.folds == fold,]
    
    MBG7[[fold]]<- earth(formula_MBG7, data=training_set, degree=1, trace=1, glm=list(family=binomial))
    
    testing_set$score_MBG7<-predict(MBG7[[fold]], testing_set, type="response")
    auc_MBG7[[i]][,fold] <- as.vector(auc(testing_set$G7,testing_set$score_MBG7))
    testing_set$pred_MBG7<-factor( ifelse(testing_set$score_MBG7 < 0.5, "0", "1") )
    conf_MBG7 <- confusionMatrix(testing_set$pred_MBG7,testing_set$G7, positive="1")
    kappa_MBG7[[i]][,fold]<-conf_MBG7[["overall"]][["Kappa"]]
  }
}

auc_MBG7<-unlist(auc_MBG7)
kappa_MBG7<-unlist(kappa_MBG7)

##### Initialize Variables for the MBG8 Model #####

MBG8<-NULL

training_set<-NULL
testing_set<-NULL
dataf<-NULL
auc_MBG8<-NULL
kappa_MBG8<-NULL

set.seed(123)


for (i in 1:5){
  df_folded <- fold(
    data = data,
    k = 5,
    id_col = "tiles",
    method = "n_dist",
  )
  
  data_f1<-df_folded[which(df_folded$.folds=="1"),]
  data_f2<-df_folded[which(df_folded$.folds=="2"),]
  data_f3<-df_folded[which(df_folded$.folds=="3"),]
  data_f4<-df_folded[which(df_folded$.folds=="4"),]
  data_f5<-df_folded[which(df_folded$.folds=="5"),]
  
  sample_f1<-stratified(data_f1, "G8",
                        nrow(data_f1[which(data_f1$G8=="1"),]))
  sample_f2<-stratified(data_f2, "G8",
                        nrow(data_f2[which(data_f2$G8=="1"),]))
  sample_f3<-stratified(data_f3, "G8",
                        nrow(data_f3[which(data_f3$G8=="1"),]))
  sample_f4<-stratified(data_f4, "G8",
                        nrow(data_f4[which(data_f4$G8=="1"),]))
  sample_f5<-stratified(data_f5, "G8",
                        nrow(data_f5[which(data_f5$G8=="1"),]))
  
  dataf[[i]]<- rbind(sample_f1,sample_f2,sample_f3,sample_f4,sample_f5)
  
  auc_MBG8[[i]]<-matrix(nrow=1, ncol=5)
  kappa_MBG8[[i]]<-matrix(nrow=1, ncol=5)
  
  for (fold in 1:5){
    training_set <- dataf[[i]][dataf[[i]]$.folds != fold,]
    testing_set <- dataf[[i]][dataf[[i]]$.folds == fold,]
    
    MBG8[[fold]]<- earth(formula_MBG8, data=training_set, degree=1, trace=1, glm=list(family=binomial))
    
    testing_set$score_MBG8<-predict(MBG8[[fold]], testing_set, type="response")
    auc_MBG8[[i]][,fold] <- as.vector(auc(testing_set$G8,testing_set$score_MBG8))
    testing_set$pred_MBG8<-factor( ifelse(testing_set$score_MBG8 < 0.5, "0", "1") )
    conf_MBG8 <- confusionMatrix(testing_set$pred_MBG8,testing_set$G8, positive="1")
    kappa_MBG8[[i]][,fold]<-conf_MBG8[["overall"]][["Kappa"]]
  }
}

auc_MBG8<-unlist(auc_MBG8)
kappa_MBG8<-unlist(kappa_MBG8)

##### Initialize Variables for the MBG9 Model #####

MBG9<-NULL

training_set<-NULL
testing_set<-NULL
dataf<-NULL
auc_MBG9<-NULL
kappa_MBG9<-NULL

set.seed(123)


for (i in 1:5){
  df_folded <- fold(
    data = data,
    k = 5,
    id_col = "tiles",
    method = "n_dist",
  )
  
  data_f1<-df_folded[which(df_folded$.folds=="1"),]
  data_f2<-df_folded[which(df_folded$.folds=="2"),]
  data_f3<-df_folded[which(df_folded$.folds=="3"),]
  data_f4<-df_folded[which(df_folded$.folds=="4"),]
  data_f5<-df_folded[which(df_folded$.folds=="5"),]
  
  sample_f1<-stratified(data_f1, "G9",
                        nrow(data_f1[which(data_f1$G9=="1"),]))
  sample_f2<-stratified(data_f2, "G9",
                        nrow(data_f2[which(data_f2$G9=="1"),]))
  sample_f3<-stratified(data_f3, "G9",
                        nrow(data_f3[which(data_f3$G9=="1"),]))
  sample_f4<-stratified(data_f4, "G9",
                        nrow(data_f4[which(data_f4$G9=="1"),]))
  sample_f5<-stratified(data_f5, "G9",
                        nrow(data_f5[which(data_f5$G9=="1"),]))
  
  dataf[[i]]<- rbind(sample_f1,sample_f2,sample_f3,sample_f4,sample_f5)
  
  auc_MBG9[[i]]<-matrix(nrow=1, ncol=5)
  kappa_MBG9[[i]]<-matrix(nrow=1, ncol=5)
  
  for (fold in 1:5){
    training_set <- dataf[[i]][dataf[[i]]$.folds != fold,]
    testing_set <- dataf[[i]][dataf[[i]]$.folds == fold,]
    
    MBG9[[fold]]<- earth(formula_MBG9, data=training_set, degree=1, trace=1, glm=list(family=binomial))
    
    testing_set$score_MBG9<-predict(MBG9[[fold]], testing_set, type="response")
    auc_MBG9[[i]][,fold] <- as.vector(auc(testing_set$G9,testing_set$score_MBG9))
    testing_set$pred_MBG9<-factor( ifelse(testing_set$score_MBG9 < 0.5, "0", "1") )
    conf_MBG9 <- confusionMatrix(testing_set$pred_MBG9,testing_set$G9, positive="1")
    kappa_MBG9[[i]][,fold]<-conf_MBG9[["overall"]][["Kappa"]]
  }
}

auc_MBG9<-unlist(auc_MBG9)
kappa_MBG9<-unlist(kappa_MBG9)

#### Models performance statistics ####

kappa <-data.frame("MAG0"=as.vector(kappa_MAG0),
                   "MAG1"=as.vector(kappa_MAG1),
                   "MAG2"=as.vector(kappa_MAG2),
                   "MAG3"=as.vector(kappa_MAG3),
                   "MAG4"=as.vector(kappa_MAG4),
                   "MAG5"=as.vector(kappa_MAG5),
                   "MAG6"=as.vector(kappa_MAG6),
                   "MAG7"=as.vector(kappa_MAG7),
                   "MAG8"=as.vector(kappa_MAG8),
                   "MAG9"=as.vector(kappa_MAG9),
                   "MBG0"=as.vector(kappa_MBG0),
                   "MBG1"=as.vector(kappa_MBG1),
                   "MBG2"=as.vector(kappa_MBG2),
                   "MBG3"=as.vector(kappa_MBG3),
                   "MBG4"=as.vector(kappa_MBG4),
                   "MBG5"=as.vector(kappa_MBG5),
                   "MBG6"=as.vector(kappa_MBG6),
                   "MBG7"=as.vector(kappa_MBG7),
                   "MBG8"=as.vector(kappa_MBG8),
                   "MBG9"=as.vector(kappa_MBG9))

AUC <-data.frame("MAG0"=as.vector(auc_MAG0),
                 "MAG1"=as.vector(auc_MAG1),
                 "MAG2"=as.vector(auc_MAG2),
                 "MAG3"=as.vector(auc_MAG3),
                 "MAG4"=as.vector(auc_MAG4),
                 "MAG5"=as.vector(auc_MAG5),
                 "MAG6"=as.vector(auc_MAG6),
                 "MAG7"=as.vector(auc_MAG7),
                 "MAG8"=as.vector(auc_MAG8),
                 "MAG9"=as.vector(auc_MAG9),
                 "MBG0"=as.vector(auc_MBG0),
                 "MBG1"=as.vector(auc_MBG1),
                 "MBG2"=as.vector(auc_MBG2),
                 "MBG3"=as.vector(auc_MBG3),
                 "MBG4"=as.vector(auc_MBG4),
                 "MBG5"=as.vector(auc_MBG5),
                 "MBG6"=as.vector(auc_MBG6),
                 "MBG7"=as.vector(auc_MBG7),
                 "MBG8"=as.vector(auc_MBG8),
                 "MBG9"=as.vector(auc_MBG9))


#### Susceptibility map ####

load("~/Library/CloudStorage/OneDrive-UNIPA/Papers/turkey2023/data/data.RData")

sample_G0<-NULL
for (i in 1:10){
  sample_G0[[i]]<-stratified(data, "G0",.5*nrow(data[which(data$G0==1),]))
}
sample_G0<-do.call("rbind", sample_G0)

sample_G3<-NULL
for (i in 1:10){
  sample_G3[[i]]<-stratified(data, "G3",.5*nrow(data[which(data$G3==1),]))
}
sample_G3<-do.call("rbind", sample_G3)

sample_G6<-NULL
for (i in 1:10){
  sample_G6[[i]]<-stratified(data, "G6",.5*nrow(data[which(data$G6==1),]))
}
sample_G6<-do.call("rbind", sample_G6)

sample_G9<-NULL
for (i in 1:10){
  sample_G9[[i]]<-stratified(data, "G9",.5*nrow(data[which(data$G9==1),]))
}
sample_G9<-do.call("rbind", sample_G9)


MBG6<- earth(G6 ~ ELE+SLO+GEO+PLC+DEV+EP+SPI+TWI+LSF+GORD, data=sample_G6, degree=1, trace=1, glm=list(family=binomial))

#data$score<-predict(MBG6, data, type="response")
#auc_data <- as.vector(auc(data$G6,data$score))
#data$pred_MBG6<-factor( ifelse(data$score < 0.5, "0", "1") )
#conf_data <- confusionMatrix(data$pred_MBG6,data$G6, positive="1")
#kappa_data<-conf_data[["overall"]][["Kappa"]]

library(raster)
score <- rasterFromXYZ(data[,c("x","y","score")])  
crs(score)<-"EPSG:26914"

writeRaster(score,
            file="/Users/christian/Library/CloudStorage/OneDrive-UNIPA/Papers/turkey2023/data/score",filetype="GTiff", overwrite=TRUE, NAflag=-99999)




