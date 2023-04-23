"""
Adjust the lat and lon of Gao&O'neil data from 2010 to 2100;
Gao & Oneil dataset gives the fraction of urban with respect to land area.
Calculate the fraction of urban with respect to the grid cell area.
"""
import netCDF4 as nc
import numpy as np
from numpy import dtype
import time
import gc
t0 = time.time()

# load land area file
Area=nc.Dataset('/glade/scratch/keerzhang/Temp_RawData/AncillaryLayers/land_area_km2.nc')
areaunmask=Area.variables["land_area_km2"][:,:] # this is land area, not grid area

# Treat no data value (-3.402823E38) as zero land area
area=np.ma.where(areaunmask<1e-6,0,areaunmask)

# the latitude of land area of Gao Oneil dataset is from 85~-58
lat0=Area.variables["lat"][:]

# the 0.00833*0.00833 deg grid on equator is 0.8589 km2
# the area of a grid at other latitudes is scaled by cos(lat)
weight=np.cos((lat0/180.0)*np.pi)*0.8589 
Gridarea=np.tile(weight,(43200,1)).T # grid area
print(Gridarea.shape)

# read the defaul MD urban raw data based on the Jackson 2010 dataset
DefaultMD = nc.Dataset('/glade/work/keerzhang/surfacedata/THESIS/med_dens_1km.nc')
Dlat=DefaultMD.variables["lat"][:] # 84~-90
Dlon=DefaultMD.variables["lon"][:] # 0~360
Dmd=DefaultMD.variables["md"][:]

nlat = 20880
nlon = 43200
# specify SSP scenario 
Scenario='1'

for year in np.arange(2010,2110,10):
    # create the new 1km Gao&Oneil data as NETCDF4 files
    NewMD = nc.Dataset('/glade/scratch/keerzhang/Temp_RawData/RawSSP'+Scenario+'/ssp'+Scenario+'_'+str(year)+'_ad_norm_c20220416.nc', 'w', format='NETCDF4')
    GaoMD = nc.Dataset('/glade/scratch/keerzhang/Temp_RawData/RawSSP'+Scenario+'/ssp'+Scenario+'_'+str(year)+'.nc')
    Nlat=GaoMD.variables["lat"][:]  # 85~-58
    Nlon=GaoMD.variables["lon"][:] # -180~180
    
    # create lat lon dimensions
    lat = NewMD.createDimension('lat', nlat)
    lon = NewMD.createDimension('lon', nlon)    
    
    # create variable array
    VAR = NewMD.createVariable('md','f4', ('lat', 'lon'),fill_value=1e+36) #,zlib=True,complevel=3
    VAR.long_name = "md"
    VAR.crs = "+init=epsg:4326"
    VAR.coordinates = "lon lat"
    
    # create latitude
    lat = NewMD.createVariable('lat', dtype('double').char, ('lat'))
    lat.long_name = "latitude coordinate"
    lat.standard_name = "latitude"
    lat.units = "degrees_north"

    # create longitude
    lon = NewMD.createVariable('lon', dtype('double').char, ('lon'))
    lon.long_name = "longitude coordinate"
    lon.standard_name = "longitude";
    lon.units = "degrees_east";

    Nmd0=GaoMD.variables["ssp"+Scenario+'_'+str(year)][:]
    # Treat no data value (-3.402823E38) as zero urban fraction   
    Nmd0=np.ma.where(Nmd0<1e-6,0,Nmd0) # negative values means missing value, set them as zeros
    
    Urbanarea=Nmd0*area # the true urban area is the land area * urban fraction
    Nmd=Urbanarea/Gridarea # calculate the urban fraction with repsect to the grid area (this is important for coastal grids)
    Nmd=np.ma.where(Nmd>1,1,Nmd) # urban fraction cannot exceed 1
    del Nmd0; gc.collect()

    # adjust longitude
    # new Gao & Oneil data uses -180 to 180
    # default md data uses 0 to 360
    Nlon=np.ma.where(Nlon<0,360+Nlon,Nlon) # 180~360 + 0~180    
    NmdLon=np.ma.concatenate((Nmd[:,21600:43200],Nmd[:,0:21600]),axis=1)
    Nlonlon=np.ma.concatenate((Nlon[21600:43200],Nlon[0:21600]),axis=0)
    
    # adjust latitude
    # new Gao & Oneil data uses 84.9958 to -57.9958
    # default md data uses 83.9958 to -89.9958
    AllZeros=np.zeros(Dmd[17040:20880,:].shape)
    NmdLonlat0=NmdLon[120:17160,:]
    NmdLonlat=np.ma.concatenate((NmdLonlat0,AllZeros),axis=0)
    
    Nlatlat0=Nlat[120:17160]
    Nlatlat=np.ma.concatenate((Nlatlat0,Dlat[17040:20880]),axis=0)
    
    # the ocean grid has 0 values; apply the landmask in the default data
    NmdLonlat=np.ma.array(NmdLonlat,mask=Dmd.mask)   
    VAR[:] = NmdLonlat
    lat[:] = Nlatlat
    lon[:] = Nlonlon
    NewMD.close()
                       
    t1 = time.time()
    total = t1-t0
    print('---------Finished year '+str(year)+'------------------')  
