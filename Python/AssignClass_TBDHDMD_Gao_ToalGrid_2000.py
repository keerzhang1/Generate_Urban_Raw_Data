"""
Divide the total urban area into three density classes (TBD, MD, and HD)
"""
import numpy as np
import netCDF4 as nc 
import time
t0 = time.time()

# load the default CLM surface data (Jackson 2010)
CLM2000=nc.Dataset('/glade/p/cesm/cseg/inputdata/lnd/clm2/rawdata/mksrf_urban_0.05x0.05_simyr2000.c170724.nc')
REGION_IDNew=CLM2000.variables["REGION_ID"][:]
REGION_IDNew3=np.tile(REGION_IDNew,(3,1,1))
UrbanDensityCLM=CLM2000.variables["PCT_URBAN"][:,:,:]

# calculate the mean HD/MD Ratio for each urban region
HMRatio=np.zeros(33)
for j in range(1,34):
    UrbanDensitySub=np.ma.masked_where(REGION_IDNew3!=j,UrbanDensityCLM)
    HD=np.nansum(UrbanDensitySub[1])
    MD=np.nansum(UrbanDensitySub[2])
    HMRatio[j-1]=HD/MD
    
NewRaw2020OriThisDecade = nc.Dataset('/glade/scratch/keerzhang/THESIS_pcturb_grid/urban_properties_180622_release/Output/SSP1/urban_properties_data.05deg_2000_c20220910.nc')#   
PCT_URBAN_NewThisDecade=np.nansum(NewRaw2020OriThisDecade.variables["PCT_URBAN"][:,:,:],axis=0)
PCT_URBAN_NewThisDecadeNoDown=PCT_URBAN_NewThisDecade

# Divide the total Gao & Oneil urban areas into TBD, HD, MD
# Assumption: if there is urban area in both Jackson data and Gao & Oneil data, the
# Gao & Oneil urban area is partitioned using the same TBD/HD/MD ratio of that grid in the J2010 dataset
# if there is urban area in Gao & Oneil data but no urban area in Jackson data, use the mean HD/MD ratio of that urban region
# if there is no urban in Gao & Oneil data, PCT_URBAN = 0   
UrbanDensityNew=np.zeros(UrbanDensityCLM.shape)
ANum=np.zeros(4)
for j in range(0,UrbanDensityCLM.shape[1]):
    for k in range(0,UrbanDensityCLM.shape[2]):
        GridCLM=UrbanDensityCLM[:,j,k]
        TotalNew=PCT_URBAN_NewThisDecadeNoDown[j,k]
        if np.sum(GridCLM)>0 and TotalNew>0:
            UrbanDensityNew[:,j,k]=TotalNew*GridCLM/np.sum(GridCLM)
            ANum[0]=ANum[0]+1
        elif np.sum(GridCLM)==0 and TotalNew>0:
            ID=REGION_IDNew[j,k]
            UrbanDensityNew[1,j,k]=TotalNew*HMRatio[ID-1]/(HMRatio[ID-1]+1)
            UrbanDensityNew[2,j,k]=TotalNew/(HMRatio[ID-1]+1)
            ANum[1]=ANum[1]+1
        elif np.sum(GridCLM)>0 and TotalNew==0:
            ANum[2]=ANum[2]+1
            UrbanDensityNew[:,j,k]=0.0
        elif np.sum(GridCLM)==0 and TotalNew==0:
            ANum[3]=ANum[3]+1
            UrbanDensityNew[:,j,k]=0.0
            
# create this new file by copying the 0.05 deg urban raw data where all urban areas are considered MD
NewRaw2020Adjusted = nc.Dataset('/glade/scratch/keerzhang/THESIS_pcturb_grid/urban_properties_180622_release/Output/SSP1/urban_properties_GaoOneil_05deg_ThreeClass_2000_c20220910.nc', 'r+', format='NETCDF4_CLASSIC')# 
NewRaw2020Adjusted.variables["PCT_URBAN"][:,:,:]=UrbanDensityNew[:,:,:]
NewRaw2020Adjusted.title="0.05 deg urban raw data (Gao & O'neil)"
NewRaw2020Adjusted.comment="PCT_URBAN was calculated with repect to the total grid cell."
NewRaw2020Adjusted.close()

t1 = time.time()
total = t1-t0             
print('---------------At Year 2000----------------')
print(ANum)
print(total)
