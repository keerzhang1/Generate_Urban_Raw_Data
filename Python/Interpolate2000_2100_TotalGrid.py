"""
Interpolate the decadal urban raw data to get annual urban data
"""
import netCDF4 as nc 

# Interpolate the decadal urban raw data under different SSP scenarios
for k in range(1,6):
    Scenario=str(k)
    for i in range(2,11,1):
        NewRawThisDecade = nc.Dataset('/glade/scratch/keerzhang/THESIS_pcturb_grid/urban_properties_180622_release/Output/SSP'+Scenario+'/urban_properties_GaoOneil_05deg_ThreeClass_ssp'+Scenario+'_2'+str(i).zfill(2)+'0_c20220910.nc')#   
        NewRawLastDecade = nc.Dataset('/glade/scratch/keerzhang/THESIS_pcturb_grid/urban_properties_180622_release/Output/SSP'+Scenario+'/urban_properties_GaoOneil_05deg_ThreeClass_ssp'+Scenario+'_2'+str(i-1).zfill(2)+'0_c20220910.nc')

        PCT_URBAN_ThisDecade=NewRawThisDecade.variables["PCT_URBAN"][:,:,:] # this decade
        PCT_URBAN_LastDecade=NewRawLastDecade.variables["PCT_URBAN"][:,:,:] # last decade
        if i==2: #calculate urban raw data from 2015 to 2020
            for j in range(5,10):
                PCT_URBAN_Inter=(PCT_URBAN_LastDecade*(10-j)+PCT_URBAN_ThisDecade*j)/10.0
                NewRaw20XX_Inter = nc.Dataset('/glade/scratch/keerzhang/THESIS_pcturb_grid/urban_properties_180622_release/Output/SSP'+Scenario+'/urban_properties_GaoOneil_05deg_ThreeClass_ssp'+Scenario+'_20'+str(10*i-10+j)+'_c20220910.nc', 'r+', format='NETCDF4_CLASSIC')#   
                NewRaw20XX_Inter.variables["PCT_URBAN"][:]=PCT_URBAN_Inter
                NewRaw20XX_Inter.close()
        else:
            for j in range(1,10):
                PCT_URBAN_Inter=(PCT_URBAN_LastDecade*(10-j)+PCT_URBAN_ThisDecade*j)/10.0
                NewRaw20XX_Inter = nc.Dataset('/glade/scratch/keerzhang/THESIS_pcturb_grid/urban_properties_180622_release/Output/SSP'+Scenario+'/urban_properties_GaoOneil_05deg_ThreeClass_ssp'+Scenario+'_20'+str(10*i-10+j)+'_c20220910.nc', 'r+', format='NETCDF4_CLASSIC')# 
                NewRaw20XX_Inter.variables["PCT_URBAN"][:]=PCT_URBAN_Inter
                NewRaw20XX_Inter.close()
        print('Finished interpolating ssp'+Scenario+' files from year 20'+str(10*i-10)+' to year 20'+str(10*i))

# Interpolate historical files
NewRaw2000=nc.Dataset('/glade/scratch/3keerzhang/THESIS_pcturb_grid/urban_properties_180622_release/Output/Historical/urban_properties_GaoOneil_05deg_ThreeClass_2000_c20220910.nc')
PCT_URBAN2000=NewRaw2000.variables["PCT_URBAN"][:,:,:]

# calculate the average urban land cover across five SSP scenarios in years 2010 and 2020
for k in range(1,6):
    Scenario=str(k)
    NewRaw2010=nc.Dataset('/glade/scratch/keerzhang/THESIS_pcturb_grid/urban_properties_180622_release/Output/SSP'+Scenario+'/urban_properties_GaoOneil_05deg_ThreeClass_ssp'+Scenario+'_2010_c20220910.nc')
    NewRaw2020=nc.Dataset('/glade/scratch/keerzhang/THESIS_pcturb_grid/urban_properties_180622_release/Output/SSP'+Scenario+'/urban_properties_GaoOneil_05deg_ThreeClass_ssp'+Scenario+'_2020_c20220910.nc')
    PCT_URBAN2010=NewRaw2010.variables["PCT_URBAN"][:,:,:] 
    PCT_URBAN2020=NewRaw2020.variables["PCT_URBAN"][:,:,:]   
    if k==1:
        PCT_URBAN2010Mean=PCT_URBAN2010/5.0
        PCT_URBAN2020Mean=PCT_URBAN2020/5.0
    if k>1:
        PCT_URBAN2010Mean=PCT_URBAN2010Mean+PCT_URBAN2010/5.0
        PCT_URBAN2020Mean=PCT_URBAN2020Mean+PCT_URBAN2020/5.0
        
#generate historical urban raw data from 2001 to 2010     
for j in range(1,11): 
    PCT_URBAN_Inter=(PCT_URBAN2000*(10-j)+PCT_URBAN2010Mean*j)/10.0
    NewRaw20XX_Inter = nc.Dataset('/glade/scratch/keerzhang/THESIS_pcturb_grid/urban_properties_180622_release/Output/Historical/urban_properties_GaoOneil_05deg_ThreeClass_20'+str(j).zfill(2)+'_c20220910.nc', 'r+', format='NETCDF4_CLASSIC')# 
    NewRaw20XX_Inter.variables["PCT_URBAN"][:]=PCT_URBAN_Inter
    NewRaw20XX_Inter.close()   
    
#generate historical urban raw data from 2011 to 2014     
for j in range(1,5):
    PCT_URBAN_Inter=(PCT_URBAN2010Mean*(10-j)+PCT_URBAN2020Mean*j)/10.0
    NewRaw20XX_Inter = nc.Dataset('/glade/scratch/keerzhang/THESIS_pcturb_grid/urban_properties_180622_release/Output/Historical/urban_properties_GaoOneil_05deg_ThreeClass_20'+str(10+j)+'_c20220910.nc', 'r+', format='NETCDF4_CLASSIC')#   
    NewRaw20XX_Inter.variables["PCT_URBAN"][:]=PCT_URBAN_Inter
    NewRaw20XX_Inter.close()
print('Finished interpolating historical urban raw data from year 2000 to 2014')
