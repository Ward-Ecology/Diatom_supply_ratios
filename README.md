Global data grids are not included in this repository

  **World Ocean Atlas 2023**
  
    Need to get these files...
    
    woa23_nitrate_annual.nc
    woa23_silicate_annual.nc
    
    and place in Original_data/WOA2023/.
  
  **MITgcm 'Darwin' ecosystem model**
  
    Need to get these files...
    
    3d.0000026160.nc
    3d.0000026400.nc
    3d.0000026640.nc
    3d.0000026880.nc
    3d.0000027120.nc
    3d.0000027360.nc
    3d.0000027600.nc
    3d.0000027840.nc
    3d.0000028080.nc
    3d.0000028320.nc
    3d.0000028560.nc
    3d.0000028800.nc
    songmiao_mineral_mth-2d.bin
    
    and place in Original_data/MITgcm/RUN33_6_NEWEDES/.

**Then run the following scripts in Collate_data**

   Collate_MITgcm.m
   Calculate_fluxes_explicit.m
   Collate_Horstmann.m
   Collate_WOA.m

Then run the following in Plot_Figures

   Figure_1.m
   Figure_2.m
   Figure_3.m
   Figure_4.m
   Figure_5.m
   Figure_6.m
   Figure_7.m
   Figure_8.m
   Figure_B1.m
