clear
clc

WOA_path = '../Original_data/WOA2023/';

WOA.nitrate.data  =         ncread([WOA_path 'woa23_nitrate_annual.nc'],'n_an' ) ;
WOA.nitrate.lat    = double(ncread([WOA_path 'woa23_nitrate_annual.nc'],'lat'  ));
WOA.nitrate.lon    = double(ncread([WOA_path 'woa23_nitrate_annual.nc'],'lon'  ));
WOA.nitrate.depth  = double(ncread([WOA_path 'woa23_nitrate_annual.nc'],'depth'));

WOA.nitrate.surface = WOA.nitrate.data(:,:,1);
i100 = find(WOA.nitrate.depth==100);
WOA.nitrate.z100    = WOA.nitrate.data(:,:,i100);

WOA.silicate.data  =        ncread([WOA_path 'woa23_silicate_annual.nc'],'i_an' ) ;
WOA.silicate.lat   = double(ncread([WOA_path 'woa23_silicate_annual.nc'],'lat'  ));
WOA.silicate.lon   = double(ncread([WOA_path 'woa23_silicate_annual.nc'],'lon'  ));
WOA.silicate.depth = double(ncread([WOA_path 'woa23_silicate_annual.nc'],'depth'));

WOA.silicate.surface = WOA.silicate.data(:,:,1);
i100 = find(WOA.silicate.depth==100);
WOA.silicate.z100    = WOA.silicate.data(:,:,i100);



save(['../Data/WOA.mat'],'WOA','-v7.3')