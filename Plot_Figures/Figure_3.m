clear 
clc

% simulation_name = 'RUN33_6_NEWEDES_HAM';
sim_name = 'RUN33_6_NEWEDES';

% data path 
data_path = '../Data';

%% load data

load([data_path '/' sim_name '_Ecosystem.mat'])
load([data_path '/' sim_name '_Physics.mat'])
load([data_path '/WOA.mat'])
load([data_path '/Horstmann.mat'])

clear data_path sim_name

%% Get WOA Si*
Sistar.WOA.surface = WOA.silicate.surface - WOA.nitrate.surface;
Sistar.WOA.z100    = WOA.silicate.z100    - WOA.nitrate.z100   ;
Sistar.WOA.lat     = WOA.nitrate.lat;
Sistar.WOA.lon     = linspace(-180,180,numel(WOA.nitrate.lon));

%% Get MITgcm Si*


i_above = find(EcosystemData.Z<=100,1,'last');
z_above = EcosystemData.Z(i_above);
i_below = find(EcosystemData.Z>=100,1,'first');
z_below = EcosystemData.Z(i_below);

d1 = 100 - z_above;    % distance from above to target
d2 = z_below - 100;    % distance from target to below
d3 = d1 + d2;          % total depth difference

% interpolate Si
Si_above = EcosystemData.SiO2.Annual(:,:,i_above);
Si_below = EcosystemData.SiO2.Annual(:,:,i_below);

Si100 = (Si_above * d2 + Si_below * d1) ./ d3;

% interpolate N
N_above = EcosystemData.NO3.Annual(:,:,i_above);
N_below = EcosystemData.NO3.Annual(:,:,i_below);

N100 = (N_above * d2 + N_below * d1) ./ d3;


Sistar.MITgcm.surface = EcosystemData.SiO2.Annual(:,:,1) ...
                      - EcosystemData.NO3.Annual(:,:,1);
Sistar.MITgcm.z100    = Si100 - N100;
Sistar.MITgcm.lat     = EcosystemData.Y;
Sistar.MITgcm.lon     = linspace(0,360,numel(EcosystemData.X));

%% calculate provinces in global model
Provinces = zeros(size(Sistar.MITgcm.surface)).*NaN;

EcosystemData.phiNFe.Annual(EcosystemData.phiNFe.Annual==0)=NaN;
EcosystemData.phiSiFe.Annual(EcosystemData.phiSiFe.Annual==0)=NaN;
EcosystemData.phiSiN.Annual(EcosystemData.phiSiN.Annual==0)=NaN;

Provinces(EcosystemData.phiNFe.Annual<1 & EcosystemData.phiSiN.Annual >1) = 1;
Provinces(EcosystemData.phiNFe.Annual<1 & EcosystemData.phiSiN.Annual <1) = 2;
Provinces(EcosystemData.phiNFe.Annual>1 & EcosystemData.phiSiFe.Annual<1) = 3;
Provinces(EcosystemData.phiNFe.Annual>1 & EcosystemData.phiSiFe.Annual>1) = 4;



% figure(1)
% imagesc(Provinces')
% hold on
% contour(Provinces',[1:4],'k','LineW',1)
% axis xy,colorbar
% colormap(parula(5)),caxis([0 5])
%% Plot Figure 3a

f3 = figure(3);
clf
f3.Position = [74 394 1662 876];


t = tiledlayout(2,2, 'TileSpacing','compact',...
                     'Padding','compact', ...
                     'TileIndexing', 'columnmajor');
nexttile
axh = axesm('MapProjection','mollweid',...
                'plinelocation',90,...
                'mlinelocation',Inf,...
                'MapLatLimit',[-90 90],...
                'MapLonLimit',[1 359]+130,...
                'Frame', 'on',...
                'Grid', 'off');
pcolorm(Sistar.WOA.lat,Sistar.WOA.lon,Sistar.WOA.z100')
colormap(redblue)
climits = [-30 30];
caxis(climits)
box off
landareas = shaperead('landareas.shp','UseGeoCoords',true);
geoshow(landareas,...
        'FaceColor',[0.75 0.75 0.75],...
        'EdgeColor','k',...
        'LineWidth',0.5);
title('(a) WOA 100m Si^*')
set(gcf,'Color','w')
axis off

%% Plot Figure 3b
nexttile
axh = axesm('MapProjection','mollweid',...
                'plinelocation',90,...
                'mlinelocation',Inf,...
                'MapLatLimit',[-90 90],...
                'MapLonLimit',[1 359]+130,...
                'Frame', 'on',...
                'Grid', 'off');
pcolorm(Sistar.WOA.lat,Sistar.WOA.lon,Sistar.WOA.surface')
hold on
plotm(HPLC.pigments.Latitude,HPLC.pigments.Longitude,'k','LineW',2)
colormap(redblue)
climits = [-30 30];
caxis(climits)
box off
landareas = shaperead('landareas.shp','UseGeoCoords',true);
geoshow(landareas,...
        'FaceColor',[0.75 0.75 0.75],...
        'EdgeColor','k',...
        'LineWidth',0.5);
title('(b) WOA surface Si^*')
set(gcf,'Color','w')
axis off

%% Plot Figure 3c
nexttile
axh = axesm('MapProjection','mollweid',...
                'plinelocation',90,...
                'mlinelocation',Inf,...
                'MapLatLimit',[-90 90],...
                'MapLonLimit',[1 359]+130,...
                'Frame', 'on',...
                'Grid', 'off');
pcolorm(Sistar.MITgcm.lat,Sistar.MITgcm.lon,Sistar.MITgcm.z100')
hold on
contourm(Sistar.MITgcm.lat,Sistar.MITgcm.lon,Provinces',1.5:3.5,'k','LineW',1)

gamma = EcosystemData.gamma.Annual(:,:,1);
gamma(gamma==0)=NaN;
contourm(Sistar.MITgcm.lat,Sistar.MITgcm.lon,gamma',0.75,'m--','LineW',1)

colormap(redblue)
climits = [-30 30];
caxis(climits)
box off
landareas = shaperead('landareas.shp','UseGeoCoords',true);
geoshow(landareas,...
        'FaceColor',[0.75 0.75 0.75],...
        'EdgeColor','k',...
        'LineWidth',0.5);
title('(c) Model 100m Si^*')
axis off

% textm(lat,lon,'txt')
textm(20,-64,'i')
textm(-22,-25,'i')
textm(-26,72,'i')
textm(5,-30,'ii')
textm(30,-180,'ii')
textm(-35,-170,'ii')
textm(0,60,'ii')
textm(0,-110,'iii')
textm(-45,-105,'iii')
textm(-9,-150,'iv')

%% Plot Figure 3d
nexttile
axh = axesm('MapProjection','mollweid',...
                'plinelocation',90,...
                'mlinelocation',Inf,...
                'MapLatLimit',[-90 90],...
                'MapLonLimit',[1 359]+130,...
                'Frame', 'on',...
                'Grid', 'off');
pcolorm(Sistar.MITgcm.lat,Sistar.MITgcm.lon,Sistar.MITgcm.surface')
hold on
contourm(Sistar.MITgcm.lat,Sistar.MITgcm.lon,Provinces',1.5:3.5,'k','LineW',1)
contourm(Sistar.MITgcm.lat,Sistar.MITgcm.lon,gamma',0.75,'m--','LineW',1)
lat = linspace(min(HPLC.pigments.Latitude'),max(HPLC.pigments.Latitude'),100);
lon = ones(1,100).*mean(HPLC.pigments.Longitude);
plotm(lat,lon,'k','LineW',2)

colormap(redblue)
climits = [-30 30];
caxis(climits)
box off
landareas = shaperead('landareas.shp','UseGeoCoords',true);
geoshow(landareas,...
        'FaceColor',[0.75 0.75 0.75],...
        'EdgeColor','k',...
        'LineWidth',0.5);
title('(d) Model surface Si^*')
axis off

% textm(lat,lon,'txt')
textm(20,-64,'i')
textm(-22,-25,'i')
textm(-26,72,'i')
textm(5,-30,'ii')
textm(30,-180,'ii')
textm(-35,-170,'ii')
textm(0,60,'ii')
textm(0,-110,'iii')
textm(-45,-105,'iii')
textm(-9,-150,'iv')

ch = colorbar;
ch.Position = [0.5 0.3 0.0172 0.4229];
ch.TickLabels{end} = '30 µM';
% Export figure to file
%%

set(findall(gcf,'-property','FontSize'),'FontSize',14)


exportgraphics(gcf,['../Figures/Figure_3.png'],'Resolution',600)