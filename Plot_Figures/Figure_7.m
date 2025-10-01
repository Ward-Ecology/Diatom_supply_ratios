clear 
clc

% simulation_name = 'RUN33_6_NEWEDES_HAM';
sim_name = 'RUN33_6_NEWEDES';

% data path 
data_path = '../Data';

%% load data

load([data_path '/' sim_name '_Ecosystem.mat'])
load([data_path '/' sim_name '_Physics.mat'])

clear data_path sim_name

Lat = EcosystemData.Y;
Lon = linspace(0,360,numel(EcosystemData.X));

%% calculate provinces in global model
Provinces = zeros(size(EcosystemData.phiNFe.Annual)).*NaN;
Provinces(EcosystemData.phiNFe.Annual<1 & EcosystemData.phiSiN.Annual >1) = 1;
Provinces(EcosystemData.phiNFe.Annual<1 & EcosystemData.phiSiN.Annual <1) = 2;
Provinces(EcosystemData.phiNFe.Annual>1 & EcosystemData.phiSiFe.Annual<1) = 3;
Provinces(EcosystemData.phiNFe.Annual>1 & EcosystemData.phiSiFe.Annual>1) = 4;
%%
r = EcosystemData.gammaFe.Annual(:,:,1);
g = EcosystemData.gammaSi.Annual(:,:,1);
b = EcosystemData.gammaN.Annual(:,:,1);

rgb = (1-cat(3,r',g',b')).^0.5;

rgb = 0.8.*rgb./max(rgb,[],3);

clf
image(rgb)
axis xy


%%

f6 = figure(6);
f6.Position = [74 801 1057 469];
clf

t = tiledlayout(1,5, 'TileSpacing','compact', ...
                     'Padding','compact');
nexttile([1 4])

axh = axesm('MapProjection','mollweid',...
                'plinelocation',90,...
                'mlinelocation',Inf,...
                'MapLatLimit',[-80 80],...
                'MapLonLimit',[1 359]+130,...
                'Frame', 'on',...
                'Grid', 'off');

sh=surfm(Lat,Lon,zeros(size(rgb,[1 2])))';
sh.CData = rgb;



contourm(Lat,Lon,Provinces',1.5:3.5,'k','LineW',2)
hold on

Browning_Incubations

gamma = EcosystemData.gamma.Annual(:,:,1);
gamma(gamma==0)=NaN;
contourm(Lat,Lon,gamma',[0.75 0.75],'w--','LineW',2)



box off
landareas = shaperead('landareas.shp','UseGeoCoords',true);
geoshow(landareas,...
        'FaceColor',[0.75 0.75 0.75],...
        'EdgeColor','k',...
        'LineWidth',0.5);
title('Limiting factors')
set(gcf,'Color','w')
axis off

nexttile
colortriangle(101)
axis off

%%

set(findall(gcf,'-property','FontSize'),'FontSize',14)
exportgraphics(gcf,['../Figures/Figure_7.png'],'Resolution',600)



