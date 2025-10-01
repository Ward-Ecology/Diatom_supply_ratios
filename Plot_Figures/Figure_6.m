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

f6 = figure(6);
f6.Position = [74 748 1179 522];
clf

% colormap % https://coolors.co/palettes/popular/5%20colors
clrs = [ 61  64  91 ;...  % prov i   ('#3d405b')
        129 178 154 ;...  % prov ii  ('#81b29a')
        242 204 143 ;...  % prov iii ('#f2cc8f')
        224 122  95];     % prov iv   ('#e07a5f')
clrs=clrs./255;

% https://coolors.co/palette/f4f1de-e07a5f-3d405b-81b29a-f2cc8f

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

axh = axesm('MapProjection','mollweid',...
                'plinelocation',90,...
                'mlinelocation',Inf,...
                'MapLatLimit',[-80 80],...
                'MapLonLimit',[1 359]+130,...
                'Frame', 'on',...
                'Grid', 'off');
contourfm(Lat,Lon,Provinces',1.5:3.5,'LineW',2)
hold on
gamma = EcosystemData.gamma.Annual(:,:,1);
gamma(gamma==0)=NaN;
contourm(Lat,Lon,gamma',[0.75 0.75],'w--','LineW',2)

colormap(clrs)
caxis([0.5 4.5])
ch = colorbar;
ch.Ticks = 1:4;
ch.TickLabels = {'Province i \newline{• \phi_{Si:N}>\theta_{Si:N,D}} \newline{• \phi_{N:Fe}<\theta_{N:Fe,D}}',...
                 'Province ii \newline{• \phi_{Si:N}<\theta_{Si:N,D}} \newline{• \phi_{N:Fe}<\theta_{N:Fe,D}}',...
                 'Province iii \newline{• \phi_{Si:Fe}<\theta_{Si:Fe,D}} \newline{• \phi_{N:Fe}>\theta_{N:Fe,D}}',...
                 'Province iv \newline{• \phi_{Si:Fe}>\theta_{Si:Fe,D}} \newline{• \phi_{N:Fe}>\theta_{N:Fe,D}}'};
set(ch, 'YDir', 'reverse' );
box off
landareas = shaperead('landareas.shp','UseGeoCoords',true);
geoshow(landareas,...
        'FaceColor',[0.75 0.75 0.75],...
        'EdgeColor','k',...
        'LineWidth',0.5);
title('Provinces')
set(gcf,'Color','w')
axis off

%%

set(findall(gcf,'-property','FontSize'),'FontSize',14)
exportgraphics(gcf,['../Figures/Figure_6.png'],'Resolution',600)



