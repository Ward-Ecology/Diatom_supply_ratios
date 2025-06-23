clear 
clc

% simulation_name = 'RUN33_6_NEWEDES';
sim_name = 'RUN33_6_NEWEDES';

% WOA path 
data_path = '../Data';

%% load data

load([data_path '/' sim_name '_Ecosystem.mat'])
load([data_path '/' sim_name '_Physics.mat'])
load([data_path '/WOA.mat'])

clear data_path sim_name

% get WOA data
WOA_lat   = WOA.nitrate.lat;
WOA_lon   = WOA.nitrate.lon;
WOA_lon   = linspace(-180,180,numel(WOA_lon));
WOA_depth = WOA.nitrate.depth;
Si_star   = WOA.silicate.data - WOA.nitrate.data;

% get model data
DiatomC = EcosystemData.Diatom.Annual;
PtotalC = EcosystemData.Prokaryote.Annual ...
        + EcosystemData.Picoeukaryote.Annual ...
        + EcosystemData.Coccolithophore.Annual ...
        + EcosystemData.Diazotroph.Annual ...
        + EcosystemData.Diatom.Annual ...
        + EcosystemData.Mixotroph.Annual;
drF     = EcosystemData.drF;
mod_lon = EcosystemData.X;
mod_lat = EcosystemData.Y;

DiatomC_integrated = tensorprod(DiatomC,drF,3,1);
PtotalC_integrated = tensorprod(PtotalC,drF,3,1);

DiatFrac_integrated = DiatomC_integrated./PtotalC_integrated;
%% Plot Figure 1

f1 = figure(1);
clf
f1.Position = [74 576 1097 694];

% background color for areas with no data
rectangle('Position',[-2 -2 4 4].*0.98,...
          'Curvature',1,...
          'FaceColor',[0.25 0.25 0.25])
hold on
% set up axes
axh = axesm('MapProjection','eqaazim',...
            'Origin',[-46 -86],... % Antipode of continental pole of inaccessibility
            'plinelocation',90,...
            'mlinelocation',0,...
            'Frame', 'on', 'Grid', 'off');
box off

[~,iz] = min(abs(WOA_depth-100));
pcolorm(WOA_lat-0.5,WOA_lon-1,Si_star(:,:,iz)')

[~,ch] = contourm(mod_lat,mod_lon,DiatomC_integrated',20:20:100,'LineW',1);

landareas = shaperead('landareas.shp','UseGeoCoords',true);
geoshow(landareas,...
        'FaceColor',[0.75 0.75 0.75],...
        'EdgeColor','k',...
        'LineWidth',0.5);

colormap(redblue)

climits = [-50 50];
caxis(climits)
h =colorbar;
h.Position = [[0.8 0.3 0.03 0.4]];
clim(climits)
h.Limits=climits;

axis off




title({'WOA Si^* (mmol m^{-3}) at 100 m and',...
    'Modelled depth-integrated diatom biomass (mmol C m^{-2})'})
set(gcf,'Color','w')

drawnow
for i=1:5
    ch.Children(i).Color = [1 1 1].*(i-1)/5 ;
end

lh=legend(ch.Children,'100 mm C m^{-2}','80','60','40','20');
legend('boxoff')



set(findall(gcf,'-property','FontSize'),'FontSize',14)
lh.Position = [0.771 0.725 0.1889 0.2026];

exportgraphics(gcf,['../Figures/Figure_1.png'],'Resolution',600)