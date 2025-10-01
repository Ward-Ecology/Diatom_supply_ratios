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

%% Map WOA data onto MIT grid

% pad WOA data with duplicate arrays either side
V_nitrate  = repmat(WOA.nitrate.surface',1,3);
V_silicate = repmat(WOA.silicate.surface',1,3);

% Pad longitudes (X)
X = linspace(min(WOA.nitrate.lon)-360.5,max(WOA.nitrate.lon)+360.5,size(V_nitrate,2));
Y = WOA.nitrate.lat;

% WOA grid
[X ,Y ] = meshgrid(X,Y);

% MIT grid
Xq = linspace(0,360,numel(EcosystemData.X));
Yq = EcosystemData.Y;
[Xq,Yq] = meshgrid(Xq,Yq);

% Interpolant WOA onto MIT grid
WOA_nitrate  = interp2(X,Y,V_nitrate,Xq,Yq);
WOA_silicate = interp2(X,Y,V_silicate,Xq,Yq);

WOA_sistar = WOA_silicate - WOA_nitrate;

%% Get MITgcm Si*

MIT_nitrate  = EcosystemData.NO3.Annual(:,:,1)';
MIT_silicate = EcosystemData.SiO2.Annual(:,:,1)';

MIT_sistar = MIT_silicate - MIT_nitrate;

%% calculate provinces in global model
Provinces = zeros(size(EcosystemData.phiNFe.Annual)).*NaN;
Provinces(EcosystemData.phiNFe.Annual<1 & EcosystemData.phiSiN.Annual >1) = 1;
Provinces(EcosystemData.phiNFe.Annual<1 & EcosystemData.phiSiN.Annual <1) = 2;
Provinces(EcosystemData.phiNFe.Annual>1 & EcosystemData.phiSiFe.Annual<1) = 3;
Provinces(EcosystemData.phiNFe.Annual>1 & EcosystemData.phiSiFe.Annual>1) = 4;


%%
figure(10)
clf

plotdata{1} =    MIT_sistar   - WOA_sistar;
plotdata{2} =    MIT_silicate - WOA_silicate;
plotdata{3} = - (MIT_nitrate  - WOA_nitrate);

titles{1} = 'Surface Si^* error (MITgcm-WOA)';
titles{2} = 'Surface Silicate error (MITgcm-WOA)';
titles{3} = 'Negative Surface Nitrate error (MITgcm-WOA)';

t = tiledlayout(3,1, 'TileSpacing', 'compact', 'Padding', 'compact');

for i=1:3
    ax = nexttile;

    axh = axesm('MapProjection','mollweid',...
        'plinelocation',90,...
        'mlinelocation',Inf,...
        'MapLatLimit',[-80 80],...
        'MapLonLimit',[1 360]+130,...
        'Frame', 'on',...
        'Grid', 'off');


    pcolorm(Yq,Xq,plotdata{i})
    hold on
    contourm(Yq,Xq,Provinces',1.5:3.5,'k','LineW',1)

    clim = [-10:1:10];
    caxis([clim(1) clim(end)])
    colormap(redblue(numel(clim)-1))

    ch = colorbar;
    ch.Ticks = clim;

    box off
    landareas = shaperead('landareas.shp','UseGeoCoords',true);
    geoshow(landareas,...
        'FaceColor',[0.75 0.75 0.75],...
        'EdgeColor','k',...
        'LineWidth',0.5);
    
    title(titles{i})
    set(gcf,'Color','w')
    axis off


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
 end