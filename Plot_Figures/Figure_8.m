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

Sistar.MITgcm.surface(Sistar.MITgcm.surface==0)=NaN;
Sistar.MITgcm.z100(Sistar.MITgcm.z100==0)=NaN;

gamma_Fe = EcosystemData.gammaFe.Annual(:,:,1);
gamma_Si = EcosystemData.gammaSi.Annual(:,:,1);
gamma_N  = EcosystemData.gammaN.Annual(:,:,1);

%% Plot Figure 3a

f8 = figure(8);
clf
f8.Position = [74 437 1279 833];


t = tiledlayout(2,3, 'TileSpacing','compact',...
                     'Padding','compact');

%% MITgcm
nexttile
[~,i] = sort(gamma_N(:),'descend');
scatter(Sistar.MITgcm.z100(i),Sistar.MITgcm.surface(i),10,gamma_N(i),'filled')
caxis([0 1])
xline(0),yline(0),box on
hold on
axis([-20 80 -20 80])
axis square
colormap(winter)
xlabel('Si^* 100m (µM)')
ylabel('Si^* surface (µM)')
title('(a) \gamma_N (Model)')

nexttile
[~,i] = sort(gamma_Si(:),'descend');
scatter(Sistar.MITgcm.z100(i),Sistar.MITgcm.surface(i),10,gamma_Si(i),'filled')
caxis([0 1])
xline(0),yline(0),box on
axis([-20 80 -20 80])
axis square
colormap(winter)
xlabel('Si^* 100m (µM)')
ylabel('Si^* surface (µM)')
title('(b) \gamma_{Si} (Model)')

nexttile
[~,i] = sort(gamma_Fe(:),'descend');
scatter(Sistar.MITgcm.z100(i),Sistar.MITgcm.surface(i),10,gamma_Fe(i),'filled')
caxis([0 1])
xline(0),yline(0),box on
axis([-20 80 -20 80])
axis square
colormap(winter)
xlabel('Si^* 100m (µM)')
ylabel('Si^* surface (µM)')
title('(c) \gamma_{Fe} (Model)')

%% WOA
KNO3 = 0.31039;
KSi  = 0.93118;

gamma_N  = WOA.nitrate.surface  ./ (WOA.nitrate.surface  + KNO3);
gamma_Si = WOA.silicate.surface ./ (WOA.silicate.surface + KSi );

nexttile
[~,i] = sort(gamma_N(:),'descend');
scatter(Sistar.WOA.z100(i),Sistar.WOA.surface(i),10,gamma_N(i),'filled')
caxis([0 1])
xline(0),yline(0),box on
hold on
axis([-20 80 -20 80])
axis square
colormap(winter)
xlabel('Si^* 100m (µM)')
ylabel('Si^* surface (µM)')
title('(d) \gamma_N (WOA)')

nexttile
[~,i] = sort(gamma_Si(:),'descend');
scatter(Sistar.WOA.z100(i),Sistar.WOA.surface(i),10,gamma_Si(i),'filled')
caxis([0 1])
xline(0),yline(0),box on
axis([-20 80 -20 80])
axis square
colormap(winter)
xlabel('Si^* 100m (µM)')
ylabel('Si^* surface (µM)')
title('(e) \gamma_{Si} (WOA)')


ch = colorbar;

ch.Position = [0.68 0.0540 0.012 0.4118];
%%


set(findall(gcf,'-property','FontSize'),'FontSize',14)
exportgraphics(gcf,['../Figures/Figure_8.png'],'Resolution',600)
