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

%% Get model data at Horstmann locations

NO3   = EcosystemData.NO3.Annual(:,:,1);
SiO2  = EcosystemData.SiO2.Annual(:,:,1);
% Diatom fraction of photosynthetic community (C biomass)
Dfrac = EcosystemData.Diatom.Annual(:,:,1) ...
     ./ ( EcosystemData.Prokaryote.Annual(:,:,1) ...
         +EcosystemData.Picoeukaryote.Annual(:,:,1) ...
         +EcosystemData.Coccolithophore.Annual(:,:,1) ...
         +EcosystemData.Diazotroph.Annual(:,:,1) ...
         +EcosystemData.Diatom.Annual(:,:,1) ...
         +EcosystemData.Mixotroph.Annual(:,:,1) );

% start and end points of Horstmann transect
x = [mean(HPLC.pigments.Longitude) mean(HPLC.pigments.Longitude)];
y = [  min(HPLC.pigments.Latitude)   max(HPLC.pigments.Latitude)];

% correct x > 180 to negative east
x(x>180) = x(x>180)-360;

% interpolate n points along transect
n = 1000;
x = linspace(x(1),x(2),n);
y = linspace(y(1),y(2),n);

[LonG, LatG] = meshgrid(EcosystemData.X, EcosystemData.Y);

NO3_vals  = interp2(LonG, LatG, NO3',x,y,'linear');
SiO2_vals = interp2(LonG, LatG, SiO2',x,y,'linear');
Dfrac_vals = interp2(LonG, LatG, Dfrac',x,y,'linear');



%% plot Horstmann et al. HPLC and nutrient data 
f4 = figure(4);
f4.Position = [74 676 887 594];
clf
mrkrsz = 50;

cmap = colororder;

t = tiledlayout(2,1, 'TileSpacing','compact', 'Padding','compact');

NO3  = HPLC.pigments.nitrate;
SiO2 = HPLC.pigments.silicate;
Si_star = SiO2 - NO3;

nexttile
plot(-y,NO3_vals,'color',cmap(2,:),'LineW',2),hold on
plot(-y,SiO2_vals,'color',cmap(1,:),'LineW',2)
plot(-y,SiO2_vals-NO3_vals,'color',cmap(3,:),'LineW',2)

scatter(-HPLC.pigments.Latitude,NO3,mrkrsz,'filled',...
    'MarkerFaceColor',cmap(2,:),'MarkerEdgeColor', 'k')
scatter(-HPLC.pigments.Latitude,SiO2,mrkrsz,'filled',...
    'MarkerFaceColor',cmap(1,:),'MarkerEdgeColor', 'k')
scatter(-HPLC.pigments.Latitude,Si_star,mrkrsz,'filled',...
    'MarkerFaceColor',cmap(3,:),'MarkerEdgeColor', 'k')

legend('N','Si','Si^*', 'AutoUpdate', 'off')
legend boxoff
set(gca,'XDir','reverse')
% xlabel('Latitude (^\circC)')
ylabel('Concentration (µM)')
axis tight
ylim([-20 40])
box on
yline(0)
title('(a) Surface nutrients and Si^*')
set(gca,'TitleHorizontalAlignment','left');


Diatom_frac = HPLC.PFT_fractions.("Diatom fraction");

nexttile
plot(-y,Dfrac_vals,'color',cmap(1,:),'LineW',2),hold on
scatter(-HPLC.PFT_fractions.Latitude,Diatom_frac,mrkrsz,'filled',...
    'MarkerFaceColor',cmap(1,:),'MarkerEdgeColor', 'k')
set(gca,'XDir','reverse')
xlabel('Latitude (^\circC)')
ylabel('Diatom Fraction')
axis tight
ylim([0 1])
box on
title('(b) Diatom biomass fraction')
set(gca,'TitleHorizontalAlignment','left');

%%


set(findall(gcf,'-property','FontSize'),'FontSize',14)


exportgraphics(gcf,['../Figures/Figure_4.png'],'Resolution',600)



