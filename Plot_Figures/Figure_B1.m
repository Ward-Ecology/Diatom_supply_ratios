clear 
clc

% sim_name = 'RUN33_6_NEWEDES_HAM';
sim_name = 'RUN33_6_NEWEDES';

% data path 
data_path = '../Data';

%% load data

load([data_path '/' sim_name '_Ecosystem.mat'])
load([data_path '/' sim_name '_Physics.mat'])

Lat = EcosystemData.Y;
Lon = linspace(0,360,numel(EcosystemData.X));

phiSiN_proxy  = EcosystemData.phiSiN.Annual;
phiSiFe_proxy = EcosystemData.phiSiFe.Annual;
phiNFe_proxy  = EcosystemData.phiNFe.Annual;

%%
load([data_path '/explicit_fluxes.mat'])

 N_supply = double(mean(cat(3, total_flux{:,1}), 3));
Fe_supply = double(mean(cat(3, total_flux{:,3}), 3));
Si_supply = double(mean(cat(3, total_flux{:,4}), 3));

phiSiN_explicit  = Si_supply./ N_supply./3;
phiSiFe_explicit = Si_supply./Fe_supply./48000;
phiNFe_explicit  =  N_supply./Fe_supply./16000;


%%

fB1 = figure(9);
fB1.Position = [50 940 1408 397];
clf

t = tiledlayout(1,3, 'TileSpacing','compact', ...
                     'Padding','compact');

phi_explicit{1} = phiSiN_explicit;
phi_proxy{1}    = phiSiN_proxy;
phi_explicit{2} = phiSiFe_explicit;
phi_proxy{2}    = phiSiFe_proxy;
phi_explicit{3} = phiNFe_explicit;
phi_proxy{3}    = phiNFe_proxy;

lab_suffix = {'_{Si:N}','_{Si:Fe}','_{N:Fe}'};

clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:3
    nexttile;
    plot([1e-3 1e3],[1e-3 1e3],'r','LineW',1)
    hold on
    scatter(phi_explicit{i},phi_proxy{i},5,'k','filled', ...
                                        'MarkerFaceAlpha', 0.05, ...  % 30% opaque face
                                        'MarkerEdgeAlpha', 0.05)
    set(gca,'XScale','Log','YScale','Log')
    box on
    axis([1e-3 1e3 1e-3 1e3])
    axis square
    x= phi_explicit{i}(find(PhysicsData.HFacC(:,:,1)==1));
    y=    phi_proxy{i}(find(PhysicsData.HFacC(:,:,1)==1));

    x_ = x(~isnan(x) & ~isnan(y));
    y_ = y(~isnan(x) & ~isnan(y));
    x=x_;
    y=y_;

    accuracy = 100 .* (numel(find(x>1 & y>1 )) + numel(find(x<1 & y<1 ))) ./ numel(x);

    r=corrcoef(log10(x),log10(y), 'Rows', 'complete');
    title({['r^2 = ' num2str(r(2).^2,'%.2f')],...
           ['accuracy = ' num2str(accuracy,'%.1f') '%']})
    xlabel(['Explicit \phi' lab_suffix{i}])
    ylabel(['Proxy \phi' lab_suffix{i}])
    set(gca,'XTick',10.^(-3:3),'YTick',10.^(-3:3))
    xline(1)
    yline(1)
end

%%

set(findall(gcf,'-property','FontSize'),'FontSize',14)
exportgraphics(gcf,['../Figures/Figure_B1.png'],'Resolution',600)


