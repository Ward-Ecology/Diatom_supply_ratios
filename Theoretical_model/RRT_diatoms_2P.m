


% number of plankton populations
parameters(1).npopn = 2;        % Set number of non-diatom phytoplankton populations
parameters(2).npopn = 1;        % Set number of diatom phytoplankton populations

%% initialise model parameters
RRT_diatom_parameters           % Load default model parameters for diatom and non-diatom species

% non-default parameters (other parameters set in (RRT_diatom_parameters))
parameters(1).kN  = [0.10 1.00];       % N  half-saturation constant - non-diatoms
parameters(1).kS  = [0.00 0.00];       % Si half-saturation constant - non-diatoms (0 = not limiting)
parameters(1).kF  = [1e-4 1e-3];       % Fe half-saturation constant - non-diatoms
parameters(2).S2N = 3.00;       % Si:N ratio - diatoms
parameters(2).kN  = 1.00;       % N  half-saturation constant - diatoms
parameters(2).kS  = 5.00;       % Si half-saturation constant - diatoms
parameters(2).kF  = 1e-3;       % Fe half-saturation constant - diatoms

% simulation length (days)
tspanmax = 500;

% resolution of solution grid (n by n)
n = 200; 

% ode solver options
options = odeset('RelTol',1e-6,...
                 'AbsTol',1e-6);  % Set relative and absolute tolerance for ODE solver
tspan = [0 tspanmax];          % Time span for ODE simulation

% Define colormap for plotting using color hex codes
% colormap % https://coolors.co/palettes/popular/5%20colors
clrs = [244 241 222 ;... % blank ...
        224 122  95 ;... % prov i   ('#3d405b')
        242 204 143 ;... % prov ii  ('#81b29a')
        129 178 154 ;... % prov iii ('#f2cc8f')
         61  64  91];... % prov iv  ('#e07a5f')
             
clrs=clrs./255;

% extract frequently used parameters
thetaSiD = parameters(2).S2N;  % Diatom Si:N ratio
thetaFeD = parameters(2).F2N;  % Diatom Fe:N ratio

% environmental parameters
S0 = 10;                       % Fixed initial silicate concentration
N0 = logspace(log10(0.1),log10(100),n); % Nitrogen gradient (log spaced)
F0 = N0.*mean(thetaFeD);       % Iron scaled to nitrogen using mean Fe:N ratio

phiFeSi = F0./S0;              % Iron to silicate ratio
phiNSi  = N0./S0;              % Nitrogen to silicate ratio

phiSiFe = 1./phiFeSi;          % Silicate to iron ratio
phiSiN  = 1./phiNSi;           % Silicate to nitrogen ratio

% Preallocate output arrays
k_transect = zeros(n,n,parameters(3).npopn+3);
R = zeros(n,n);

%% iterate along transect
for i=1:n                      % Loop over Fe:Si ratios
    for j=1:n                  % Loop over N:Si ratios
        clc,disp([i j])        % Display current indices

        parameters(3).F0 = F0(i); % Set initial iron concentration
        parameters(3).N0 = N0(j); % Set initial nitrogen concentration
        parameters(3).S0 = S0;    % Set initial silicate concentration

        % initialise state variables
        variables = [parameters(1).B0';         % Non-diatom biomass
                     parameters(2).B0';         % Diatom biomass
                     parameters(3).N0;          % N
                     parameters(3).S0;          % Si
                     parameters(3).F0];         % Fe

        % Run ODE solver
        [t,y] = ode45( ...
                      @(t,y) model_equations(t,y,parameters), ... 
                      tspan, ...
                      variables, ...
                      options ...
                      );

        % Save final state of simulation
        k_transect(j,i,:) = y(end,:); 
    end
end

%% Collate model outputs
nP = parameters(1).npopn;      % Number of non-diatom populations
nB = parameters(2).npopn;      % Number of diatom populations

BP = sum(k_transect(:,:,1:nP),3);        % Sum non-diatom biomass
BD = sum(k_transect(:,:,nP+1:nP+nB),3);  % Sum diatom biomass
N  = k_transect(:,:,end-2);              % Extract nitrogen values
Si = k_transect(:,:,end-1);              % Extract silicate values
Fe = k_transect(:,:,end);                % Extract iron values

ttls  ={'Provinces','N','Diatom fraction','Si','Si^*','Fe'}; % Subplot titles
unts  ={'',[char(181) 'M'],'',[char(181) 'M'],[char(181) 'M'],[char(181) 'M']}; % Units for colorbars

colormap(turbo(256))            % Set default colormap

                       % Clear figure

%% Province colour legend

subplot(321);
% Plot 4 polygon regions
% Province 1
ps = polyshape([max(phiSiFe) min(phiSiFe) min(phiSiFe) thetaSiD./thetaFeD max(phiSiFe)] ...
    ,[max(phiSiN) max(phiSiN) thetaSiD thetaSiD max(phiSiN)]); 
pg = plot(ps);
pg.FaceAlpha = 1;
pg.FaceColor = clrs(5,:);
hold on
% Province 2
ps = polyshape([thetaSiD./thetaFeD min(phiSiFe) min(phiSiFe) thetaSiD./thetaFeD],...
    [thetaSiD min(phiSiN) thetaSiD thetaSiD]);
pg = plot(ps);
pg.FaceAlpha = 1;
pg.FaceColor = clrs(4,:);
% Province 3
ps = polyshape([thetaSiD./thetaFeD min(phiSiFe) thetaSiD./thetaFeD thetaSiD./thetaFeD],...
    [thetaSiD min(phiSiN) min(phiSiN) thetaSiD]);
pg = plot(ps);
pg.FaceAlpha = 1;
pg.FaceColor = clrs(3,:);
% Province 4
ps = polyshape([max(phiSiFe) thetaSiD./thetaFeD thetaSiD./thetaFeD max(phiSiFe) max(phiSiFe)],...
    [max(phiSiN) thetaSiD min(phiSiN) min(phiSiN) max(phiSiN)]);
pg = plot(ps);
pg.FaceAlpha = 1;
pg.FaceColor = clrs(2,:);

% annotate provinces
tx = (3.3e-4.*[1 12 2.4 0.2]).^-1; % Annotation markers - x coordinates
ty = (0.06  .*[1 12 60 5]).^-1;    % Annotation markers - y coordinates
ph1 = plot(tx,ty,'w-','LineWidth',2);      % White lines
ph2 = plot(tx,ty,'wo');                    % White circle markers
ph2.MarkerFaceColor = 'w';
ph2.MarkerEdgeColor = 'k';
ph2.MarkerSize = 20;
ph2.LineWidth = 1;
th = text(tx,ty,{'i','ii','iii','iv'},'HorizontalAlignment','center'); % Label points

%% N concentration
subplot(322);
imagesc(phiSiFe,phiSiN,N)         
set(gca,'ColorScale','log') % Log colour scale
% annotate
text(0.009.^-1,0.6.^-1,{'N = R^*_{N,P}'},'color','k','rotation',0,'HorizontalAlignment','left')
text(0.0005.^-1,0.05.^-1,{'N = R^*_{N,D}'},'color','k','rotation',0,'HorizontalAlignment','center')

%% Diatom fraction
ax1 = subplot(323);
imagesc(phiSiFe,phiSiN,BD./(BP+BD))
set(gca,'ColorScale','lin') % Linear colour scale
clim([0 1]) % data limited 0 to 1
% annotate
text(5e-5.^-1,0.6.^-1,{'Diatoms win','(Fe limited)'},'color','w','rotation',0,'HorizontalAlignment','center')
text(0.0005.^-1,0.05.^-1,{'Diatoms win','(N limited)'},'color','w','rotation',0,'HorizontalAlignment','center')
text(0.0004.^-1,6.^-1,{'Coexistence','(Si/Fe limited)'},'color','w','rotation',0,'HorizontalAlignment','right')
text(0.009.^-1,0.6.^-1,{'Coexistence','(Si/N limited)'},'color','w','rotation',0,'HorizontalAlignment','left')
cmap = flipud(winter(256));
colormap(ax1,cmap)

%% Si concentration
subplot(324);
imagesc(phiSiFe,phiSiN,Si)
set(gca,'ColorScale','log') % Log colour scale
% annotate
text(0.009.^-1,0.6.^-1,{'Si = R^*_{Si,D}'},'color','k','rotation',0,'HorizontalAlignment','left')
text(0.0004.^-1,5.^-1,{'Si = R^*_{Si,D}'},'color','k','rotation',0,'HorizontalAlignment','right')

%% Si*
ax2 = subplot(325);
imagesc(phiSiFe,phiSiN,Si-N)
set(gca,'ColorScale','lin') % Linear colour scale
clim([-10 10])
% red-blue colour scale for Si*
cmap = redblue(256); 
colormap(ax2,cmap)

%% Fe concentration
subplot(326);
imagesc(phiSiFe,phiSiN,Fe)
set(gca,'ColorScale','log') % Log colour scale
% annotate
text(0.0004.^-1,5.^-1,{'Fe = R^*_{Fe,P}'},'color','k','rotation',0,'HorizontalAlignment','right')
text(5e-5.^-1,0.4.^-1,{'Fe = R^*_{Fe,D}'},'color','k','rotation',0,'HorizontalAlignment','center')


%% Axis configuration and annotations for all subplots
lnclr={'k','w','w','w','k','w'}; % Line colors for overlay plots
for i=1:6
    oldAx = subplot(3,2,i); % get original axes
    set(oldAx,'XDir','normal','YDir','normal') % set axis directions
    axis square
    % label axes
    xlabel('\phi_{Si:Fe}\rightarrow','Position',[min(phiSiFe).*5.0 min(phiSiN)./1.5]);
    ylabel('\phi_{Si:N}\rightarrow' ,'Position',[min(phiSiFe)./2.5 min(phiSiN).*4.0],'Rotation',90);
    % log x and y scales
    set(gca,'XScale','log','YScale','log')
    
    % call function to plot province boundaries
    plot_SupplyRatioSpace(parameters,N0,S0,F0,i,lnclr{i})

    % Place tick marks at transition ratios
    set(oldAx,'XTick',unique([min(phiSiFe) thetaSiD./thetaFeD max(phiSiFe)]),...
              'YTick',unique([min(phiSiN) 1 thetaSiD max(phiSiN)]),...
              'XTickLabelRotation',0)
    % label transitions
    oldAx.XTickLabel{2} = '\theta_{Si:Fe,D}';
    oldAx.YTickLabel{3} = '\theta_{Si:N,D}';
    % set axis limits
    axis([min(phiSiFe) max(phiSiFe) min(phiSiN) max(phiSiN) ])

    % new axes for inverse stoichiometries
    newAx = axes('Position', oldAx.Position);
    axis square
    % configure new axes
    set(newAx,...
        'color', 'none',... % blank so can see original 
        'XDir','reverse','YDir','reverse',... % reverse axes
        'XScale','log','XTick',fliplr(1./oldAx.XTick),... % log scale and ticks
        'YScale','log','YTick',fliplr(1./oldAx.YTick),.... % log scale and ticks
        'XLim', fliplr(1./oldAx.XLim), 'XAxisLocation', 'top',... % Flip limits and axis L/R
        'YLim', fliplr(1./oldAx.YLim), 'YAxisLocation', 'right'); % Flip limits and axis T/B
    % label transitions
    newAx.XTickLabel{2} = '\theta_{Fe:Si,D}';
    newAx.YTickLabel{2} = '\theta_{N:Si,D}';
    % label axes
    xlabel('\leftarrow\phi_{Fe:Si}','Position',[min(phiFeSi).*5 min(phiNSi)./1.5]);
    ylabel('\leftarrow\phi_{N:Si}' ,'Position',[min(phiFeSi)./1.5 min(phiNSi).*5],'Rotation',90);

    % generate subplot titles
    ci = ceil(i/2)+rem(i-1,2)*3;
    ttl = title(['(' char(96+ci) ') ' ttls{i}]);
    % adjust position
    ttl.Units = 'Normalize';
    ttl.Position(1:2) = [-0.25 1.15];
    ttl.HorizontalAlignment = 'left';

    % set color scale limits
    if strcmp(oldAx.ColorScale,'lin')
        clim(oldAx,[floor(oldAx.CLim(1)) ceil(oldAx.CLim(2))])
    elseif strcmp(oldAx.ColorScale,'log')
        clim(oldAx,10.^[floor(log10(oldAx.CLim(1))) ceil(log10(oldAx.CLim(2)))])
    end

    % generate and place colorbar
    if i~=1
        hcb = colorbar(newAx);
        hcb.Position(1) = hcb.Position(1) + 0.05;
        colormap(newAx,colormap(oldAx))
        clim(newAx,clim(oldAx))
        newAx.ColorScale = oldAx.ColorScale;
        set(newAx,'Position',oldAx.Position);
        set(oldAx,'Position',oldAx.Position);
        colorTitleHandle = get(hcb,'Title');
        set(colorTitleHandle ,'String',unts{i});
    end
end
%%
uistack(ph1,"top") % Bring line to top
uistack(ph2,"top") % Bring points to top
uistack(th,"top")  % Bring text to top

