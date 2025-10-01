clear
clc

fh = figure(2);                 % Create figure
fh.Position = [1 88 1084 1249];
% Set figure size
clf      

run('../Theoretical_model/RRT_diatoms.m')

% Set font size for all elements
set(findall(gcf,'-property','FontSize'),'FontSize',14) 
% Set figure background color to white
set(gcf,'Color','w')           

% Export figure to file

exportgraphics(gcf,['../Figures/Figure_2.png'],'Resolution',600)

return
%% Plot transects

% grid the supply ratios
[xx yy] = meshgrid(phiSiFe,phiSiN);

% find points in grid where phiSiFe == min(phiSiFe)
ii = find(xx==min(phiSiFe))';
% find points in grid where phiSiN == min(phiSiN)
jj = flipud(find(yy==min(phiSiN )))';

ij = unique([ii jj],'stable');

%%
BD_t = BD(ij);
BP_t = BP(ij);
Fe_t = Fe(ij);
 N_t =  N(ij);
Si_t = Si(ij);


figure(2)
clf
plot(BD_t./(BD_t+BP_t),'LineW',2)
hold on
plot(N_t,'LineW',2)
plot(Si_t,'LineW',2)
plot(Fe_t,'LineW',2)
% plot(Si_t-N_t,'LineW',2)
set(gca,'YScale','log')
        

legend('Diatom fraction','N','Si','Fe','Si^*',...
       'Location','SouthWest',...
       'AutoUpdate', 'off');

[~, idx1] = min(abs(yy(ij) - parameters(2).S2N));
xline(idx1)
[~, idx2] = min(abs(yy(ij)./xx(ij) - parameters(2).F2N));
xline(idx2)
[~, idx3] = min(abs(xx(ij) - parameters(2).S2N./parameters(2).F2N));
xline(idx3)

set(gca,'YTick',10.^[-4:2])
% set(gca,'XTick',[1 idx1 idx2 idx3 numel(ij)],...
%         'XTickLabel',{' ',...
%                       '\theta_{Si:N,D}',...
%                       '\theta_{N:Fe,D}',...
%                       '\theta_{Si:Fe,D}',...
%                       ''})
set(gca,'XTick',[])


xlim([1 numel(ij)])
ylim([3e-5 2e2])

cmap = colororder;
colororder(cmap([1 3 4 2],:))

% Set font size for all elements
set(findall(gcf,'-property','FontSize'),'FontSize',14) 
exportgraphics(gcf,['/Users/baw1d17/Library/CloudStorage/OneDrive-UniversityofSouthampton/Meetings and Presentations/2025 Traits VII/Transect.png'],...
                   'Resolution',600)

