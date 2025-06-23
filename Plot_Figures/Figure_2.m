clear
clc

fh = figure(2);                 % Create figure
fh.Position = [1 88 1084 1249]; % Set figure size
clf      

run('../Theoretical_model/RRT_diatoms.m')

% Set font size for all elements
set(findall(gcf,'-property','FontSize'),'FontSize',14) 
% Set figure background color to white
set(gcf,'Color','w')           

% Export figure to file

exportgraphics(gcf,['../Figures/Figure_2.png'],'Resolution',600)