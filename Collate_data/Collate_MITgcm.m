clear
clc

simulation_name = 'RUN33_6_NEWEDES';
% simulation_name = 'RUN33_6_NEWEDES_HAM';

%% collate physics data

% list monthly physics files
Monthly_Physics_Files = dir('../Original_data/MITgcm/Physics/phys*.mat');

% get paths and files for each file
matFiles = fullfile({Monthly_Physics_Files.folder}, {Monthly_Physics_Files.name});

% collate physics data
PhysicsData = collateMonthlyPhysics(matFiles);

%% collate ecosystem data

% list monthly darwin ecosystem files
Monthly_Ecosystem_Files = dir(['../Original_data/MITgcm/' simulation_name '/3d.00*.nc']);

% get paths and files for each file
matFiles = fullfile({Monthly_Ecosystem_Files.folder}, {Monthly_Ecosystem_Files.name});

% collate ecosystem data
EcosystemData = collateMonthlyEcosystem(matFiles);


%% load grid data
grid = load('../Original_data/MITgcm/Grid/grid_ecco2.mat');

PhysicsData.Depth = grid.Depth;
PhysicsData.HFacC = grid.HFacC;
PhysicsData.X     = grid.X;
PhysicsData.Y     = grid.Y;
PhysicsData.Z     =-grid.Z;
PhysicsData.Depth = grid.Depth;
PhysicsData.drC   = grid.drC;
PhysicsData.drF   = grid.drF;
PhysicsData.dxG   = grid.dxG;
PhysicsData.dyG   = grid.dyG;
PhysicsData.rA    = grid.rA;

EcosystemData.Z     =-grid.Z;
EcosystemData.Depth = grid.Depth;
EcosystemData.drF   = grid.drF;
EcosystemData.rA    = grid.rA;

%% Calculate normalised supply ratios (integrating across surface 100m)

% Integrate surface 100m
EcosystemData.TN.Monthly  = tensorprod(EcosystemData.TN.Monthly( :,:,1:7,:),EcosystemData.drF(1:7),3,1);
EcosystemData.TP.Monthly  = tensorprod(EcosystemData.TP.Monthly( :,:,1:7,:),EcosystemData.drF(1:7),3,1);
EcosystemData.TFe.Monthly = tensorprod(EcosystemData.TFe.Monthly(:,:,1:7,:),EcosystemData.drF(1:7),3,1);
EcosystemData.TSi.Monthly = tensorprod(EcosystemData.TSi.Monthly(:,:,1:7,:),EcosystemData.drF(1:7),3,1);
EcosystemData.TN.Annual  = tensorprod(EcosystemData.TN.Annual( :,:,1:7),EcosystemData.drF(1:7),3,1);
EcosystemData.TP.Annual  = tensorprod(EcosystemData.TP.Annual( :,:,1:7),EcosystemData.drF(1:7),3,1);
EcosystemData.TFe.Annual = tensorprod(EcosystemData.TFe.Annual(:,:,1:7),EcosystemData.drF(1:7),3,1);
EcosystemData.TSi.Annual = tensorprod(EcosystemData.TSi.Annual(:,:,1:7),EcosystemData.drF(1:7),3,1);

% Calculate normalised monthly supply ratios
EcosystemData.phiSiN.Monthly  = EcosystemData.TSi.Monthly ./ EcosystemData.TN.Monthly ...
                             ./ 3;
EcosystemData.phiSiFe.Monthly = EcosystemData.TSi.Monthly ./ EcosystemData.TFe.Monthly ...
                             ./ 48000;
EcosystemData.phiNFe.Monthly  = EcosystemData.TN.Monthly ./ EcosystemData.TFe.Monthly ...
                             ./ 16000;

% Calculate normalised annual supply ratios
EcosystemData.phiSiN.Annual  = EcosystemData.TSi.Annual ./ EcosystemData.TN.Annual ...
                            ./ 3;
EcosystemData.phiSiFe.Annual = EcosystemData.TSi.Annual ./ EcosystemData.TFe.Annual ...
                            ./ 48000;
EcosystemData.phiNFe.Annual  = EcosystemData.TN.Annual ./ EcosystemData.TFe.Annual ...
                            ./ 16000;

%% Calculate nutrient limitations

% half saturations 10 micron diatom
KNO3 = 0.31039;
KPO4 = 0.0194;
KSi  = 0.93118;
KFe  = 1.94e-05;

EcosystemData.gammaN.Monthly  = EcosystemData.NO3.Monthly  ./ (EcosystemData.NO3.Monthly  + KNO3);
EcosystemData.gammaN.Annual   = EcosystemData.NO3.Annual   ./ (EcosystemData.NO3.Annual   + KNO3);

EcosystemData.gammaP.Monthly  = EcosystemData.PO4.Monthly  ./ (EcosystemData.PO4.Monthly  + KPO4);
EcosystemData.gammaP.Annual   = EcosystemData.PO4.Annual   ./ (EcosystemData.PO4.Annual   + KPO4);

EcosystemData.gammaFe.Monthly = EcosystemData.FeT.Monthly  ./ (EcosystemData.FeT.Monthly  + KFe );
EcosystemData.gammaFe.Annual  = EcosystemData.FeT.Annual   ./ (EcosystemData.FeT.Annual   + KFe );

EcosystemData.gammaSi.Monthly = EcosystemData.SiO2.Monthly ./ (EcosystemData.SiO2.Monthly + KSi );
EcosystemData.gammaSi.Annual  = EcosystemData.SiO2.Annual  ./ (EcosystemData.SiO2.Annual  + KSi );

EcosystemData.gamma.Monthly  =  min(...
                                    cat(5,...
                                    EcosystemData.gammaN.Monthly ,...
                                    EcosystemData.gammaP.Monthly ,...
                                    EcosystemData.gammaFe.Monthly,...
                                    EcosystemData.gammaSi.Monthly...
                                    ),...
                                    [],5);

EcosystemData.gamma.Annual   =  min(...
                                    cat(5,...
                                    EcosystemData.gammaN.Annual ,...
                                    EcosystemData.gammaP.Annual ,...
                                    EcosystemData.gammaFe.Annual,...
                                    EcosystemData.gammaSi.Annual...
                                    ),...
                                    [],5);
%% save data

save(['../Data/' simulation_name '_Ecosystem.mat'],'EcosystemData','-v7.3')

save(['../Data/' simulation_name '_Physics.mat'],'PhysicsData','-v7.3')

