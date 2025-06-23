%% set model parameters in structural array

parameters(3).npopn = sum([parameters(1).npopn parameters(2).npopn]); % total

% Names of PFTs
parameters(1).pft ='P';
parameters(2).pft ='D';

% N.B. linspace uses second value if npopn = 1

% plankton initial populations
parameters(1).B0 = 1e-6; % initial populations of generic phytoplankton
parameters(2).B0 = 1e-6; % initial populations of diatoms

% maximum growth rates
parameters(1).mumax = 1.0; % maximum growth rate of generic phytoplankton
parameters(2).mumax = 2.0; % maximum growth rate of diatoms

% N half-saturation concentrations
parameters(1).kN = 1.0e+0; % N  half-saturation concentration of generic phytoplankton
parameters(2).kN = 1.0e+0; % N  half-saturation concentration of diatoms

% Si half-saturation concentrations
parameters(1).kS = 0.0e+0; % Si half-saturation concentration of generic phytoplankton
parameters(2).kS = 1.0e+0; % Si half-saturation concentration of diatoms

% Fe half-saturation concentrations
parameters(1).kF = 1.0e-3; % Fe  half-saturation concentration of generic phytoplankton
parameters(2).kF = 1.0e-3; % Fe  half-saturation concentration of diatoms

% Fe:N stoichiometries
parameters(1).F2N = 1.0e-3; % Fe:N stoichiometries of generic phytoplankton
parameters(2).F2N = 1.0e-3; % Fe:N stoichiometries of diatoms

% Si:N stoichiometries
parameters(1).S2N = 0.0e+0; % Si:N stoichiometries of generic phytoplankton
parameters(2).S2N = 1.0e+0; % Si:N stoichiometries of diatoms


% linear mortality rates
parameters(1).m = linspace(1.0e-1,1.0e-1,parameters(1).npopn); % mortality rate of generic phytoplankton
parameters(2).m = linspace(0.95e-1,0.95e-1,parameters(2).npopn); % mortality rate of diatoms

% environmental parameters
parameters(3).kappa = 0.1;    % chemostat dilution rate
parameters(3).N0    = 10;    % Deep-water N
parameters(3).F0    = 1.0e-4; % Deep-water Fe
parameters(3).S0    = 10;    % Deep-water Si
parameters(3).f_atm = 0.0;    % Iron deposition rate