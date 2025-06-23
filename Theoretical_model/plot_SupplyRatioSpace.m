function plot_SupplyRatioSpace(parameters,N0,S0,F0,pp,linecolor)

linewidth = 1;      % Set line width for plots
axscale = 3;        % Axis scaling factor (unused in current code)

nP = parameters(1).npopn;  % Number of non-diatom phytoplankton
nD = parameters(2).npopn;  % Number of diatom phytoplankton

isphytop = 1;       % Index for non-diatom phytoplankton
isdiatom = 2;       % Index for diatom phytoplankton

% Validate that there is exactly one of each phytoplankton type
if nP ~= 1 || nD ~= 1 || nD + nP ~= 2
    error('Can only handle 1 P and 1 D');
end

% Initialise vectors to store R* values
Rstar.N_phyto = [];
Rstar.S_phyto = [];
Rstar.F_phyto = [];

% Initialise vectors for stoichiometric ratios and mortality rates
thetaSiN = [];
thetaFeN = [];
m = [];

% Loop over phytoplankton types to compute R* values and traits
for i = 1:2
    % Calculate R* for each nutrient
    Rstar.N_phyto = [Rstar.N_phyto, parameters(i).kN * parameters(i).m / (parameters(i).mumax - parameters(i).m)];
    Rstar.S_phyto = [Rstar.S_phyto, parameters(i).kS * parameters(i).m / (parameters(i).mumax - parameters(i).m)];
    Rstar.F_phyto = [Rstar.F_phyto, parameters(i).kF * parameters(i).m / (parameters(i).mumax - parameters(i).m)];
    
    % Store stoichiometric ratios and mortality rates
    thetaSiN = [thetaSiN, parameters(i).S2N];
    thetaFeN = [thetaFeN, parameters(i).F2N];
    m = [m, parameters(i).m];
end

kappa = parameters(3).kappa;               % Nutrient uptake rate constant
thetaFeSi = thetaFeN ./ thetaSiN;          % Fe:Si stoichiometric ratio

n = 10000;  % Number of points for smooth curves

% Generate log-spaced vectors for nutrient concentrations
N0 = logspace(log10(min(N0)), log10(max(N0)), n);
S0 = logspace(log10(min(S0)), log10(max(S0)), n);
F0 = logspace(log10(min(F0)), log10(max(F0)), n);

% Ensure uniqueness
N0 = unique(N0);
S0 = unique(S0);
F0 = unique(F0);

hold on

% Calculate supply ratios
x = S0 ./ F0;  % Silicate to Iron ratio
y = S0 ./ N0;  % Silicate to Nitrogen ratio

%% Transition i - Diatom N to Si limitation transitions
% diatoms N/Si co-limited

% Retrieve R* values for diatoms
N = Rstar.N_phyto(isdiatom);
Si = Rstar.S_phyto(isdiatom);

% Compute transition functions
PhiSiN{1} = phiSiN_fun(thetaSiN(isdiatom), kappa, Si, N, N0);
PhiSiFe{1} = x;

%% Transition ii - Diatom Fe to Si limitation transitions
% diatoms F/Si co-limited

% Retrieve R* values for diatoms
Fe = Rstar.F_phyto(isdiatom);
Si = Rstar.S_phyto(isdiatom);

% Compute transition functions
PhiSiN{2} = y;
PhiSiFe{2} = phiSiFe_fun(thetaSiN(isdiatom),thetaFeN(isdiatom),kappa,Si,Fe,F0);

%% Transition iii - Diatom Fe to N limitation transitions
% diatoms F/N co-limited

% Retrieve R* values for diatoms
Fe = Rstar.F_phyto(isdiatom);
N  = Rstar.N_phyto(isdiatom);

% Compute transition functions
phiNFe = phiNFe_fun(thetaFeN(isdiatom),thetaFeN(isphytop),kappa,N,Fe,F0,m(isphytop),0);
PhiSiFe{3} = x;
PhiSiN{3}  = x./phiNFe;

%% Transition iv  - minimum N:Si required by non diatoms
% Si-limited diatoms and N-limited non-diatoms

% Retrieve R* values
N = Rstar.N_phyto(isphytop);
Si = Rstar.S_phyto(isdiatom);

% Compute transition functions
PhiSiN{4} = phiSiN_fun(thetaSiN(isdiatom),kappa,Si,N,N0);
PhiSiFe{4} = x;

%% Transition v  - minimum Fe:Si required by non diatoms
% Si-limited diatoms and Fe-limited non-diatoms

% Retrieve R* values
Fe = Rstar.F_phyto(isphytop);
Si = Rstar.S_phyto(isdiatom);

% Compute transition functions
PhiSiN{5} = y;
PhiSiFe{5} = phiSiFe_fun(thetaSiN(isdiatom),thetaFeN(isdiatom),kappa,Si,Fe,F0);

%% Transition vi - Non-diatom Fe to N limitation transition
% Si-limited diatoms, N/Fe co-limited non-diatoms
% Biomass matters

% Retrieve R* values
Fe = Rstar.F_phyto(isphytop);
N  = Rstar.N_phyto(isphytop);
Si = Rstar.S_phyto(isdiatom);

% get mortality rates
mD = m(isdiatom);
mP = m(isphytop);

% get biomasses
BD = kappa.*(S0-Si)./(thetaSiN(isdiatom).*mD);
BP = (kappa.*(N0-N) - mD.*BD) ./ mP;

% Compute transition functions
phiNFe = phiNFe_fun(thetaFeN(isdiatom),thetaFeN(isphytop),kappa,N,Fe,F0,m(isphytop),BP);

PhiSiN{6}  = x./phiNFe;
PhiSiFe{6} = x;

%% Transition vii - Negative Si* 
% Si-limited diatoms
% Biomass matters
Fe = Rstar.F_phyto(isdiatom);
 
% get mortality rates
mD = m(isdiatom);
mP = m(isphytop);

% get biomass
BD =  kappa.*(F0-Fe)./thetaFeN(isdiatom)./mD; % Fe limited diatoms

% Compute transition functions
PhiSiN{7} = (1 + (mD.*BD.*(1-thetaSiN(isdiatom))./kappa./S0)).^-1;
PhiSiFe{7} = x;

%%

% Plot diatom solutions
Fereplete = x < PhiSiFe{2};
plot(PhiSiFe{1}(Fereplete), PhiSiN{1}(Fereplete), '-', 'LineWidth', linewidth, 'Color', linecolor);

Nreplete = y < PhiSiN{1};
plot(PhiSiFe{2}(Nreplete), PhiSiN{2}(Nreplete), '-', 'LineWidth', linewidth, 'Color', linecolor);

Sireplete = x > PhiSiFe{2} & y > PhiSiN{1};
plot(PhiSiFe{3}(Sireplete), PhiSiN{3}(Sireplete), '-', 'LineWidth', linewidth, 'Color', linecolor);

% Plot non-diatom solutions
Fereplete = x < PhiSiFe{5};
plot(PhiSiFe{4}(Fereplete), PhiSiN{4}(Fereplete), '-', 'LineWidth', linewidth, 'Color', linecolor);

Nreplete = y < PhiSiN{4};
plot(PhiSiFe{5}(Nreplete), PhiSiN{5}(Nreplete), '-', 'LineWidth', linewidth, 'Color', linecolor);

Sireplete = x < PhiSiFe{5} & y < PhiSiN{4};
plot(PhiSiFe{6}(Sireplete), PhiSiN{6}(Sireplete), '-', 'LineWidth', linewidth, 'Color', linecolor);

% Optional: Plot negative Si* transition if pp == 5
if pp == 5
    Felimited = x > PhiSiFe{2};
    plot(PhiSiFe{7}(Felimited), PhiSiN{7}(Felimited), ':', 'LineWidth', linewidth, 'Color', linecolor);
end



end

%% Local Functions
% Compute the Si:N transition function
function [phiSiN] = phiSiN_fun(thetaSiND,Kappa,Si,N,N0)
    phiSiN = thetaSiND ... 
           + Kappa.*(Si-N.*thetaSiND)./(Kappa.*N0);
end

% Compute the Si:Fe transition function
function [phiSiFe] = phiSiFe_fun(thetaSiND,thetaFeND,Kappa,Si,Fe,F0)
    phiSiFe = thetaSiND./thetaFeND ...
            + Kappa.*(Si-Fe.*thetaSiND./thetaFeND)./(Kappa.*F0);
end

% Compute the N:Fe transition function
function [phiNFe] = phiNFe_fun(thetaFeND,thetaFeNP,Kappa,N,Fe,F0,mP,BP)
    phiNFe = 1./thetaFeND ...
            + (Kappa.*(N-Fe./thetaFeND) + mP.*BP.*(1-thetaFeNP./thetaFeND))./(Kappa.*F0);
end
