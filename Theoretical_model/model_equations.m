function [dVdt] = RCT_equations(t,variables,parameters)

%% retrieve model variables from structural array

B  = variables(1:end-3)'; % N biomass of populations
N  = variables(end-2);    % dissolved nitrogen
S  = variables(end-1);    % dissolved silicate
F  = variables(end  );    % dissolved iron

%% retrieve model parameters from structural array

% maximum growth rates
mumax = [parameters(1).mumax ... % maximum growth rate of generic phytoplankton
         parameters(2).mumax]  ; % maximum growth rate of diatoms

% N half-saturation concentrations
kN = [parameters(1).kN ... % N  half-saturation concentration of generic phytoplankton
      parameters(2).kN]  ; % N  half-saturation concentration of diatoms

% Fe half-saturation concentrations
kF = [parameters(1).kF ... % Fe half-saturation concentration of generic phytoplankton
      parameters(2).kF]  ; % Fe half-saturation concentration of diatoms

% Si half-saturation concentrations
kS = [parameters(1).kS ... % Si half-saturation concentration of generic phytoplankton
      parameters(2).kS]  ; % Si half-saturation concentration of diatoms

% Mortality rates
m  = [parameters(1).m ... % Mortality rates of generic phytoplankton
      parameters(2).m]  ; % Mortality rates of diatoms

% Fe:N stoichiometries
F2N = [parameters(1).F2N ... % Fe:N stoichiometries of generic phytoplankton
       parameters(2).F2N]  ; % Fe:N stoichiometries of diatoms

% Fe:N stoichiometries
S2N = [parameters(1).S2N ... % Si:N stoichiometries of generic phytoplankton
       parameters(2).S2N]  ; % Si:N stoichiometries of diatoms

% environmental parameters
kappa = parameters(3).kappa; % chemostat dilution rate
N0    = parameters(3).N0;    % Deep-water N
F0    = parameters(3).F0;    % Deep-water Fe
S0    = parameters(3).S0;    % Deep-water Si
f_atm = parameters(3).f_atm; % Iron deposition rate


%% auxiliary variables 

% nutrient limitation factors
gamma = min([ N./(N+kN) ; S./(S+kS) ; F./(F+kF) ]);

% population growth rates
growth = mumax .* gamma .* B; 

% population mortality rates
mortal = m .* B; % linear mortality rates

%% plankton PDEs

% generic phytoplankton (mmol N m^-3 d^-1)
dB_dt = + growth ... % growth
        - mortal   ; % mortality
 
%% nutrient PDEs

% nitrogen (mmol N m^-3 d^-1)
dN_dt = - sum(growth)    ... % uptake by all phytoplankton
        + kappa .* (N0 - N); % mixing

% Silicate (mmol Si m^-3 d^-1)
dS_dt = - sum(growth.*S2N) ... % uptake by all phytoplankton
        + kappa .* (S0 - S)  ; % mixing

% Iron (mmol Fe m^-3 d^-1)
dF_dt = - sum(growth.*F2N)  ... % uptake by all phytoplankton
        + kappa .* (F0 - F) ... % mixing
        + f_atm               ; % iron deposition

%% collate rates of change in output vector

dVdt = [dB_dt';
        dN_dt;
        dS_dt;
        dF_dt];
end