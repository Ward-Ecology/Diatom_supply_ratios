clear
clc
%% read Horstman data

data_path = '../Original_data/Horstmann/';

% % nutrient data
% T = readtable('MD206_TS1.tab','FileType','text','VariableNamingRule','preserve');
% Tmean = varfun(@mean,T,'GroupingVariables','Event','InputVariables', @isnumeric);

% pigment data
HPLC.pigments = readtable([data_path 'MD206_pigments.tab'],'FileType','text','HeaderLines',51,'VariableNamingRule','preserve');
environmental_data  = readtable([data_path 'MD206_TS1.tab'],'FileType','text','HeaderLines',48,'VariableNamingRule','preserve');

% Option 1
% filter by timeslice
idx = HPLC.pigments.Timeslice == 0;
HPLC.pigments = HPLC.pigments(idx,:);

% Option 2
% % mean data by event
% Pmean = varfun(@mean,P,'GroupingVariables','Event','InputVariables', @isnumeric);
% for i = 1:length(Pmean.Event)
%     Pmean.Event(i) = strrep(Pmean.Event(i),'MD206_','');
% end

% % sort data as Horstman et al.
% Pmean = Pmean([10 2 1 15 13 14 12 3 4 5 11 9 6 8 7],:);


clc
for i=1:size(HPLC.pigments,1)
    ii = find(matches(environmental_data.(1),HPLC.pigments.(1)(i)));        

    % temperature data
    if ~isempty(environmental_data.(8)(ii))
        HPLC.pigments.temperature(i) = unique(environmental_data.(8)(ii));
    else
        HPLC.pigments.temperature(i) = NaN;
    end      

    % nitrate data
    if ~isempty(environmental_data.(8)(ii))
        HPLC.pigments.nitrate(i) = mean(environmental_data.(11)(ii));
    else
        HPLC.pigments.nitrate(i) = NaN;
    end      

    % silicate data
    if ~isempty(environmental_data.(8)(ii))
        HPLC.pigments.silicate(i) = mean(environmental_data.(15)(ii));
    else
        HPLC.pigments.silicate(i) = NaN;
    end

    
end

clear i idx ii H
%% Assign variables

% assign pigments
Fuco     = HPLC.pigments.("Fuco [µg/l]");
Perid    = HPLC.pigments.("Perid [µg/l]");
Hex      = HPLC.pigments.("Hex-fuco [µg/l]");
Allo     = HPLC.pigments.("Allo [µg/l]");
But      = HPLC.pigments.("But-fuco [µg/l]");
Chl_b    = HPLC.pigments.("Chl b [µg/l]");
Zea      = HPLC.pigments.("Zea [µg/l]");
Chl_a    = HPLC.pigments.("Chl a [µg/l]");
DV_chl_a = HPLC.pigments.("DV chl a [µg/l]");

% diagnostic pigment concentrations and weights 
% Vectorised to align with Brewin 2010
W = [1.41 1.41  1.27 0.6  0.35 1.01  0.86];
P = [Fuco Perid Hex  Allo But  Chl_b Zea ];

% weighted pigments
WP = W.*P;
% total chlorophyll concentration
C = sum(WP,2);

% Uitz et al. assign all hex to nanoplankton
% but Brewin (2010) adjusted this 
% so that more hex assigned to picoplankton 
% in ultra-oligotrophic conditions
% F_pico_hex is fraction of hex assigned to pico
% Xn = 1-Yp is fraction assigned to nano
% if C <= 0.001, Yp = 1 Xn = 0
% if C >= 0.080, Yp = 0 Xn = 1
% linear transition for 0.001<C<0.08

F_pico_hex = (C-0.08)./(0.001-0.08);
F_pico_hex = min(max(F_pico_hex,0),1);
F_nano_hex = 1 - F_pico_hex;

% Only stations 16 and 18 qualify as ultra-oligotrophic, 
% so in most cases, all hex should go to nano

% Note that it looks like Hortsmann's analysis 
% assigned hex to pico in most cases,
% massively over-estimating pico biomassover-estimated

%% Hirata PFTs 
microplankton   = 1.41.*(Fuco + Perid)./C;
nanoplankton    = (F_nano_hex.*1.27.*Hex + 0.60.*Allo + 0.35.*But + 1.01.*Chl_b)./C;
picoplankton    = (F_pico_hex.*1.27.*Hex + 0.86.*Zea)./C;

diatoms         = 1.41.*Fuco./C;
dinoflagellates = 1.41.*Perid./C;
green_algae     = 1.01.*Chl_b./C;
prokaryotes     = 0.86.*Zea./C;
haptophytes     = (F_nano_hex.*1.27.*Hex + 0.60.*Allo + 0.35.*But)./C;
pico_eukaryotes = (F_pico_hex.*1.27.*Hex)./C;

% Prochlorococcus excluded as it appears to be a subset of prokaryotes
% prochlorococcus = 0.74.*DV_chl_a./Chl_a;

% make table object for size fcations
HPLC.size_fractions = table(HPLC.pigments.Latitude, ...
                            HPLC.pigments.Longitude, ...
                            C, ...
                            picoplankton, ...
                            nanoplankton, ...
                            microplankton, ...
                            'VariableNames', {'Latitude', ...
                                              'Longitude', ...
                                              'Total Chlorophyll', ...
                                              'Picoplankton fraction', ...
                                              'Nanoplankton fraction', ...
                                              'Microplankton fraction'});


% make table object for PFTs
HPLC.PFT_fractions = table(HPLC.pigments.Latitude, ...
                           HPLC.pigments.Longitude, ...
                           C, ...
                           prokaryotes, ...
                           pico_eukaryotes, ...
                           green_algae, ...
                           haptophytes, ...
                           dinoflagellates, ...
                           diatoms, ...
                           'VariableNames', {'Latitude', ...
                                             'Longitude', ...
                                             'Total Chlorophyll', ...
                                             'Prokaryote fraction', ...
                                             'Pico eukaryote fraction', ...
                                             'Green algae fraction', ...
                                             'Haptophyte fraction', ...
                                             'Dinoflagellate fraction', ...
                                             'Diatom fraction'});



clear data_path C DV_chl_a environmental_data W WP P F_nano_hex F_pico_hex Allo but Chl_a Chl_b flip_correction Fuco Hex Perid Zea
clear microplankton diatoms dinoflagellates nanoplankton green_algae haptophytes picoplankton prokaryotes pico_eukaryotes prochlorococcus 

save('../Data/Horstmann.mat','HPLC')