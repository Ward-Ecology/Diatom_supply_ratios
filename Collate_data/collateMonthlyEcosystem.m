function EcosystemData = collateMonthlyEcosystem(ncFiles)
    % Check input
    if ~iscell(ncFiles) || isempty(ncFiles)
        error("Input must be a cell array of .nc file names.");
    end
    
    % Initialize an empty structure
    EcosystemData = struct();

    % Loop through each .mat file
    for i = 1:numel(ncFiles)
        file = ncFiles{i};

        if ~isfile(file)
            warning("File not found: %s. Skipping...", file);
        end

        % get netcdf info
        nci = ncinfo(file);

        npop = 0;
        
        for itr=1:numel(nci.Variables)

            % get variable names
            longname{itr}  = nci.Variables(itr).Attributes(1).Value;
            shortname{itr} = nci.Variables(itr).Name;
            
            % set single whitespaces to underscores
            longname{itr} = strrep(longname{itr},' ','_');
            
            % remove 'concentration' suffix
            longname{itr} = strrep(longname{itr},'_concentration','');
            
            % get variable units
            if numel(nci.Variables(itr).Attributes)>1
                varunits{itr} = nci.Variables(itr).Attributes(2).Value;
            else
                varunits{itr} = '-';
            end

            % Change names of population tracers
            if strcmp(longname{itr}(1),'c')
                % get population number
                pop_ID = str2num(longname{itr}(2:end));
                % change name 
                longname{itr} = strrep(longname{itr},'c','Population_');
            
                % Assign PFTs (hard-coded using Follet 2020 PFTs)
                if ismember(pop_ID,01:02); PFT{pop_ID} = 'Prokaryote';      end
                if ismember(pop_ID,03:04); PFT{pop_ID} = 'Picoeukaryote';   end
                if ismember(pop_ID,05:09); PFT{pop_ID} = 'Coccolithophore'; end
                if ismember(pop_ID,10:14); PFT{pop_ID} = 'Diazotroph';      end
                if ismember(pop_ID,15:23); PFT{pop_ID} = 'Diatom';          end
                if ismember(pop_ID,24:31); PFT{pop_ID} = 'Mixotroph';       end
                if ismember(pop_ID,32:47); PFT{pop_ID} = 'Zooplankton';     end
                if ismember(pop_ID,48:50); PFT{pop_ID} = 'Bacteria';        end
            end

            % for each tracer
            TRACER_id = shortname{itr};
            data = double(ncread(file,TRACER_id));
            if ismember(TRACER_id,{'X','Y'}) 
                % if time-independent metadata (lat/lon)
                EcosystemData.(TRACER_id) = data;
            elseif ismember(TRACER_id,{'T','iter'}) 
                % if time-dependent metadata
                EcosystemData.(longname{itr})(i) = data;
            elseif contains(TRACER_id,{'TRAC'}) 
                % if time-dependent data
                if contains(longname{itr},'Population')
                    % if population
                    if isfield(EcosystemData,'TotalBiomass')==0
                        EcosystemData.TotalBiomass  = zeros([size(data) numel(ncFiles)]);
                    end
                    if isfield(EcosystemData,PFT{pop_ID})==0
                        EcosystemData.(PFT{pop_ID}) = zeros([size(data) numel(ncFiles)]);
                    end
                    EcosystemData.TotalBiomass(:,:,:,i)  = EcosystemData.TotalBiomass(:,:,:,i) + data;
                    EcosystemData.(PFT{pop_ID})(:,:,:,i) = EcosystemData.(PFT{pop_ID})(:,:,:,i) + data;
                else
                    % if other BGC tracer
                    EcosystemData.(longname{itr})(:,:,:,i) = data;
                end
            end
        end
    end

    % Calculate total bioavailable nitrogen
    EcosystemData.TN  = EcosystemData.NO3 ...
                      + EcosystemData.NH4 ...
                      + EcosystemData.NO2 ...
                      + EcosystemData.DON ...
                      + EcosystemData.PON ...
                      + EcosystemData.Prokaryote      .* 0.1333 ...
                      + EcosystemData.Picoeukaryote   .* 0.1333 ...
                      + EcosystemData.Coccolithophore .* 0.1333 ...
                      + EcosystemData.Diazotroph      .* 0.1333 ...
                      + EcosystemData.Diatom          .* 0.1333 ...
                      + EcosystemData.Mixotroph       .* 0.1333 ...
                      + EcosystemData.Zooplankton     .* 0.1333 ...
                      + EcosystemData.Bacteria        .* 0.1333 ;

    % Calculate total bioavailable phosphate
    EcosystemData.TP  = EcosystemData.PO4 ...
                      + EcosystemData.DOP ...
                      + EcosystemData.POP ...
                      + EcosystemData.Prokaryote      .* 8.333e-3 ...
                      + EcosystemData.Picoeukaryote   .* 8.333e-3 ...
                      + EcosystemData.Coccolithophore .* 8.333e-3 ...
                      + EcosystemData.Diazotroph      .* 8.333e-3 ...
                      + EcosystemData.Diatom          .* 8.333e-3 ...
                      + EcosystemData.Mixotroph       .* 8.333e-3 ...
                      + EcosystemData.Zooplankton     .* 8.333e-3 ...
                      + EcosystemData.Bacteria        .* 8.333e-3 ;

    % Calculate total bioavailable iron
    EcosystemData.TFe = EcosystemData.FeT ...
                      + EcosystemData.DOFe ...
                      + EcosystemData.POFe ...
                      + EcosystemData.Prokaryote      .* 8.33e-6 ...
                      + EcosystemData.Picoeukaryote   .* 8.33e-6 ...
                      + EcosystemData.Coccolithophore .* 8.33e-6 ...
                      + EcosystemData.Diazotroph      .* 3.30e-5 ...
                      + EcosystemData.Diatom          .* 8.33e-6 ...
                      + EcosystemData.Mixotroph       .* 8.33e-6 ...
                      + EcosystemData.Zooplankton     .* 8.33e-6 ...
                      + EcosystemData.Bacteria        .* 8.33e-6 ;

    % Calculate total bioavailable silica
    EcosystemData.TSi = EcosystemData.SiO2 ...
                      + EcosystemData.POSi ...
                      + EcosystemData.Diatom          .* 0.4000 ;

    % Compute mean across months
    varNames = fieldnames(EcosystemData);
    for j = 1:numel(varNames)
        varName = varNames{j};
        nd = ndims(EcosystemData.(varName));
        if nd>=4
            EcosystemData.(varName) = struct('Monthly',     EcosystemData.(varName), ...
                                             'Annual', mean(EcosystemData.(varName), nd, 'omitnan'));
        end
    end

end