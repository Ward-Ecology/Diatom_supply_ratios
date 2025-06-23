function PhysicsData = collateMonthlyPhysics(matFiles)
    % Check input
    if ~iscell(matFiles) || isempty(matFiles)
        error("Input must be a cell array of .mat file names.");
    end
    
    % Initialize an empty structure
    PhysicsData = struct();

    % Loop through each .mat file
    for i = 1:numel(matFiles)
        file = matFiles{i};

        if ~isfile(file)
            warning("File not found: %s. Skipping...", file);
        end
        
        % Load variables from the file
        fileData = load(file);
        varNames = fieldnames(fileData);

        % Process each variable in the file
        for j = 1:numel(varNames)
            varName = varNames{j};
            data = fileData.(varName);

            % Store monthly data in a structure field
            if ~isfield(PhysicsData, varName)
                PhysicsData.(varName) = [];
            end

            % Append current month's data
            PhysicsData.(varName)(:,:,:,i) = data;
        end
    end

    % Compute mean across months
    varNames = fieldnames(PhysicsData);
    for j = 1:numel(varNames)
        varName = varNames{j};
        PhysicsData.(varName) = struct('Monthly', PhysicsData.(varName), ...
                                       'Annual', mean(PhysicsData.(varName), 4, 'omitnan'));
    end
end