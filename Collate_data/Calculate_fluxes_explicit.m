clear
clc

run_id = 'RUN33_6_NEWEDES';

% find incoming fluxes to surface 100 m (6 layers)
n_layers = 6;

% scaling for quiver plots
sc = 0.25;

% load grid data
load('../Original_data/MITgcm/Grid/grid_ecco2.mat')

% dust data
fileID = fopen(['../Original_data/MITgcm/' run_id '/songmiao_mineral_mth-2d.bin']);
dust = fread(fileID,[1 Inf],'single','ieee-be');
dust = reshape(dust,360,160,12);

for i_month = 1:12

    % load physics data
    load(['../Original_data/MITgcm/Physics/phys_' num2str(i_month,'%02i') '.mat'])

    fname = ['../Original_data/MITgcm/' run_id '/3d.0000026400.nc'];
    % load nutrient data
    for itr = 2:7
        nutr(itr,:,:,:) = ncread(fname,['TRAC' num2str(itr,'%02i')]);
    end
    for itr = 1:3
        DOM{itr} = ncread(fname,['TRAC' num2str(itr+8,'%02i')]);
    end
    DOM{4} = DOM{3}.*0;
    for itr = 1:4
        POM{itr} = ncread(fname,['TRAC' num2str(itr+12,'%02i')]);
    end
    % load('run33_6_newedes_coccop_y10_car_nut_gr.mat')

    % ALSO NEED ATMOSPHERIC DEPOSITION

    % u, v and w are located at western, southern and lower boundaries of each grid box

    % let cx, cy and cz be coordinates at grid centres
    cx = X; % degrees
    cy = Y; % degrees
    cz = Z; % metres

    nx=length(X);
    ny=length(Y);
    nz=length(Z);

    % differences along uniform grid
    deltaX = mean(diff(cx));
    deltaY = mean(diff(cy));
    deltaZ = drF;

    % get face areas and cell volumes
    for i=1:23
        % area of E/W face = y distance * depth
        u_area(:,:,i) = dyG(2:end,:).*drF(i);
        % area of N/S face = x distance * depth
        v_area(:,:,i) = dxG(:,2:end).*drF(i);
        % cell volume
        volume(:,:,i) = rA.*drF(i);
    end
    % area of U\L face = rA
    w_area = repmat(rA,[1,1,23]);

    % re-mesh centres of grid
    [cxx,cyy,czz] = ndgrid(cx,cy,cz,3);

    % % find coordinates at edges of grid boxes
    % x = [cx - deltaX./2 ; cx(end)+deltaX./2];
    % y = [cy - deltaY./2 ; cy(end)+deltaY./2];
    % z = [0;-cumsum(drF)];

    % get coordinates at grid edges
    x = cx-deltaX./2;
    y = cy-deltaY./2;


    %%
    % titles, data limits, indices
    nut_ttls = {'DIN','PO_4','Fe','Si'};
    nut_clms = [1e-3 1e2];
    nut_inds = {2:4,5,6,7};
    N2X      = { 1, 16, 16000,  3};

    for ip=1:4

        % get nutrient data
        nutrient{ip}=squeeze(sum(nutr(nut_inds{ip},:,:,:),1));

        % apply masks
        nutrient{ip}(HFacC==0)=NaN;
        nutrient{ip}(nutrient{ip}<0)=0;

        % multiply concentrations by volume
        DOM{ip} = DOM{ip}.*volume;
        DOM{ip}(HFacC==0)=NaN;
        DOM{ip}(DOM{ip}<0)=0;

        POM{ip} = POM{ip}.*volume;
        POM{ip}(HFacC==0)=NaN;
        POM{ip}(POM{ip}<0)=0;
        
    end

    %% calculate bolus velocity from GM terms

    % Following assumes gmkwx and gmkwy are located at grid centres

    psibx=gmkwx./2; % m^2/s % (W.point, X.dir) of GM-Redi tensor % WHY DIVIDE BY TWO?
    psiby=gmkwy./2; % m^2/s % (W.point, Y.dir) of GM-Redi tensor % WHY DIVIDE BY TWO?

    % differentiate wrt z
    for k=1:nz-1
        vbc(:,:,k) = -(psiby(:,:,k)-psiby(:,:,k+1)) ./ drF(k);
        ubc(:,:,k) = -(psibx(:,:,k)-psibx(:,:,k+1)) ./ drF(k);
    end
    vbc(:,:,nz)=(psiby(:,:,nz)-0)/drF(nz);
    ubc(:,:,nz)=(psibx(:,:,nz)-0)/drF(nz);

    % interpolate to grid edges (same as uvel and vvel)
    vbg(:,1:ny-1,:) = (vbc(:,1:ny-1,:)+vbc(:,2:ny,:))/2;
    vbg(:,  ny  ,:) = (vbc(:,  ny  ,:)+vbc(:,1   ,:))/2;
    ubg(1:nx-1,:,:) = (ubc(1:nx-1,:,:)+ubc(2:nx,:,:))/2;
    ubg(  nx  ,:,:) = (ubc(  nx  ,:,:)+ubc(1   ,:,:))/2;


    %% KPP

    for k=1:nz
        % m^3 s^-1        = m^2 s^-1    * m^2 / m
        kpp_volume(:,:,k) = kpp(:,:,k) .* rA ./ drF(k);
    end

    %% calculate volume transport from product of velocities and face areas
    u_transport = (uvel + ubg) .* u_area; % m s^-1 * m^2 = m^3 s^-1
    v_transport = (vvel + vbg) .* v_area;
    w_transport = (wvel      ) .* w_area;

    % get surface data
    u_surf = u_transport(:,:,1);
    v_surf = v_transport(:,:,1);

    % interpolate to find velocities at centre of grid boxes
    % (only for sanity check)
    u_centre = interp2(y,x,u_surf,cy,cx');
    v_centre = interp2(y,x,v_surf,cy,cx');

    % set land to zero
    u_centre(HFacC(:,:,1)==0)=NaN;
    v_centre(HFacC(:,:,1)==0)=NaN;

    %% Find incoming fluxes to each box from upstream concentrations and fluxes

    for ip=1:4
        % use ocean mask as tracer for now
        tracer = nutrient{ip};
        tracer(isnan(tracer))=0;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % eastwards (i.e. from west to east)
        eastward{ip} = zeros(nx,ny,nz);
        % positive fluxes at western boundary
        upstream = + max(u_transport,0);
        % flux at western edge of domain based on tracer at eastern edge
        eastward{ip}(1,:,:) = upstream(1,:,:) .* tracer(nx,:,:); % m^3 s^-1 * mmol m^-3 = mmol s^-1
        % eastward supply into layer based on flux at western boundary and tracer one to the west
        eastward{ip}(2:nx,:,:) = upstream(2:nx,:,:) .* tracer(1:nx-1,:,:);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % westwards (i.e. from east to west)
        westward{ip} = zeros(nx,ny,nz);
        % negative fluxes at western boundary
        upstream = - min(u_transport,0);
        % westward supply into layer based on flux and tracer one to the east
        westward{ip}(1:nx-1,:,:) = upstream(2:nx,:,:) .* tracer(2:nx,:,:);
        % flux at eastern edge of domain based on flux and tracer at western edge
        westward{ip}(nx,:,:) = upstream(1,:,:) .* tracer(1,:,:);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % northwards (i.e. from south to north)
        northward{ip} = zeros(nx,ny,nz);
        % positive fluxes at southern boundary
        upstream = + max(v_transport,0);
        % flux at southern edge always zero
        northward{ip}(:,1,:) = 0;
        % northward supply into layer based on flux at southern boundary and tracer one to the south
        northward{ip}(:,2:ny,:) = upstream(:,2:ny,:) .* tracer(:,1:ny-1,:);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % southwards (i.e. from north to south)
        southward{ip} = zeros(nx,ny,nz);
        % negative fluxes at southern boundary
        upstream = - min(v_transport,0);
        % southward supply into layer based on flux and tracer one to the north
        southward{ip}(:,1:ny-1,:) = upstream(:,2:ny,:) .* tracer(:,2:ny,:);
        % flux at northern edge always zero
        southward{ip}(:,ny,:) = 0;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % upwards (i.e. from below to above)
        upward{ip} = zeros(nx,ny,nz);
        % positive fluxes at lower boundary
        upstream = + max(w_transport,0);
        % upward supply into layer based on flux at lower boundary and tracer one below
        upward{ip}(:,:,1:nz-1) = upstream(:,:,1:nz-1) .* tracer(:,:,2:nz);
        % upward flux into bottom layer always zero
        upward{ip}(:,:,nz) = 0;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % downwards (i.e. from above to below)
        downward{ip} = zeros(nx,ny,nz);
        % negative fluxes at lower boundary
        upstream = - min(w_transport,0);
        % downward flux into upper layer always zero
        downward{ip}(:,:,1) = 0;
        % downward supply into layer based on flux and tracer one above
        downward{ip}(:,:,2:nz) = upstream(:,:,1:nz-1) .* tracer(:,:,1:nz-1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%

        westward{ip}  = sum(westward{ip} (:,:,1:n_layers),3);
        eastward{ip}  = sum(eastward{ip} (:,:,1:n_layers),3);
        northward{ip} = sum(northward{ip}(:,:,1:n_layers),3);
        southward{ip} = sum(southward{ip}(:,:,1:n_layers),3);
        upward{ip}    = upward{ip}(:,:,n_layers);

        horizontal_100{ip} = westward{ip}   ...
                           + eastward{ip}   ...
                           + northward{ip}  ...
                           + southward{ip};

        vertical_100{ip} = upward{ip};

        % calculate remineralisation (currently ingnoring temperature)
        remineralisation = (DOM{ip}.*0.0333 + POM{ip} .* 0.033); % mmol day^-1
        % convert from mmoles/d to mmol m^-2 s^-1
        remineralisation = remineralisation ./ 86400 ./ rA;
        remin_flux{i_month,ip} =  sum(remineralisation(:,:,1:n_layers),3); % mmol m^-2 s^-1

        % add together all positive incoming fluxes
        % mmol s^-1
        advective_flux{i_month,ip} = horizontal_100{ip} + vertical_100{ip}; % mmol s^-1

        % get flux per unit area (mmol s^-1 / m^2 = mmol m^-2 s^-1)
        advective_flux{i_month,ip} = advective_flux{i_month,ip} ./rA;

        %  m^2 s^-1 * mmol m^-3 / m = mmol m^-2
        diffusive_flux{i_month,ip} = kpp(:,:,n_layers+1) .* (tracer(:,:,n_layers+1) - 0.*tracer(:,:,n_layers)) ./ drC(n_layers);

        % Mask land with zeros
        advective_flux{i_month,ip}(HFacC(:,:,1)==0)=0;
        diffusive_flux{i_month,ip}(HFacC(:,:,1)==0)=0;
        % Mask land with Infs for plotting
        mask_advective{i_month,ip} = advective_flux{i_month,ip};
        mask_diffusive{i_month,ip} = diffusive_flux{i_month,ip};
        mask_advective{i_month,ip}(HFacC(:,:,1)==0)=Inf;
        mask_diffusive{i_month,ip}(HFacC(:,:,1)==0)=Inf;

        % Atmospheric Fe
        if ip==3
            dust_flux{i_month,ip} = dust(:,:,i_month).*1000.*0.04;
        else
            dust_flux{i_month,ip} = zeros(360,160);
        end

        total_flux{i_month,ip} = advective_flux{i_month,ip} ...
                               + diffusive_flux{i_month,ip} ...
                               +     remin_flux{i_month,ip} ...
                               +      dust_flux{i_month,ip};
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here endeth the calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% save fluxes

save('../Data/explicit_fluxes.mat','advective_flux','diffusive_flux','remin_flux','dust_flux','total_flux')


