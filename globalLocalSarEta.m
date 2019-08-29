function [gSar_v,lSar_m,eta_v] = globalLocalSarEta(wfull_m,B1plus_m, etaE, sarE, voxelizedMesh, Indices, flag)
    % calculates local and global sar given the shim data
    % also calculates eta
    % Input: wfull_m - array of shims(Number of coils, some var)
    %        B1plus_m - specific array of B1plus_m values.
    %        Edata - etaE and sarE to calculate eta and SAR
    %        voxelizedMesh - struct stuff needed to calculate eta and SAR
    %        Indices - struct containing indexMesh_x,y,z
    %        flag - true if you want 10 gram SAR but it might take a while
    %        Output: local, global and Eta for the array
    
    % load the structs
    load('voxelizedMesh','sigma_m', 'rho_m','frankMask');
    load('Indices', 'indexMesh_x','indexMesh_y','indexMesh_z');
    
    [Nc,numvars] = size(wfull_m);
    
    % take only rho and sigma at points of sensor.
    % when dealing with E field values the field had to be linearly
    % interpolated, so the mesh points need to be adjusted a little.
    % they have already been multiplied by the mask in voxelation function
    sigma_m = sigma_m(indexMesh_z(1:end-1),indexMesh_y(1:end-1),indexMesh_x(1:end-1));
    rho_m = rho_m(indexMesh_z(1:end-1),indexMesh_y(1:end-1),indexMesh_x(1:end-1));
    
    sensorFrankMask = frankMask(indexMesh_z(1:end-1),indexMesh_y(1:end-1),indexMesh_x(1:end-1));
    sigma_v = sigma_m(sensorFrankMask);
    rho_v = rho_m(sensorFrankMask);
    
    repSigma_v = repmat(sigma_v,[1,Nc]);
    repRho_v = repmat(rho_v,[1,Nc]);
    
    % Pre-allocate
    eta_v = zeros(numvars,1);
    gSar_v = zeros(numvars,1);
    lSar_m = zeros(size(sensorFrankMask,1),size(sensorFrankMask,2),size(sensorFrankMask,3),numvars);
    
   % set the number of times to loop through and calculate
    for i = 1:numvars 
        
        w = wfull_m(:,i);
        % First Calculate ETA
        % aveB1plusSquared
        Nc = size(B1plus_m,2);
        Np = size(B1plus_m,1);
        G = (1/Np)*(B1plus_m'*B1plus_m);
        ab1ps = w'*G*w;
        clear G
        
        fprintf('Calculating Eta for the %d element array Shim\n',Nc)
        C = (.5*repSigma_v).*etaE;
        phi_m = C'*etaE*(.001^3);
        P = w'*phi_m*w;
        eta = ab1ps/P;
        eta_v(i) = eta;

        % calculate global SAR
        fprintf('Calculating global SAR data for the %d element array Shim\n',Nc)
        % apply the shim to all the E
        sarEshim = sarE*((abs(w).^2));
        C = .5*(1/size(sarE,1)*.001^3)*(.001^3);
        SAR = C*((1./rho_v).*sigma_v.*sarEshim);
        % finally, this means the density is 0 or undefined or air in some
        % parts of my mask, this makes sense, so at these points, the SAR value
        % is entirely ignored...
        SAR(isnan(SAR)) = 0;
        gSar = sum(SAR);
        gSar_v(i) = gSar;
        
        if flag == true
            % calculate local 10 gram SAR . . .
            fprintf('Calculating local SAR data for the %d element array Shim\n',Nc)
            wE_m = 0*sensorFrankMask;
            C = sigma_v./rho_v;
            wE_m(sensorFrankMask) = C.*sarEshim;
            wE_m(isnan(wE_m)) = 0;
            % calls NYU n gram SAR average tool.
            lSar = SARavesphnp(rho_m, wE_m, .001,.001,.001,10);
            lSar_m(:,:,:,i) = lSar;
        end
        fprintf('Calculation set %d/%d done...\n',i,numvars)
    end



end