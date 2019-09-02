function [local_SAR, global_SAR] = localGlobalSAR(wfull_m, sarE, voxelizedMesh,Indices, flag )
    % 
    % Calculates local and global SAR
    % Benjamin M Hardy 8-29-2019
    % inputs:
    %        wfull_m: shim values
    %        sarE: electrical field values catered to sar calculation
    %        voxelizedMesh: struct containing all the masks and conduct.
    %        Indices: the indices of the mesh that correspond to the sensor
    %        flag, true for the local SAR as well
    % Outputs: the string of SAR values Watts/kg
    
    
    % load info
    sigma_m = voxelizedMesh.sigma_m;
    rho_m = voxelizedMesh.rho_m;
    frankMask = voxelizedMesh.frankMask;
    indexMesh_x = Indices.indexMesh_x;
    indexMesh_y = Indices.indexMesh_y;
    indexMesh_z = Indices.indexMesh_z;
    
    
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
    
    % Pre-allocate
    global_SAR = zeros(numvars,1);
    local_SAR = zeros(size(sensorFrankMask,1),size(sensorFrankMask,2),size(sensorFrankMask,3),numvars);

% set the number of times to loop through and calculate
    for i = 1:numvars 
        
        w = wfull_m(:,i);
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
        global_SAR(i) = gSar;
        
        
        if flag == true
            % calculate local 10 gram SAR . . .
            fprintf('Calculating local SAR data for the %d element array Shim\n',Nc)
            wE_m = 0*sensorFrankMask;
            C = sigma_v./rho_v;
            wE_m(sensorFrankMask) = C.*sarEshim;
            wE_m(isnan(wE_m)) = 0;
            % calls NYU n gram SAR average tool.
            lSar = SARavesphnp(rho_m, wE_m, .001,.001,.001,10);
            local_SAR(:,:,:,i) = lSar;
        end
        fprintf('Calculation set %d/%d done...\n',i,numvars)
        
        
    end



end






