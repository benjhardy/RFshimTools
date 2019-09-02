function eta_v = calculateEta(wfull_m, ab1ps, etaE, voxelizedMesh, Indices)
    % 
    % Calculates ETA
    % Benjamin M Hardy 8-29-2019
    % inputs:
    %        wfull_m: shim values
    %        B1plus_m: B1plus_m values
    %        etaE electrical field for ETA
    %        voxelizedMesh: struct containing all the masks and conduct.
    %        Indices: the indices of the mesh that correspond to the sensor
    % Outputs: the string of eta values in uT^2/Watt
    
    % load structs
    sigma_m = voxelizedMesh.sigma_m;
    frankMask = voxelizedMesh.frankMask;
    indexMesh_x = Indices.indexMesh_x;
    indexMesh_y = Indices.indexMesh_y;
    indexMesh_z = Indices.indexMesh_z;
    
    % number of coils and number of vars
    [Nc,numvars] = size(wfull_m);
    
    % this just ensures the sigma is the same size as the Edata
    sigma_m = sigma_m(indexMesh_z(1:end-1),indexMesh_y(1:end-1),indexMesh_x(1:end-1));
    sensorFrankMask = frankMask(indexMesh_z(1:end-1),indexMesh_y(1:end-1),indexMesh_x(1:end-1));
    sigma_m = sigma_m(sensorFrankMask);
    
    %repSigma_v = repmat(sigma_m,[1,Nc]);
    
    % Pre-allocate
    eta_v = zeros(numvars,1);
    
   % set the number of times to loop through and calculate
    for i = 1:numvars 
        
        w = wfull_m(:,i);
        fprintf('Calculating Eta for the %d element array Shim\n',Nc)
        phi_m = ((.001^3)*(.5*sigma_m).*etaE)'*etaE;
        P = w'*phi_m*w;
        eta = ab1ps(i)/P;
        clear phi_m
        eta_v(i) = eta;

    end


end