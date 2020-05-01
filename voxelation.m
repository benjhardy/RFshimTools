function voxelation(fullheadMesh)
% Input: fullheadMesh matlab object that contains various mesh information
% Output: save a file to the workspace that contains the following
%       rho_m - the density matrix giving the density of the materials at
%       the mesh locations found in x_mesh, y_mesh, z_mesh
%       sigma_m epsilon_m - conductivity and permittivity in mesh
%       frankMask, brainMask, - total face mask, and total brain mask in
%       mesh
%       x_mesh, y_mesh, z_mesh hold the mesh locations
       
    MeshExDensity = fullheadMesh.MeshExDensity;                             
    MeshExEpsilon_r = fullheadMesh.MeshExEpsilon_r;                           
    MeshExSigma = fullheadMesh.MeshExSigma;                           
    MeshEyDensity = fullheadMesh.MeshEyDensity;                          
    MeshEyEpsilon_r = fullheadMesh.MeshEyEpsilon_r;                             
    MeshEySigma = fullheadMesh.MeshEySigma;                     
    MeshEzDensity = fullheadMesh.MeshEzDensity;                       
    MeshEzEpsilon_r = fullheadMesh.MeshEzEpsilon_r;                        
    MeshEzSigma = fullheadMesh.MeshEzSigma;                   
    grid_X = fullheadMesh.grid_X;                
    grid_Y = fullheadMesh.grid_Y;                       
    grid_Z = fullheadMesh.grid_Z;
   
    %load('C:\Users\benja\OneDrive - Vanderbilt\Documents\MATLAB\RFshimTools\densityValsBrain.mat', 'densityValsBrain')
    load('F:\5634\yuruiDensityVals.mat','yuruiDensityVals')
    
    densityValsBrain = yuruiDensityVals;
    
    % Density Mesh creation...
    xD = MeshExDensity;
    yD = MeshEyDensity;
    zD = MeshEzDensity;
    % xD(isnan(xD))=0;
    % yD(isnan(yD))=0;
    % zD(isnan(zD))=0;

    rho_m = zeros(length(grid_Z)-1,length(grid_Y)-1,length(grid_X)-1);

    for i=1:length(grid_Z)-1
       for j = 1:length(grid_Y)-1
           for k = 1:length(grid_X)-1

           % for every combination
           rho_m(i,j,k) = ...
               (xD(i,j,k) + xD(i,j+1,k)+ xD(i,j,k+1)+ xD(i,j+1,k+1)...
                + yD(i,j,k) + yD(i+1,j,k) + yD(i,j,k+1)+ yD(i+1,j,k+1)...
                + zD(i,j,k) + zD(i+1,j,k) + zD(i,j+1,k)+ zD(i+1,j+1,k))/12;



           end
       end
    end
    
    %rho_m(isnan(rho_m)) = 0;
    % select the density values of interest to create the brainMask
    % figure
    % imagesc(squeeze(rho_m(120,:,:)))
    % % find the vals in densityValsBrain
    % a = roipoly;
    % colorbar
    % %
    % vals = unique(a.*squeeze(rho_m(150,:,:)));
    % vals(isnan(vals)) = 0;
    % vals = unique(vals);

    % ****** LOAD this first
    vals = densityValsBrain;
    i = ismember(rho_m,vals);
    brainMask = rho_m.*i;
    brainMask(isnan(brainMask)) = 0;
    %
    brainMask = bwareaopen(brainMask,10000);
    brainMask = imfill(brainMask,'holes');

    %cancel the bottom z slices
    %temp = 0*brainMask;
    %temp(176:end,:,:) = 1;
    %brainMask = brainMask.*temp;

    zMovie(brainMask,250)

    %
    % frankMask is anything but freespace minus the shield?
    frankMask = ~isnan(rho_m);
    frankMask = imfill(frankMask,'holes');

    zMovie(frankMask,250);


    % Sigma and Epsilon...
    % epsilon
    xE = MeshExEpsilon_r;
    yE = MeshEyEpsilon_r;
    zE = MeshEzEpsilon_r;
    % sigma
    xS = MeshExSigma;
    yS = MeshEySigma;
    zS = MeshEzSigma;

    sigma_m = zeros(length(grid_Z)-1,length(grid_Y)-1,length(grid_X)-1);
    epsilonr_m = zeros(length(grid_Z)-1,length(grid_Y)-1,length(grid_X)-1);

    for i=1:length(grid_Z)-1
       for j = 1:length(grid_Y)-1
           for k = 1:length(grid_X)-1

           % for every combination
           % sigma
           sigma_m(i,j,k) = ...
               (xS(i,j,k) + xS(i,j+1,k)+ xS(i,j,k+1)+ xS(i,j+1,k+1)...
                + yS(i,j,k) + yS(i+1,j,k) + yS(i,j,k+1)+ yS(i+1,j,k+1)...
                + zS(i,j,k) + zS(i+1,j,k) + zS(i,j+1,k)+ zS(i+1,j+1,k))/12;
           % for every combination
           % epsilonr
           epsilonr_m(i,j,k) = ...
               (xE(i,j,k) + xE(i,j+1,k)+ xE(i,j,k+1)+ xE(i,j+1,k+1)...
                + yE(i,j,k) + yE(i+1,j,k) + yE(i,j,k+1)+ yE(i+1,j,k+1)...
                + zE(i,j,k) + zE(i+1,j,k) + zE(i,j+1,k)+ zE(i+1,j+1,k))/12;




           end
       end
    end
    %
    epsilonr_m = frankMask.*epsilonr_m;
    zMovie(epsilonr_m,250)
    
    %
    sigma_m = sigma_m.*frankMask;
    zMovie(sigma_m,250)
    
    % Recreate the Grid...
    % 
    % gz = repmat(grid_Z,[1,length(grid_Y),length(grid_X)]);
    % gy = repmat(grid_Y,[1,length(grid_X),length(grid_Z)]);
    % gy = permute(gy,[3,1,2]);
    % gx = repmat(grid_X,[1,length(grid_Y),length(grid_Z)]);
    % gx = permute(gx,[3,2,1]);
    % 
    % for i=1:length(grid_Z)-1
    %    for j = 1:length(grid_Y)-1
    %        for k = 1:length(grid_X)-1
    %            
    %        % for every combination
    %        % same format as E field, new points for the revoxelized data.
    %        gridx(i,j,k) = ...
    %            (gx(i,j,k) + gx(i,j+1,k)+ gx(i,j,k+1)+ gx(i,j+1,k+1))/4;
    %        gridy(i,j,k) = ...
    %            (gy(i,j,k) + gy(i+1,j,k) + gy(i,j,k+1)+ gy(i+1,j,k+1))/4;
    %        gridz(i,j,k) = ...
    %            (gz(i,j,k) + gz(i+1,j,k) + gz(i,j+1,k)+ gz(i+1,j+1,k))/4;
    %             
    %        
    %         
    %        
    %        
    %        
    %        end
    %    end
    % end
    % x_mesh = unique(gridy);
    % y_mesh = unique(gridz);
    % z_mesh = unique(gridx);

    % % % recreate the Sensor Grid
    % % load SolidSensorDef.mat
    % gz = repmat(Z_Dimension_1,[1,length(Y_Dimension_2),length(X_Dimension_3)]);
    % gy = repmat(Y_Dimension_2,[1,length(X_Dimension_3),length(Z_Dimension_1)]);
    % gy = permute(gy,[3,1,2]);
    % gx = repmat(X_Dimension_3,[1,length(Y_Dimension_2),length(Z_Dimension_1)]);
    % gx = permute(gx,[3,2,1]);
    % 
    % for i=1:length(Z_Dimension_1)-1
    %    for j = 1:length(Y_Dimension_2)-1
    %        for k = 1:length(X_Dimension_3)-1
    %            
    %        % for every combination
    %        % same format as E field, new points for the revoxelized data.
    %        gridx(i,j,k) = ...
    %            (gx(i,j,k) + gx(i,j+1,k)+ gx(i,j,k+1)+ gx(i,j+1,k+1))/4;
    %        gridy(i,j,k) = ...
    %            (gy(i,j,k) + gy(i+1,j,k) + gy(i,j,k+1)+ gy(i+1,j,k+1))/4;
    %        gridz(i,j,k) = ...
    %            (gz(i,j,k) + gz(i+1,j,k) + gz(i,j+1,k)+ gz(i+1,j+1,k))/4;
    %             
    %        
    %         
    %        
    %        
    %        
    %        end
    %    end
    % end
    % y = unique(gridy);
    % z = unique(gridz);
    % x = unique(gridx);
    x_mesh = grid_X;
    y_mesh = grid_Y;
    z_mesh = grid_Z;
    close all
    clc
save('voxelizedMesh.mat','rho_m','sigma_m','epsilonr_m','frankMask','brainMask','x_mesh','y_mesh','z_mesh')