

function [B1plus_m, FOX] = B1plusReform(B1plus_m, voxelizedMesh,Indices)
    % Purpose: Restructures the B1plus_m array into an array where only 
    % the brain is included
    %
    % Input: the original B1plus_m
    %        voxelizedMesh -  struct containing the brain and frankMask
    %        Indices - the indices of mesh where the sensor is located
    %        flag to indicate whether or not the operation should be done
    %        or not, may be redundant, yet feels like a safety feature to
    %        me
    % Output: the B1plus_m only present for the brainMask
    
    % first check to see if the B1plus_m has already been restructured
    load('voxelizedMesh','frankMask')
    load('Indices')
    sensorfrankMask = frankMask(indexMesh_z,indexMesh_y,indexMesh_x);
    ptsfm = length(nonzeros(sensorfrankMask));
    pts = size(B1plus_m,1);
    flag = true;
    
    if (pts ~= ptsfm)
       fprintf('B1plus_m has already been restructed and has length %d\n',pts) 
       flag = false;
    end
    
    load('voxelizedMesh', 'brainMask');
    FOX = logical(brainMask(indexMesh_z, indexMesh_y, indexMesh_x));
    
    if (flag == true)
        
        temp = sensorfrankMask + FOX;
        temp = nonzeros(temp(:));
        brainIn = find(temp == 2);
        B1plus_m = B1plus_m(brainIn,:);
    
    end
    
end