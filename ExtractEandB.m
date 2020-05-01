function ExtractEandB(coilNum, voxelizedMesh)
% EXTRACT E AND B FIELD into variables
      
    % load voxelized mesh and sensorDims
    load('voxelizedMesh','x_mesh','y_mesh','z_mesh','frankMask')
    % brainMask, frankMask, sigma, epsilon,
    Bpath = fullfile(pwd,'Bfield');
    Epath = fullfile(pwd,'Efield');
    etaEpath = fullfile(Epath,'etaEvals');
    sarEpath = fullfile(Epath,'sarEvals');
    

    for m=1:coilNum

        M = m-1;
        filename = [num2str(M) 'Rho.mat'];
        fullpath = fullfile(Bpath, filename);
        % check to see if it exists...
        if ~exist(fullpath,'file')
           fprintf('%s doesnt exist...\n', filename)
        end
        load(fullpath)

        % find the shared locations and their index in the mesh
        tol = 10e-4/max(abs([xdim(:);x_mesh(:)]));
        [A,B] = ismembertol(x_mesh,xdim,tol,'OutputAllIndices',true);
        indexMesh_x = find(A==1);
        indexDim_x = nonzeros(cell2mat(B));
        tol = 10e-4/max(abs([ydim(:);y_mesh(:)]));
        [A,B] = ismembertol(y_mesh,ydim,tol,'OutputAllIndices',true);
        indexMesh_y = find(A==1);
        indexDim_y = nonzeros(cell2mat(B));
        tol = 10e-4/max(abs([zdim(:);z_mesh(:)]));
        [A,B] = ismembertol(z_mesh,round(zdim),tol,'OutputAllIndices',true);
        indexMesh_z = find(A==1);
        indexDim_z = nonzeros(cell2mat(B));
        % add these to the matrices
%         meshx_m(:,m) = indexMesh_x;
%         meshy_m(:,m) = indexMesh_y;
%         meshz_m(:,m) = indexMesh_z;
%         dimx_m(:,m) = indexDim_x;
%         dimy_m(:,m) = indexDim_y;
%         dimz_m(:,m) = indexDim_z;


        % Bfield nonzero index
        nonZB_v = find(frankMask(indexMesh_z,indexMesh_y,indexMesh_x));
        %a = brainMask(indexMesh_z,indexMesh_y,indexMesh_x);
        % Save the Bfield where the brainMask is!
        B1plus_m(:,m) = 10e6*(BxField(nonZB_v)+1i*ByField(nonZB_v))/2;
        fprintf('Coil %d loaded and placed in B1plus_m...\n', M)



    end
    save('B1plus_m.mat','B1plus_m','-v7.3')
    
    % preallocate for E since it might be a little bigger
        M=1;
        filename = [num2str(M) 'Rho.mat'];
        fullpath = fullfile(Epath, filename);
        % check to see if it exists...
        if ~exist(fullpath,'file')
           fprintf('%s doesnt exist...\n', filename)
        end
        load(fullpath,'xdim','ydim','zdim')

        % maintain Dimensionality...
        % find the shared locations and their index in the mesh
        tol = 10e-4/max(abs([xdim(:);x_mesh(:)]));
        A = ismembertol(x_mesh,xdim,tol,'OutputAllIndices',true);
        indexMesh_x = find(A==1);
        tol = 10e-4/max(abs([ydim(:);y_mesh(:)]));
        A = ismembertol(y_mesh,ydim,tol,'OutputAllIndices',true);
        indexMesh_y = find(A==1);
        tol = 10e-2/max(abs([zdim(:);z_mesh(:)]));
        A = ismembertol(z_mesh,round(zdim),tol,'OutputAllIndices',true);
        indexMesh_z = find(A==1);
    
    sensorFM = frankMask(indexMesh_z(1:numel(indexMesh_z)-1),indexMesh_y(1:numel(indexMesh_y)-1),indexMesh_x(1:numel(indexMesh_x)-1));
    nonZE_v = find(sensorFM);
    sarE_rho = zeros(length(nonZE_v),coilNum);
    etaE_rho = zeros(length(nonZE_v),coilNum);

    % E field stuff
    for m=1:coilNum

        M = m-1;
        filename = [num2str(M) 'Rho.mat'];
        fullpath = fullfile(Epath, filename);
        % check to see if it exists...
        if ~exist(fullpath,'file')
           fprintf('%s doesnt exist...\n', filename)
        end
        load(fullpath,'ExField','EyField','EzField','xdim','ydim','zdim')

        % maintain Dimensionality...
        % find the shared locations and their index in the mesh
        tol = 10e-4/max(abs([xdim(:);x_mesh(:)]));
        [A,B] = ismembertol(x_mesh,xdim,tol,'OutputAllIndices',true);
        indexMesh_x = find(A==1);
        indexDim_x = nonzeros(cell2mat(B));
        tol = 10e-4/max(abs([ydim(:);y_mesh(:)]));
        [A,B] = ismembertol(y_mesh,ydim,tol,'OutputAllIndices',true);
        indexMesh_y = find(A==1);
        indexDim_y = nonzeros(cell2mat(B));
        tol = 10e-2/max(abs([zdim(:);z_mesh(:)]));
        [A,B] = ismembertol(z_mesh,round(zdim),tol,'OutputAllIndices',true);
        indexMesh_z = find(A==1);
        indexDim_z = nonzeros(cell2mat(B));
        % add these to the matrices
%         meshx_m(:,m) = indexMesh_x;
%         meshy_m(:,m) = indexMesh_y;
%         meshz_m(:,m) = indexMesh_z;
%         dimx_m(:,m) = indexDim_x;
%         dimy_m(:,m) = indexDim_y;
%         dimz_m(:,m) = indexDim_z;


        % Revoxelization
        [zd,yd,xd] = size(ExField);
        sarE = zeros(zd-1,yd-1,xd-1);
        etaE = zeros(zd-1,yd-1,xd-1);
        
%         fprintf('Caclulating Efield/etaE/sarE for %s\n',filename)
%         for i=1:zd-1
%            for j = 1:yd-1
%                for k = 1:xd-1
% 
%                % for every combination
%                sarE(i,j,k) = ...
%                    abs((ExField(i,j,k) + ExField(i,j+1,k)+ ExField(i,j,k+1)+ ExField(i,j+1,k+1))/4)^2 ...
%                     + abs((EyField(i,j,k) + EyField(i+1,j,k) + EyField(i,j,k+1)+ EyField(i+1,j,k+1))/4)^2 ...
%                     + abs((EzField(i,j,k) + EzField(i+1,j,k) + EzField(i,j+1,k)+ EzField(i+1,j+1,k))/4)^2;
%                % for every combination
%                 etaE(i,j,k) = ...
%                    sqrt(((ExField(i,j,k) + ExField(i,j+1,k)+ ExField(i,j,k+1)+ ExField(i,j+1,k+1))/4)^2 ...
%                     +((EyField(i,j,k) + EyField(i+1,j,k) + EyField(i,j,k+1)+ EyField(i+1,j,k+1))/4)^2 ...
%                     +((EzField(i,j,k) + EzField(i+1,j,k) + EzField(i,j+1,k)+ EzField(i+1,j+1,k))/4)^2);
% 
% 
% 
%                end
%            end
%         end

%         % designate the nonZE_v
%         % sarE_rho(:,m) = sarE(nonZE_v);
%         % etaE_rho(:,m) = etaE(nonZE_v);
%         name = sprintf('%dsarE.mat',M);
%         full = fullfile(sarEpath,name);
%         save(full,'sarE')
%         name = sprintf('%detaE.mat',M);
%         full = fullfile(etaEpath,name);
%         save(full,'etaE')

    end


    %
    save('Indices.mat','indexMesh_x','indexMesh_y','indexMesh_z')
    %save('Edata.mat','sarE_rho','etaE_rho','-v7.3')

end





