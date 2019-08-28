function [wfull_m, RSD_v, maxI_v] = tikhonovShim(B1plus_m, betas, stopWhen, FOX, voxelizedMesh, Indices, phase, d)
% Tikhonov Regularization Function
% source Code: https://vuiis.vumc.org/~grissowa/teaching.html
% *click pulse design code
% OUTPUT : [wfull_m, RSD_v, maxI_v]
% Input: B1plus_m - the array of Coils B1plus complex field values
%        Indices - the points in the field that the field has been
%        measured, or the points in the mesh that the sensor box measured
%        the B and E field
%        beta -  a string of beta values to shim the array over
%        stopWhen - a value that specifies when to stop the iterative phase
%        optimization, when the abs. difference between the previous cost
%        function and present cost function hasn't changed
%        FOX - the field of Excitation, matrix where the Bfield corresponds
%        in  a linearized fashion
%        voxelizedMesh - Struct containing the mesh points (may need to
%        change later)
%        Indices - the points within the mesh where the sensor is located
%        phase - a string specifying what the initial phase should be...
%        choices 'random', 'zeroSum', 'phs_cp'
%        d - target field strength in uT

    % General Parameters
    % to minimize ||Ax-d||^2 + beta/2 * ||x||^2
    
    % number of points to homogenize over, try frankMask at first
    Np = size(B1plus_m,1);
    Nc = size(B1plus_m,2); % number of coils

    % size of matrix
    [zdim,ydim,xdim] = size(FOX); % in mm

    % step 1: setup design grid and kx,ky,kz search grid
    load('voxelizedMesh','z_mesh','y_mesh','x_mesh'); 
    load('Indices');
    [z,y,x] = ndgrid(z_mesh(indexMesh_z),y_mesh(indexMesh_y),x_mesh(indexMesh_x));
%     deltamax = 1; % mm, max res of traj
%     kmax = 1/deltamax; %/mm max spatial freq of traj
%     [kzs,kys,kxs] = ndgrid(-kmax/2:1/zdim:kmax/2-1/zdim, -kmax/2:1/ydim:kmax/2-1/ydim, -kmax/2:1/xdim:kmax/2-1/xdim);
%     kxs = kxs(:); kys = kys(:); kzs = kzs(:);
%     dc = find((kxs == 0) & (kys == 0) & (kzs == 0));
%     kxs = kxs([1:dc-1 dc+1:end]); kys = kys([1:dc-1 dc+1:end]); kzs = kzs([1:dc-1 dc+1:end]); % remove dc points

    % step 2 design weights
    kx = [0];
    ky = [0];
    kz = [0];
    kxlin=kx;kylin=ky;kzlin=kz;

    % STEP 3 build system matrix
    A = exp(1i*2*pi*(x(FOX)*kx(:)' + y(FOX)*ky(:)' +z(FOX)*kz(:)'));
    Afull = zeros(Np,Nc);
    for j = 1:Nc
        % multiply b1maps and fourier kernel
        Afull(:,j) = bsxfun(@times,B1plus_m(:,j),A);
        fprintf('Fourier Kernel multiplied for coil %d\n',j)
    end
    
    % Step 4:
    % Set initial Phase
    if strcmp(phase,'zeroSum')
        phs = zeros(Np,1);
    end
    if strcmp(phase,'phs_cp')
        phs = zeroCenterPhase(B1plus_m);
    end
    if strcmp(phase, 'random')
        phs = 2*pi*rand(Np,1) - pi; 
    end
    
    % preallocation for speed:
    num = length(betas);
    RSD_v = zeros(num,1);
    maxI_v = zeros(num,1);
    wfull_m = zeros(Nc,num);
    
    % Loop through each beta value and perform the shim...
    for i = 1:num
        
        beta = betas(i);
        % STEP 5:
        % initial parameters for iterative phase solution
        wfull = inv(Afull'*Afull + beta*eye(Nc))*(Afull'*(d*exp(1i*phs)));
        err = Afull*wfull - d*exp(1i*phs);
        cost = real(err'*err +beta*wfull'*wfull);
        costold = 10*cost;

        % loop through until cost stops changing
        while abs(costold-cost) > stopWhen*costold
               costold = cost;
               phs = angle(Afull*wfull);
               wfull = inv(Afull'*Afull + beta*eye(Nc))*(Afull'*(d*exp(1i*phs)));
               err = Afull*wfull - d*exp(1i*phs);
               cost = real(err'*err + beta*wfull'*wfull);
        end

        RSD_v(i) = std(abs(Afull*wfull))/mean(abs(Afull*wfull));
        maxI_v(i) = abs(max(wfull));
        wfull_m(:,i) = wfull;
        fprintf('Beta value %d shimmed, %d/%d shims complete.\n',beta, i, num)
    end

end