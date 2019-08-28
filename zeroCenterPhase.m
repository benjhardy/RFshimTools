function phs_cp = zeroCenterPhase(B1plus_m)
 
    % The phi[s]_cp is the spatial phase pattern of the sum of the coils’ B1+ maps 
    % after aligning the phases of all the coils at some voxel or ROI of voxels 
    % in the middle of the volume of interest. - Dr. Will Grissom
    % Input: B1plus_m
    % Output: phs array at every B1plus_m point

    % Return: a spatial phase map the same size as the B1plus_m space

    % find center voxel, coil number
    centerVoxel = round(size(B1plus_m,1)/2);

    % for every coil find the phase that results in zero phase at voxel
    % center

    phs = angle(B1plus_m(centerVoxel,:));
    
    phaseWeight = exp(-phs*1i).';
    
    % check the summed phase at the center of the voxel is indeed zero
    % diff = angle(B1plus_m(centerVoxel,:)*phaseWeight);
    
    % now apply the phase to the mapping and return the phase map...
    phs_cp = angle(B1plus_m*phaseWeight);





end