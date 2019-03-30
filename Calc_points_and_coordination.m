
count_difference_this_member = population(initialisations_converged + 1,(sequence_length+1):(2*sequence_length))' ; % used internally.
        
% Recalculate best z and output results:
x_cords = zeros(1,original_sum + sum(count_difference_this_member)) ; % These are initialised inside the loop to allow for changes in the number of atoms.
y_cords = zeros(1,original_sum + sum(count_difference_this_member)) ;
z_cords =  zeros(1,original_sum + sum(count_difference_this_member)) ;
atom_count = 0 ;
for column = 1:column_count
    for i = 1:(count_difference_this_member(column) + count_bank(column))
        atom_count = atom_count + 1 ; % Increase the number of atoms by one.
        x_cords(atom_count) = point_coordinates(column,1) * pixelWidth ;
        y_cords(atom_count) = point_coordinates(column,2) * pixelWidth ;
        %         z_cords(atom_count) = (i + population(member,column) + z_offsets_bank(column)) * Rz ;
        z_cords(atom_count) = (i + population(initialisations_converged + 1,column) ) * Rz ;
        
    end
end

% Create simple CN data:
observations = [x_cords ; y_cords ; z_cords ]' ; % Create xyz by n matrix:
atom_neighbour_count = zeros(atom_count,1) ;
if atom_count < max_atoms_to_draw
    all_distances = squareform(pdist(observations)) ;
    all_distances(all_distances==0) = 10 ; % Force self distance to be large.
    atom_count = size(observations,1) ;
    for i = 1:atom_count
        this_atom_distances = all_distances(:,i) ;
        this_atom_neighbours = (this_atom_distances<=1.207*R0) ;
        atom_neighbour_count(i,1) = sum(this_atom_neighbours) ;
    end
else
    progress = 0 ;
    h = waitbar(progress,'Calculating coordination numbers...') ;
    for test_atom = 1:atom_count
        test_observations = bsxfun(@minus, observations, observations(test_atom,:)) ;
        test_distances = sum( test_observations.^2 , 2) ; % This is kept as the squared euclidean distance!!
        atom_neighbour_count(test_atom,1) = sum((test_distances < ((1.3*R0).^2) )) -1 ;
        if mod(test_atom,100) == 0
            progress = test_atom / atom_count ;
            waitbar(progress,h)
        end
    end
    delete(h)
end


figure('Position',[100 100 500 500]);
disp('Drawing best current solution...')
clf
run('C:\Users\Lewys\Box Sync\DPhil\MatLab Software\Absolute Integrator\3D Reconstruction\SingleDraw.m')
title(num2str(generation))
zlim([-1.5 1.5])
G(generation) = getframe(gcf);
close(gcf)