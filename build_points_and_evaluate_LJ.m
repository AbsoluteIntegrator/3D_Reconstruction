% Calculate cost function for new population.
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
        z_cords(atom_count) = (i + population(member,column) ) * Rz ;
    end
end

max_neighbours_in_cutoff = 42 ;
X = (1/R0) * [x_cords ; y_cords ; z_cords]' ;
[IDX,fractional_distances] = knnsearch(X,X,'K',max_neighbours_in_cutoff) ;

site_energy = sum( LJ_A*(fractional_distances(:,2:end).^-LJ_a) + LJ_B*(fractional_distances(:,2:end).^-LJ_b) , 2) ;

% 5 columns added for total-energy, total atoms, energy-per-atom, cost-function, and model-probability.
population(member,(2*sequence_length)+1) = 0.5 * sum(site_energy) ; % Divide 2 becuase distances appear twice.
population(member,(2*sequence_length)+2) = atom_count ;
population(member,(2*sequence_length)+3) = population(member,(2*sequence_length)+1)/population(member,(2*sequence_length)+2) ;

% This is the cost_function:

% Loop over the columns and get the global model probability:
model_global_probability = 1 ; % Initialise with proability 1.
for column = 1:column_count
    % Extract probability of each column's atom count and include in running product:
    % FOR "proBMatrix_final" that begins at 0 atom, hence +1 in the index.
    model_global_probability = model_global_probability * probMatrix_final((count_difference_this_member(column) + count_bank(column) + 1) , column) ;
end
% Normalise this to reflect the number of columns:
model_global_probability = model_global_probability .^ (1/column_count) ;
% model_global_probability = -log(model_global_probability) ;

% Create cost-funtion:
population(member,(2*sequence_length)+4) = model_global_probability * population(member,(2*sequence_length)+3) ; % Calculate the neergy per atom without new summations.

% Store probability:
population(member,(2*sequence_length)+5) = model_global_probability ;

