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

% %***************************
max_neighbours_in_cutoff = 42 ;
X = 10 * [x_cords ; y_cords ; z_cords]' ; % Note conversion here to ANGSTROM.
[IDX,angstrom_distances] = knnsearch(X,X,'K',max_neighbours_in_cutoff) ;
% angstrom_distances = pdist (X) ; % This is NOT squared distances but Euclidean.
% angstrom_distances = squareform(angstrom_distances) ; % Reshape this into a square, now represents the i,j neighbour distances.
% angstrom_distances = sort(angstrom_distances,2,'ascend') ;
% Find index of furthest site worth checking:
% is_close = angstrom_distances<cutoff.*(eye(atom_count)~=1) ;

% Initialise an empty vector to store the e-density:
host_electron_density = zeros(atom_count,max_neighbours_in_cutoff) ;
% Now we should build the host electron-density for each site,
% going FROM each site TO all others (within cutoff):
for from_site = 1:atom_count
    for to_site = 2:max_neighbours_in_cutoff % Start with site two, as first one is self when sorted by distance.
        if angstrom_distances(from_site,to_site) < cutoff
            inspection_distance = angstrom_distances(from_site,to_site) ;
            look_up_index = floor(inspection_distance/dr) + 1 ; % plus one because the look-up starts at zero distance.
            % Store the result of the density in a cumulatove manner:
            host_electron_density(from_site,to_site) = density_function_rho(1,look_up_index) ;
        end
    end
end
host_electron_density = sum(host_electron_density,2) ;
clear from_site
clear to_site

% Now we take these host-electron-densities and lookup the embedding
% function:
% Initialise an empty vector to store the e-density:
embedding_function = zeros(atom_count,1) ;
for site = 1:atom_count
    look_up_index = floor(host_electron_density(site,1)/drho) + 1 ; % plus one because the look-up starts at zero distance.
    embedding_function(site,1) = embedding_function_F(1,look_up_index) ;
end
clear site

% Initialise an empty vector to store the e-density:
pair_potential_term = zeros(atom_count,max_neighbours_in_cutoff) ;
% Now we should build the pair-potential term (within cutoff):
for from_site = 1:atom_count
    for to_site = 2:max_neighbours_in_cutoff % Start with site two, as first one is self when sorted by distance.
        if angstrom_distances(from_site,to_site) < cutoff
            inspection_distance = angstrom_distances(from_site,to_site) ;
            look_up_index = floor(inspection_distance/dr) + 1 ; % plus one because the look-up starts at zero distance.
            % Store the result of the density in a cumulatove manner:
            pair_potential_term(from_site,to_site) = potential_function_phi(1,look_up_index) ;
        end
    end
end
pair_potential_term = sum(pair_potential_term,2) ;

site_energy = embedding_function + 0.5*pair_potential_term ;

%***************************
member_dummy_masses = repmat((count_difference_this_member+count_bank)',[column_count,1]) ;
member_dummy_masses = member_dummy_masses .* column_is_neighbour ;
member_dummy_masses = round(sum(member_dummy_masses,2) ./ sum(member_dummy_masses~=0,2)) ;

if utlise_neighbours == 1
    % Build neighbour PM:
    neighbour_PM = zeros(size(probMatrix_final)) ;
    for column = 1:column_count
        neighbour_PM(:,column) = neighbour_mass_matrix_norm(:,member_dummy_masses(column)) ;
    end
    
    probMatrix_with_neighbours = probMatrix_final .* neighbour_PM ;
    for i = 1:column_count
        probMatrix_with_neighbours(:,i) = probMatrix_with_neighbours(:,i) / sum(probMatrix_with_neighbours(:,i)) ;
    end
end

% subplot(1,3,1)
% imagesc(probMatrix_final)
% title('Observation + Poisson Based PM')
% 
% subplot(1,3,2)
% imagesc(neighbour_PM)
% title('Neighbour Based PM')
% 
% subplot(1,3,3)
% imagesc(probMatrix_with_neighbours)
% title('Combined PM')

% New cost function:

% Loop over the columns and get the global model probability:
colum_probability = zeros(1,column_count) ; % Initialise with proability 1.
for column = 1:column_count
    % Extract probability of each column's atom count and include in running product:
    % FOR "proBMatrix_final" that begins at 0 atom, hence +1 in the index.
    if utlise_neighbours == 1
        colum_probability(column) = probMatrix_with_neighbours((count_difference_this_member(column) + count_bank(column) + 1) , column) ;
    else
        colum_probability(column) = probMatrix_final((count_difference_this_member(column) + count_bank(column) + 1) , column) ;
    end
end

% model_global_probability = prod(colum_probability) ; % Product based overall probability.
% model_global_probability = model_global_probability .^ (1/column_count) ; % Normalise this to reflect the number of columns:

model_global_probability = mean(colum_probability) ; % Average based overall probability.

% 5 columns added for total-energy, total atoms, energy-per-atom, cost-function, and model-probability.
population(member,(2*sequence_length)+1) = sum(site_energy) ;
population(member,(2*sequence_length)+2) = atom_count ;

if compute_skin_only == 1
    if exist('bulk_energy','var') == 0
        bulk_energy = mode(site_energy) ;
    end
    surface_site_energy = site_energy(site_energy > bulk_energy) ;
    population(member,(2*sequence_length)+3) = sum(surface_site_energy) / size(surface_site_energy,1) ;
else
    population(member,(2*sequence_length)+3) = sum(site_energy) / atom_count ;
end

% Create cost-funtion:
population(member,(2*sequence_length)+4) = population(member,(2*sequence_length)+3) * model_global_probability; % Calculate the neergy per atom without new summations.

% Store probability:
population(member,(2*sequence_length)+5) = model_global_probability ;

