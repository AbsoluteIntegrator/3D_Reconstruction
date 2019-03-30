clear all
close all
clc

F = findall(0,'type','figure','tag','TMWWaitbar'); % Find all waitbars...
delete(F); % Close all waitbars.
% clear
format compact

%****************************

element = 'Pt' ; % Choose the element, either 'Au' or 'Pt'.
radius_truncation = 2.5 ; % Typicaly value is 2.5

% Display options:
front_angle = -111 ;
lift_angle = 31 ;
max_atoms_to_draw = 5000 ;

% Initialisation Options:
vectorise_x_y = 1 ; % Move lateral coordinates onto nearest clean base vectors.
selection_multiples = 4 ;
use_analytic_mode  =  1 ;  %

% Vibration range max values (multiples of R0):
lat_range = 0.0 ; % Standard deviations of the normal distribution, units of fractional inter-atomic separation.

save_outputs = 1 ; % Indicate the model to save (within the first 20).
population_size = 500 ;
convergence_population_size = 500 ;
unique_initialisations = 25 ;
view_every = 10 ;
grab_frame = 0 ;

% target_extinction_events = 1 ;
height_mutation_range = 1 ;
height_mutation_step = 1 ; % 1 for single crystal, 0.5 if twinned.
count_mutation_range = 1 ;
mutation_density = 0.02 ;
breeding_fraction = 0.5 ;
use_EAM = 1 ;
utlise_neighbours = 1 ;
use_blockwise_breeding = 0 ; % Interbreed genome in blocks not random-parent. '1' for blockwise in x, '2' for y, '0' for random gene-mix.
compute_skin_only = 0 ;
%------------------------- End User Variables -----------------------------
%--------------------------------------------------------------------------
stagnation_period = 100 ; % How long the inbreeding is allowed to stagnate.

% Crystallographic parameters:
if strcmp(element,'Pt') == 1
% % Platinum:
R0 = 0.27748 ; % Equilibrium separation (2x VdW radius). Units of nm.
Rz = 0.27748 ;%* realsqrt(2) ; % Projected separation in z. Units of nm.
% Lennard-Jones paramerters (from Avinc et al.):
LJ_a = 9 ;
LJ_b = 6.52648 ;
LJ_A = 0.11564 ;
LJ_B = -0.14805 ;
% EAM potential data:
load('C:\Users\Lewys\Box Sync\PostDoc\Structure Solving\Pt_EAM_Input_data.mat')
elseif strcmp(element,'Au') == 1
% Gold:
R0 = 0.28837 ; % Equilibrium separation (2x VdW radius). Units of nm.
Rz = 0.28837 ;%* realsqrt(2) ; % Projected separation in z. Units of nm.
% Lennard-Jones paramerters (from Avinc et al.):
LJ_a = 9 ;
LJ_b = 6.54223 ;
LJ_A = 0.10173 ;
LJ_B = -0.13003 ;
% EAM potential data:
load('C:\Users\Lewys\Box Sync\PostDoc\Structure Solving\Au_EAM_Input_data.mat')
else
    error('Invalid element selected')
end
% Predivide potentials:
potential_function_phi = potential_function_phi ./ radial_grid ;


run 'C:\Users\Lewys\Box Sync\DPhil\MatLab Software\File Opener\File_Opener_1_5_6'

% The following two lines create inputs from Aakash's files:
% point_coordinates = position_and_thickness(:,1:2) ;
% count = round(position_and_thickness(:,3)) ;

% The following three lines create inputs from Annick's files:
point_coordinates = coordinates ;
clear coordinates
count = atomCounts ;
clear atomCounts
pixelWidth = 0.1 ;

probMatrix_final_bigger = zeros( size(probMatrix_final,1)+4 , size(probMatrix_final,2) ) ;
probMatrix_final_bigger(1:size(probMatrix_final,1),:) = probMatrix_final ;
probMatrix_final = probMatrix_final_bigger ;
clear probMatrix_final_bigger

if use_blockwise_breeding > 0
    % Append coordinates to probability matrix:
    probMatrix_final_for_sorting = [probMatrix_final',point_coordinates,count] ;
    probMatrix_final_for_sorting = (sortrows(probMatrix_final_for_sorting,(size(probMatrix_final,1)+use_blockwise_breeding)))' ;
    % Extract sorted versions:
    probMatrix_final = probMatrix_final_for_sorting(1:size(probMatrix_final,1),:) ;
    point_coordinates = probMatrix_final_for_sorting(size(probMatrix_final,1)+1:size(probMatrix_final,1)+2,:)' ;
    count = probMatrix_final_for_sorting(size(probMatrix_final,1)+3,:)' ;
    clear probMatrix_final_for_sorting
end

% Intorduce extra count uncertainty:
figure
subplot(2,1,1)
imagesc(probMatrix_final(1:20,:))
xlabel('Column Number')
ylabel('Library')
yticks([1:size(probMatrix_final,1)])
yticklabels({'0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20'})
title('Original Uncertainty')
colorbar

% H = fspecial('gaussian',[15 1],2) ;
% probMatrix_final = conv2(probMatrix_final,H,'same');

% Verify that the probabilities of each column remain equal to 1.
for i=1:size(probMatrix_final,2)
    probMatrix_final(:,i) = probMatrix_final(:,i) / sum(probMatrix_final(:,i)) ;
end

subplot(2,1,2)
imagesc(probMatrix_final)
xlabel('Column Number')
ylabel('Library')
yticks(1:size(probMatrix_final,1))
yticklabels({'0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18'})
title('Extra Uncertainty')
colorbar
set(gcf,'color','w');

figure

original_sum = sum(count) ;
count_bank = count ;

if vectorise_x_y == 1
    rotate_coordinate_set = 0 ;
    manual_centre_point = 0 ;
    run Refine_points_0_13
end

close all

% Count the number of coloumns:
column_count = size(point_coordinates,1) ;

% Solve the number of atoms in each column (improve this method later):
if exist('count' , 'var') == 0
    % Blank variable to hold the number of atoms per column:
    count = zeros(1,column_count) ;
    for column = 1:column_count
        similarities = abs( library - integratedIntensity(column) ) ;
        [~, I] = min(similarities) ;
        count(column) = I ;
        %         columnwise_residual(column) = integratedIntensity(column) - library(I) ;
        %         percentage_residual(column) = ( integratedIntensity(column) - library(I) ) / (library(I+1)-library(I)) ;
    end
end

% Feedback to user:
atom_count = sum(count) ;
disp(['Atom counting assignments complete. Total number of atoms = ', num2str(atom_count),'.'])
disp(['Greatest number of atoms in any column = ', num2str(max(count)),'.'])
disp('Starting GA now...')
tic;
% Generate blank variables for the atom coordinates:
x_offsets = zeros(1,column_count) ;
y_offsets = zeros(1,column_count) ;
z_offsets = zeros(1,column_count) ;

% Now calculate starting properties in a columnwise sense:
atom_neighbour_count = zeros(atom_count,1) ;
column_energy = zeros(1,column_count) ;

D = squareform(pdist(point_coordinates)) ;
column_is_neighbour = D<=(R0*10*1.2) ;
column_is_neighbour = column_is_neighbour - eye(size(column_is_neighbour)) ;

clear D
%****** HERE ENDS THE LINES COMMON TO ALL POSIBLE INITIALISATIONS *****
figure(1)
figure(2)
figure(3)

% Gives the starting estimate of the absolute heights to bring all to
% midplane. This is the origin of all in the population.
z_offsets_bank = -round(0.5*count') ;

% Create starting 3D coordinates:
if use_analytic_mode == 1
    if exist('analytic_heights','var') == 1
        z_offsets = z_offsets_bank + analytic_heights' ; % Adds an estimated analytic solution.
    end
else
    z_offsets = z_offsets_bank + ( rand(1,size(count,2)) ) ; % Gives the starting estimate of the absolute heights to bring all to midplane BUT also adds a random additional offset.
end

%% Variable setup:
tic
generation = 1 ;
% Setup variables for the genetic algorthm:
breeding_cutoff = floor(breeding_fraction*population_size) ; % an integer.
sequence_length = column_count ;
initialisations_converged = 0 ;
last_extinction_event = 0 ;
calculate_from = breeding_cutoff + 1 ;
% Initialise master matrix and append data saving columns:
population = zeros(population_size,(sequence_length*2) + 4) ; % 4 columns added for total-energy, total-negihbours, total atoms, energy-per-atom.
% Metadata storage variables:
energy_log = 0 ;
best_mass_log = 0 ;
energy_per_atom_log = 0 ;
cost_function_log = 0 ;
model_probability_log = 0 ;
inbreeding_log = 0 ;
%% GA setup:
for member = 1:breeding_cutoff % Create starting population for breeding:
    population(member,1:sequence_length) = z_offsets + (rand(1,sequence_length) <= 2*mutation_density).*(height_mutation_step*randi([-height_mutation_range/height_mutation_step,height_mutation_range/height_mutation_step],sequence_length,1)') ;
    count = round( count_bank + (rand(1,sequence_length) <= 2*mutation_density)'.*randi([-count_mutation_range,count_mutation_range],sequence_length,1)) ;
    % Check the proposed counts to be valid:
    count(count<0) = 0 ; % Cannot allow columns to become empty.
    count_difference_this_member = count - count_bank ;
    population(member,(sequence_length+1):(2*sequence_length)) = count_difference_this_member' ;
end
%% GA Loop:
while initialisations_converged <= unique_initialisations
    %*********** New solution generation step:
    
    for member = calculate_from:population_size
        
        % Define parent responsible for each trait:
        if use_blockwise_breeding > 0 % Genes from parents in contiguous blocks:
            trait_from = zeros(1,sequence_length) ; % Reset clean traits.
            second_parent_from = randi(sequence_length-1)+1 ;
            trait_from(second_parent_from:end) = 1 ;
            trait_from = repmat(trait_from,1,2)-1 ;
        else % Genes randomly intermixed:
            trait_from = randi(2,1,sequence_length) ; % Equal genes from both parents.
            trait_from = repmat(trait_from,1,2)-1 ;
        end
        
        % Identify parent pairs:
        if initialisations_converged < unique_initialisations % Initialisations not complete yet:
            first_parent = randi([(initialisations_converged+1) (initialisations_converged+breeding_cutoff)]) ;
            second_parent = randi([(initialisations_converged+1) (initialisations_converged+breeding_cutoff)]) ;
        else % Initialisations complete - final supermodel convergence:
            first_parent = randi([1 (initialisations_converged+breeding_cutoff)]) ;
            second_parent = randi([1 (initialisations_converged+breeding_cutoff)]) ;
        end
        
        % Pure breeding:
        population(member,1:(2*sequence_length)) = ...
            population(first_parent,1:(2*sequence_length)) .* trait_from + ...
            population(second_parent,1:(2*sequence_length)) .* (1-trait_from) ;
        
        % Then introduce mutation for heights:
        population(member,1:sequence_length) = population(member,1:sequence_length) + (rand(1,sequence_length) <= mutation_density).*(height_mutation_step*randi([-height_mutation_range/height_mutation_step,height_mutation_range/height_mutation_step],1,sequence_length)) ;
        % and for atoms counts:
        count_difference_this_member = (population(member,(sequence_length+1):(2*sequence_length)) + (rand(1,sequence_length) <= mutation_density).*randi([-count_mutation_range count_mutation_range],1,(sequence_length)))' ;
        
        % Check the proposed counts to be valid:
        count = count_difference_this_member + count_bank ; % This is what the occupancies will be
        count(count<0) = 0 ; % Cannot allow columns to become negative. Empty is NOT specifically prohibited.
        count_difference_this_member = count - count_bank ; % Rephrase in terms of differences and store...
        population(member,(sequence_length+1):(2*sequence_length)) = count_difference_this_member' ;
    end
    
    if utlise_neighbours == 1
        % Calculate here the neighbour mass matrix probabilities afresh:
        if generation == 1 % mod(generation,50) == 1
            run Calc_mass_matrix_0_1
            figure;imagesc(neighbour_mass_matrix_norm);
        end
    end
    
    %*********** Check for duplicates here:
    unique_population = unique(population(calculate_from:end,1:2*sequence_length), 'rows');
    unique_children = size(unique_population,1) ;
    calculate_to = calculate_from+unique_children-1 ;
    population(calculate_from:calculate_to,1:2*sequence_length) = unique_population ; % Fill with unique solutions
    clear unique_population
    
    %*********** Evaluation step:
    for member = calculate_from:calculate_to % The evalutaion will only go as far as the unique models, duplicates are not evaluated.
        % Extract relevant data and evaulate:
        count_difference_this_member = population(member,(sequence_length+1):(2*sequence_length))' ; % used internally.
        if use_EAM == 1
            run build_points_and_evaluate_EAM.m
        else
            run build_points_and_evaluate_LJ.m
        end
    end
    
    % Sort the best breeders by cost-function, or in a tie, then by energy...
    population = sortrows(population,[(2*sequence_length)+4 (2*sequence_length)+3]) ; 
    
    
    unique_parents = size(unique(population(initialisations_converged + 1:initialisations_converged+breeding_cutoff-1,:), 'rows'),1) ;
    % Record metadata:
    energy_log(generation)            = population(initialisations_converged + 1,((2*sequence_length)+1)) ; % Record the lowest energy.
    best_mass_log(generation)         = population(initialisations_converged + 1,((2*sequence_length)+2)) ; % Record total number of atoms.
    energy_per_atom_log(generation)   = population(initialisations_converged + 1,((2*sequence_length)+3)) ; % Record the energy per atom.
    cost_function_log(generation)     = population(initialisations_converged + 1,((2*sequence_length)+4)) ; % Record the cost function.
    model_probability_log(generation) = population(initialisations_converged + 1,((2*sequence_length)+5)) ; % Record the model probability.
    inbreeding_log(generation)        = unique_parents / breeding_cutoff ; % Record percentage unique members.
    child_diversity(generation)       = unique_children / (population_size+1-calculate_from) ;
    
    if mod(generation,view_every) == 0
        run GA_tracking_figures
    end
    
    if grab_frame == 1
       run Calc_points_and_coordination 
    end
    
    % Initialisation tracking window:
    if generation < stagnation_period
        rolling_start = 1 ;
    else
        rolling_start = generation - stagnation_period + 1 ;
    end
    
    % New initialistion condition:
    if (std(cost_function_log(rolling_start:generation))/-cost_function_log(generation) < 10^-5) && (initialisations_converged < unique_initialisations) && (generation > last_extinction_event+stagnation_period)
        
        % Incriment initialisation counter:
        initialisations_converged = initialisations_converged + 1 ;
        last_extinction_event = generation ;
        disp(['Initialisations completed = ',num2str(initialisations_converged),'.'])
        
        % Ensure active models are sorted before calculation start point is
        % incrimented and increase population siz by one:
        population(initialisations_converged:end) = sortrows(population(initialisations_converged:end),((2*sequence_length)+4)) ; % Sort the best breeders by energy...
        population_size = population_size + 1 ;
        supermodel_CF(initialisations_converged) = population(initialisations_converged,((2*sequence_length)+4)) ;
        % Preserve unique models by shifting calculation start threshold:
        calculate_from = initialisations_converged + 1 ; % Trigger a one off reinitialisation.
        
        figure(3)
        hist(supermodel_CF,floor(sqrt(population_size)),'normal')
        pause(0.1)
        title(['Supermodel Diversity. Min = ',num2str(min(supermodel_CF)),'eV. Mean = ',num2str(mean(supermodel_CF))])
        xlabel('Cost Function (ev)')
        ylabel('Frequency')
        pause(0.1)
        
        if initialisations_converged == unique_initialisations
            population(convergence_population_size+1:end,:) = [] ;
            population_size = convergence_population_size ;
%             stagnation_period = population_size ;
            breeding_cutoff = floor(breeding_fraction*population_size) ; % Update to reflect the larger population.
            calculate_from = breeding_cutoff + 1 ; % Otherwise only reinitialise from beyond breeders.
        end
        
        % Create new breeding population (erases previous breeders):
        for member = initialisations_converged+1:initialisations_converged+breeding_cutoff % Create starting population for breeding:
            population(member,1:sequence_length) = z_offsets + (rand(1,sequence_length) <= 2*mutation_density).*(height_mutation_step*randi([-height_mutation_range/height_mutation_step,height_mutation_range/height_mutation_step],sequence_length,1)') ;
            count = round( count_bank + (rand(1,sequence_length) <= 2*mutation_density)'.*randi([-count_mutation_range,count_mutation_range],sequence_length,1)) ;
            % Check the proposed counts to be valid:
            count(count<0) = 0 ; % Cannot allow columns to become empty.
            count_difference_this_member = count - count_bank ;
            population(member,(sequence_length+1):(2*sequence_length)) = count_difference_this_member ;
        end
        
    elseif initialisations_converged == unique_initialisations
        calculate_from = breeding_cutoff + 1 ; % Otherwise only reinitialise from beyond breeders.
    else
        calculate_from = breeding_cutoff + initialisations_converged + 1 ; % Otherwise only reinitialise from beyond breeders.
    end
    
    % This line will break us out of the GA when the supermodels have
    % converged:
    if (std(cost_function_log(rolling_start:generation))/-cost_function_log(generation) < 10^-5) && ...
            (initialisations_converged == unique_initialisations) && (generation > last_extinction_event+stagnation_period )
        run GA_tracking_figures
        breeding_cutoff = floor(breeding_fraction*population_size) ; % an integer.
        break
    end
    
    % Loop updating:
    generation = generation + 1 ;
    
end

toc

speed_metric = toc / 3600 /column_count / population_size / (initialisations_converged+1) ; % +1 because of the pseudo-initialisation with the super-models.
disp(['CPUh/(col*pop*init) = ',num2str(speed_metric),'.'])
speedup = round((10000/692/60) / speed_metric) ;
disp(['Speedup over Voyles = ',num2str(speedup),'x.'])
disp(['Total generations used = ',num2str(generation),'.'])
% Output family:
% Find unique population:
unique_population = unique(population, 'rows');
% Sort for best cost-function.
unique_population = sortrows(unique_population,(2*sequence_length+4)) ;
count_difference_this_member = unique_population(1,(sequence_length+1):(2*sequence_length))' ;
disp(['Best model Cost-function = ',num2str(unique_population(1,((2*sequence_length)+4))),'eV.'])

figure(4)
clf
subplot(1,2,1)
imagesc(unique_population(:,1:sequence_length))
title('Height Offset Map')
xlabel('Columns')
ylabel('Population')
colorbar

subplot(1,2,2)
imagesc(unique_population(:,(sequence_length+1):(2*sequence_length)))
title('Atom Difference Map')
xlabel('Columns')
ylabel('Population')
colorbar

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
        z_cords(atom_count) = (i + unique_population(1,column) ) * Rz ;
        
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


% Plot final annotated probability matrix:
figure(5)
clf
final_assignments = count_difference_this_member + count_bank ;
imagesc(probMatrix_final)
xlabel('Column Number')
ylabel('Library')
yticks([1:size(probMatrix_final,1)])
yticklabels({'0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18'})
title('Final Assignments')
colorbar
hold on
plot([1:column_count],final_assignments+1,'o','MarkerEdgeColor',[0.8,0,0],'MarkerFaceColor',[0.7,0,0]) % +1 shifts to include the zero.
hold off


% Plot the two models:
load('C:\Users\Lewys\Box Sync\PostDoc\Structure Solving\Test Object Simulation\New Embedded Test Object\Phantom_coordinates.mat')
% Centre models at the origin:
x_cords_phantom = x_cords_phantom - mean(x_cords_phantom) ;
y_cords_phantom = y_cords_phantom - mean(y_cords_phantom) ;
z_cords_phantom = z_cords_phantom - mean(z_cords_phantom) ;
x_cords = x_cords - mean(x_cords) ;
y_cords = y_cords - mean(y_cords) ;
z_cords = z_cords - mean(z_cords) ;

figure(6)
clf
subplot(1,2,1)
scatter3(x_cords_phantom,y_cords_phantom,z_cords_phantom,'filled')
hold on
scatter3(x_cords,y_cords,z_cords,'filled')
hold off
axis equal
axis vis3d


subplot(1,2,2)
edges = [0.5:1:12.5];
h1 = histogram(atom_neighbour_count,edges);
xlim([0.5 12.5])
grid minor
hold on
h2 = histogram(atom_neighbour_count_phantom,edges);
hold off

figure(7)
clf
for show = 1:20
    figure(7)
    subplot(4,5,show)
    
    member = show ;
    count_difference_this_member = unique_population(member,(sequence_length+1):(2*sequence_length))' ;
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
            z_cords(atom_count) = (i + unique_population(member,column) ) * Rz ;
            
        end
    end
    
    %     z_cords = -z_cords ;
    
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
    
    % Save final results:
    if member == save_outputs
        run Crystal_Maker_Export_1_0.m
    end
    
    run('C:\Users\Lewys\Box Sync\DPhil\MatLab Software\Absolute Integrator\3D Reconstruction\SingleDraw.m')
    pause(0.5)
end


