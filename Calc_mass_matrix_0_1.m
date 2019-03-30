
count_base = repmat(count_bank',[population_size,1]) ;
all_col_masses = population(:,(sequence_length+1):(2*sequence_length)) + count_base ;
neighbour_mass_matrix = zeros(size(probMatrix_final,1),size(probMatrix_final,1)) ;

for mass_member = 1:breeding_cutoff
    member_mass = all_col_masses(mass_member,:) ;
    
    for i = 1:column_count
       my_mass = member_mass(i) ;
       
       for j = 1:column_count
        if column_is_neighbour(i,j) == 1
           neighbour_mass = member_mass(j) ;
           neighbour_mass_matrix(my_mass+1,neighbour_mass+1) = neighbour_mass_matrix(my_mass+1,neighbour_mass+1) + 1 ; % Coordinates are both '+1' so that the resulting matrix is (0,0) origined.
        end
       end
    end
        
end

% Now normalise matric ready for combining with the prob_Matrix:
neighbour_mass_matrix_norm = zeros(size(neighbour_mass_matrix)) ;
for i = 1:size(neighbour_mass_matrix,2)
   f = fit([1:size(neighbour_mass_matrix,1)]',neighbour_mass_matrix(:,i),'gauss1') ;
   neighbour_mass_matrix_norm(:,i) = f([1:size(neighbour_mass_matrix,1)]) ;
   neighbour_mass_matrix_norm(:,i) = neighbour_mass_matrix_norm(:,i) / sum(neighbour_mass_matrix_norm(:,i)) ; 
end

% subplot(1,2,1)
% imagesc(neighbour_mass_matrix)
% axis image
% ylabel('Inspection Column Mass')
% xlabel('Neighbouring Masses')
% xticks([1:size(probMatrix_final,1)])
% xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18'})
% yticks([1:size(probMatrix_final,1)])
% yticklabels({'0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18'})
% colorbar
% 
% subplot(1,2,2)
% imagesc(neighbour_mass_matrix_norm)
% axis image
% ylabel('Inspection Column Mass')
% xlabel('Neighbouring Masses')
% xticks([1:size(probMatrix_final,1)])
% xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18'})
% yticks([1:size(probMatrix_final,1)])
% yticklabels({'0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18'})
% colorbar

clear all_col_masses
clear member_mass
clear mass_member
