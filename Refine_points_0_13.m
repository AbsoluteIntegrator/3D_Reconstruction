
% Setup variables:
if exist('selection_multiples' , 'var') == 0
    selection_multiples = 4
end

if exist('coordinates' , 'var') == 1
    point_coordinates = coordinates ;
end

if exist('manual_centre_point' , 'var') == 0
    manual_centre_point = 0 ; % Choose whether to let the code choose the origin.
end

if exist('input','var') == 1
    clear input
end

rattle_range = 0.02 ; % Set here the rattling range.
max_rattle_tries = 100000 ;
rotate_coordinate_set = 0 ;
lock_vectors   =  0 ;

% ****** END setup variables *********

if exist('zoom_for_selection','var') == 0
    zoom_for_selection = 1
end

if rotate_coordinate_set ~= 0
    point_coordinates = point_coordinates*[cosd(rotate_coordinate_set) -sind(rotate_coordinate_set); sind(rotate_coordinate_set) cosd(rotate_coordinate_set)] ;
end

%  point_coordinates = point_coordinates*[ 1 -0.03 ; 0.085 1 ] ;

% imagesc(input_file)
% hold on
plot(point_coordinates(:,2) , point_coordinates(:,1),'ro')
set(gca,'YDir','reverse');
axis image
grid on
hold off


acceptable_solution = 'n' ;
while strcmp(acceptable_solution,'y') == 0
    
    if manual_centre_point == 1
        [home_point(2) ,home_point(1)] = getpts ;
    else
        % Following line for auto centre of mass:
        home_point = mean(point_coordinates,1) ; % Find centroid of mass:
    end
    distance_to_points = sum( ( bsxfun(@minus,point_coordinates,home_point) ) .^2 , 2 ) ; % Find closest point to start at:
    start_node = find(distance_to_points==min(distance_to_points),1,'first') ;
    temp_origin = point_coordinates(start_node,:) ;
    dummy_point_coordinates = ( point_coordinates - repmat(temp_origin,[size(point_coordinates,1) 1]) )' ;
    
    % Graphically choose approximate vectors if undefined at start:
    if exist('manual_vectors' , 'var') == 0
        base_select = figure ;
        plot(dummy_point_coordinates(1,:) , dummy_point_coordinates(2,:),'o')
        hold on
        plot(0,0,'ro')
        set(gca,'YDir','reverse');
        axis image
        grid on
        hold off
        zoom(zoom_for_selection)
        
        % Click nearby points for base vector estimate:
        title('Choose first base vector','FontSize',13,'FontWeight','bold')
        [base_a(1) ,base_a(2)]  = getpts ;
        title('Choose second base vector','FontSize',13,'FontWeight','bold')
        [base_b(1) ,base_b(2)]  = getpts ;
        
        manual_vectors = [ base_a' , base_b' ]  / selection_multiples ;
        close(base_select)
        
    end
    
    % This creates a matrix with the form:
    % manual_vectors = [ a_x , b_x ]
    %                  [ a_y , b_y ]
    
    values = (manual_vectors \ dummy_point_coordinates)' ;
    
    subplot(1,3,3)
    plot(values(:,1) , values(:,2),'bo')
    axis equal
    set(gca,'YDir','reverse');
    grid on
    title('Solution in New Basis Set','FontSize',13,'FontWeight','bold')
    
    % Guess some vectors:
    vectors = 25*randn(2) ;
    
    % Initialise search. Start with knonw best values:
    best_vectors = manual_vectors ;
    least_error = 10e10 ;
    
    values = vectors \ dummy_point_coordinates ;
    restored_points_trial = vectors * round(values) ;
    differences = restored_points_trial - dummy_point_coordinates ;
    %     errors = sum(abs(differences(:))) ;
    errors = sum((differences(:)).^2) ;
    
    % Start with knonw best values
    best_vectors = manual_vectors ;
    h = waitbar(0,'Refining positions...') ;
    
    for updates = 1:max_rattle_tries % Do not change. This is the number of times the progress will be updated.
        
        % Try a new guess (both guessed together):
        vectors = best_vectors + rattle_range*randn(2) ;
        
        %%% Special lines for locked ratio & angle:
        if lock_vectors == 1
            vectors(:,2) = flipud(vectors(:,1))/sqrt(2) ;
            vectors(1,2) = -vectors(1,2) ;
        end
        %%% End special lines.
        
        values = (vectors \ dummy_point_coordinates)' ;
        
        restored_points_trial = (round(values) * vectors')' ;
        
        % Offsets:
        differences = restored_points_trial - dummy_point_coordinates ;
        
%         errors = sum(abs(differences(:)).^0.5)  ;
        errors = sum(abs(differences(:)))  ;
        
        if errors < least_error
%             best_vectors = vectors ;
            % Try only accepting half the change??
            best_vectors = 0.5 * ( best_vectors + vectors ) ;
            best_restored_points_trial = restored_points_trial ;
            best_differences = differences ;
            least_error = errors ;
            
            subplot(1,3,1)
            plot(dummy_point_coordinates(1,:) , dummy_point_coordinates(2,:),'r*')
            axis image
            set(gca,'YDir','reverse');
            grid on
            title('Experimental Observations','FontSize',13,'FontWeight','bold')
            hold on
            plot(restored_points_trial(1,:),restored_points_trial(2,:),'bo')
            line(selection_multiples*[0 manual_vectors(1,1)],selection_multiples*[0 manual_vectors(2,1)],'LineWidth',4,'Color',[.2 .2 .2])
            line(selection_multiples*[0 manual_vectors(1,2)],selection_multiples*[0 manual_vectors(2,2)],'LineWidth',4,'Color',[.2 .2 .2])
            line(selection_multiples*[0 best_vectors(1,1)],selection_multiples*[0 best_vectors(2,1)],'LineWidth',2,'Color',[.8 .8 .8])
            line(selection_multiples*[0 best_vectors(1,2)],selection_multiples*[0 best_vectors(2,2)],'LineWidth',2,'Color',[.8 .8 .8])
            hold off
            axis image
            
            subplot(1,3,2)
            quiver(dummy_point_coordinates(1,:) , dummy_point_coordinates(2,:),-best_differences(1,:),-best_differences(2,:))
            axis image
            set(gca,'YDir','reverse');
            title('Offsets from Grid Points','FontSize',13,'FontWeight','bold')
            
            subplot(1,3,3)
            plot(values(:,1) , values(:,2),'r*')
            hold on
            plot(round(values(:,1)) , round(values(:,2)),'bo')
            hold off
            axis image
            set(gca,'YDir','reverse');
            grid on
            legend('Measured','Gridded','Location','Best')
            title('Solution in New Basis Set','FontSize',13,'FontWeight','bold')
            
            
            pause(0.01)
        end
        
        if mod(updates-1,100) == 0 ;
            temp_origin = temp_origin - mean(best_differences,2)' ;  % Steps the tempoary origin in the direction of the global average residual.
            dummy_point_coordinates = ( point_coordinates - repmat(temp_origin,[size(point_coordinates,1) 1]) )' ;
            rattle_range = ( rattle_range * 0.99 ) + 0.00001 ;
            waitbar(updates/max_rattle_tries)
        end
        
    end
    acceptable_solution = input('Is solution acceptable (y/n)?','s') ;
    
    delete(h)
end

% Cleanup:
clear rattle_range
clear vectors

values = (best_vectors \ dummy_point_coordinates)' ;

% Now we have an acceptable solution for the values, these can be solved to
% unambiguously yiled the final vectors:
best_vectors = (values \ dummy_point_coordinates')' ;
% lattice_parameter = realsqrt(sum(best_vectors(:,2).^2))*sqrt(2) 
lattice_parameter = 0.5 * (realsqrt(sum(best_vectors(:,2).^2)) + realsqrt(sum(best_vectors(:,1).^2)) )

values = (best_vectors \ dummy_point_coordinates)' ;
restored_points_trial = (round(values) * best_vectors')' ;
% Offsets:
differences = restored_points_trial - dummy_point_coordinates ;
errors = sum(abs(differences(:)).^0.5)  ;

% Final Replotting:
subplot(1,3,1)
plot(dummy_point_coordinates(1,:) , dummy_point_coordinates(2,:),'r*')
axis image
set(gca,'YDir','reverse');
grid on
title('Experimental Observations','FontSize',13,'FontWeight','bold')
hold on
plot(restored_points_trial(1,:),restored_points_trial(2,:),'bo')
line(selection_multiples*[0 manual_vectors(1,1)],selection_multiples*[0 manual_vectors(2,1)],'LineWidth',4,'Color',[.2 .2 .2])
line(selection_multiples*[0 manual_vectors(1,2)],selection_multiples*[0 manual_vectors(2,2)],'LineWidth',4,'Color',[.2 .2 .2])
line(selection_multiples*[0 best_vectors(1,1)],selection_multiples*[0 best_vectors(2,1)],'LineWidth',2,'Color',[.8 .8 .8])
line(selection_multiples*[0 best_vectors(1,2)],selection_multiples*[0 best_vectors(2,2)],'LineWidth',2,'Color',[.8 .8 .8])
hold off
axis image

subplot(1,3,2)
quiver(dummy_point_coordinates(1,:) , dummy_point_coordinates(2,:),-best_differences(1,:),-best_differences(2,:))
axis image
set(gca,'YDir','reverse');
title('Offsets from Grid Points','FontSize',13,'FontWeight','bold')

subplot(1,3,3)
plot(values(:,1) , values(:,2),'r*')
hold on
plot(round(values(:,1)) , round(values(:,2)),'bo')
hold off
axis image
set(gca,'YDir','reverse');
grid on
legend('Measured','Gridded','Location','Best')
title('Solution in New Basis Set','FontSize',13,'FontWeight','bold')
% End of replotting.

    
% figure
% imagesc(input_file)
% colormap gray
% axis image
% axis off
% hold on
% quiver( best_restored_points_trial(1,:)+temp_origin(1),best_restored_points_trial(2,:)+temp_origin(2) ,-best_differences(1,:),-best_differences(2,:),'r')
% hold off


% Compute analytical height forms:

% [100]:
% Here the
% analytic_heights = 0.5 *  mod( round(values(:,1)) + round(values(:,2)),2) ;

% [110]:
% Here the first vector selected by the user will step up by one half unit in height:
analytic_heights = 0.5 * mod(round(values(:,1)),2) ;

if size(unique(best_restored_points_trial,'rows'),2) == size(point_coordinates,1)
    % Pass back coordinates to Solver:
    point_coordinates = best_restored_points_trial' + repmat(temp_origin,[size(point_coordinates,1) 1]) ;
    disp('Atom coordinates updated to grid points.')
else
    warning('Poor solution')
end

if rotate_coordinate_set == 1 % Select this to rotate all poitns globally
    figure
    % imagesc(input_file)
    hold on
    plot(point_coordinates(:,1) , point_coordinates(:,2),'ro')
    set(gca,'YDir','reverse');
    axis image
    grid on
    
    title('Select start point for scan direction...')
    [scan_start(1) , scan_start(2)]  = getpts ;
    plot(scan_start(1) , scan_start(2) , 'b*')
    
    title('Select end point for scan direction...')
    [scan_end(1) , scan_end(2)]  = getpts ;
    plot(scan_end(1) , scan_end(2) , 'b*')
    
    hold off
    
    scan_angle = atand ( ( scan_end(2) - scan_start(2) ) / ( scan_end(1) - scan_start(1) ) ) ;
    
    best_vectors_rotated = best_vectors' * [cosd(scan_angle) -sind(scan_angle); sind(scan_angle) cosd(scan_angle)] ;
    
    % Modify the 'perfect' vetors here to obey the crystallography:
    warning('Lattice vectors manual overide engaged!')
    best_vectors_rotated(2,1) = 0 ;
    best_vectors_rotated(1,:) =  best_vectors_rotated(1,:) * ((0.5*best_vectors_rotated(2,2))/-best_vectors_rotated(1,2)) ;
    
    point_coordinates = (round(values) * best_vectors_rotated) ;
    
    plot(point_coordinates(:,1) , point_coordinates(:,2),'ro')
    set(gca,'YDir','reverse');
    axis image
    grid on
    title('Points after scan rotation.')
    
end

% clear best_vectors
clear best_vectors_rotated
clear manual_vectors
