% Check output is not blocking:
if exist('crystal_maker','var') == 1
    clear crystal_maker
end
if exist('filename','var') == 0
    filename = ('filename.mat') ;
end


time_last_save = toc ;
% save( ['3D_reconstruction_progress_at_t_' , num2str(round(time_last_save))] ) ;
% disp( ['Workspace saved to disk at time = ',num2str(round(time_last_save)) , 's.'] )

% Output to Crystal maker:

% Build 'elemental' designations:
for atom = 1:atom_count
    if (atom_neighbour_count(atom,1)<= 1)
        crystal_maker(atom+2,:) = uint16('V ') ;
    elseif (atom_neighbour_count(atom)== 2)
        crystal_maker(atom+2,:) = uint16('Mg') ;
    elseif (atom_neighbour_count(atom)== 3)
        crystal_maker(atom+2,:) = uint16('Au') ;
    elseif (atom_neighbour_count(atom)== 4)
        crystal_maker(atom+2,:) = uint16('Na') ;
    elseif (atom_neighbour_count(atom)== 5)
        crystal_maker(atom+2,:) = uint16('Se') ;
    elseif (atom_neighbour_count(atom)== 6)
        crystal_maker(atom+2,:) = uint16('Zr') ;
    elseif (atom_neighbour_count(atom)== 7)
        crystal_maker(atom+2,:) = uint16('Lu') ;
    elseif (atom_neighbour_count(atom)== 8)
        crystal_maker(atom+2,:) = uint16('Yb') ;
    elseif (atom_neighbour_count(atom)== 9)
        crystal_maker(atom+2,:) = uint16('Al') ;
    elseif (atom_neighbour_count(atom)==10)
        crystal_maker(atom+2,:) = uint16('Np') ;
    elseif (atom_neighbour_count(atom)==11)
        crystal_maker(atom+2,:) = uint16('Ho') ;
    elseif (atom_neighbour_count(atom)>= 12)
        crystal_maker(atom+2,:) = uint16('Co') ;
    end
end
% Restore to strings:
crystal_maker = cellstr(char(crystal_maker)) ;

% Build 'elemental' designations:
for atom = 1:atom_count
    crystal_maker(atom+2,2) = cellstr(num2str(x_cords(atom)*10)) ;
    crystal_maker(atom+2,3) = cellstr(num2str(y_cords(atom)*10)) ;
    crystal_maker(atom+2,4) = cellstr(num2str(z_cords(atom)*10)) ;
    
end

% Record total number of atoms:
crystal_maker(1,1) = cellstr(num2str(atom_count)) ;

% write line-by-line
fid = fopen('Solver_3D_Reconstruction.xyz','wt');
for i=1:size(crystal_maker,1)
    fprintf(fid, '%s %s %s %s\n', crystal_maker{i,:});
end
fclose(fid);


export_string = strsplit(filename , '.' );          % Split filename string into multiple parts...
export_string = export_string(1,1) ;                % and keep first part...
export_string = strcat(export_string , '.xyz') ;    % and create new export filename.

disp('XYZ file can now be opened in CrystalMaker. Be sure to use ''Edit'' > ''Elements'' > ''VFI'' to get correct colours.')
