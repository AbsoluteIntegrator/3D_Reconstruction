my_color_spec = [2/3 0 0 ; 1 0 0 ; 1 1/3 0 ; 1 2/3 0 ; 1 1 0 ; 2/3 1 1/3 ; 1/3 1 2/3 ; 0 1 1 ; 0 2/3 1 ; 0  1/3  1 ; 0 0 1 ; 0 0 2/3] ;
    
% Generate atom colour data:
atom_colors = (atom_neighbour_count<= 1) * my_color_spec(1,:) + ...
    (atom_neighbour_count== 2) * my_color_spec(2,:) + ...
    (atom_neighbour_count== 3) * my_color_spec(3,:) + ...
    (atom_neighbour_count== 4) * my_color_spec(4,:) + ...
    (atom_neighbour_count== 5) * my_color_spec(5,:) + ...
    (atom_neighbour_count== 6) * my_color_spec(6,:) + ...
    (atom_neighbour_count== 7) * my_color_spec(7,:) + ...
    (atom_neighbour_count== 8) * my_color_spec(8,:) + ...
    (atom_neighbour_count== 9) * my_color_spec(9,:) + ...
    (atom_neighbour_count==10) * my_color_spec(10,:) + ...
    (atom_neighbour_count==11) * my_color_spec(11,:) + ...
    (atom_neighbour_count>= 12) * my_color_spec(12,:)   ; % Neighbour numbers include the self-neighbour so are one higher.


if atom_count < max_atoms_to_draw
    scatter3sph(x_cords, y_cords,z_cords , 'size' ,Rz/(2*1) ,'color',atom_colors) % Outside function - requires time to render.
else
    plot3(x_cords, y_cords,z_cords ,'.') % Internal function - faster.
end
axis equal
axis vis3d
view(front_angle,lift_angle);
xlabel('X (nm)','FontSize',11,'FontWeight','bold')
ylabel('Y (nm)','FontSize',11,'FontWeight','bold')
zlabel('Z (nm)','FontSize',11,'FontWeight','bold')
grid on