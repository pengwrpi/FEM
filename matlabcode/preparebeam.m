clear all;

num_nodes_x = 11;
num_nodes_y = 3;

dx = 10/(num_nodes_x - 1);
dy = 2/(num_nodes_y - 1);

num_nodes = num_nodes_x * num_nodes_y;
num_element = (num_nodes_x - 1) * (num_nodes_y - 1);
zero_dog = zeros(num_nodes_y + 1, 3);
load_element = zeros(num_nodes_y - 1, 4);

coord = zeros(num_nodes, 2);
index_element = zeros(num_element, 6);
for i = 1 : num_nodes
    coord(i, 1) = rem(i - 1, num_nodes_x) * dx;
    coord(i, 2) = floor((i - 1)/num_nodes_x) * dy;      
end

for i = 1 : num_element
    start_index = i + floor((i - 1)/(num_nodes_x - 1));
    index_element(i, 1) = 1;
    index_element(i, 2) = 4;
    index_element(i, 3) = start_index;
    index_element(i, 4) = start_index + 1;
    index_element(i, 5) = start_index + num_nodes_x + 1;
    index_element(i, 6) = start_index + num_nodes_x;
end

for i = 1 : num_nodes_y
    zero_dog(i, 1) = 1 + (i - 1) * num_nodes_x;
    zero_dog(i, 2) = 1;
    zero_dog(i, 3) = 0.0;
end

zero_dog(num_nodes_y + 1, 1) = 1;
zero_dog(num_nodes_y + 1, 2) = 2;
zero_dog(num_nodes_y + 1, 3) = 0.0;

for i = 1 : num_nodes_y - 1
    load_element(i, 1) = i * (num_nodes_x - 1);
    load_element(i, 2) = 2;
    load_element(i, 3) = 0.0;
    load_element(i, 4) = 0.01;
end


fid = fopen('Linear_elastic_beam_40.txt','w');
fprintf(fid, 'No._material_props:    3\n');
fprintf(fid, 'Shear_modulus:   10.\n');
fprintf(fid, 'Poissons_ratio:  0.3\n');
fprintf(fid, 'Plane_strain/stress: 1\n');
fprintf(fid, 'No._coords_per_node:   2\n');
fprintf(fid, 'No._DOF_per_node:      2\n');
fprintf(fid, 'No._nodes:             %d\n', num_nodes);
fprintf(fid, 'Nodal_coords:\n');
fprintf(fid, '%f %f\n', coord');
fprintf(fid, 'No._elements:                       %d\n', num_element);
fprintf(fid, 'Max_no._nodes_on_any_one_element:   4\n');
fprintf(fid, 'element_identifier; no._nodes_on_element; connectivity:\n');
fprintf(fid, '%d %d %d %d %d %d\n', index_element');
fprintf(fid, 'No._nodes_with_prescribed_DOFs:  %d\n', num_nodes_y + 1);
fprintf(fid, 'Node_#, DOF#, Value:\n');
fprintf(fid, '%d %d %f\n', zero_dog');
fprintf(fid, 'No._elements_with_prescribed_loads: %d\n', num_nodes_y - 1);
fprintf(fid, 'Element_#, Face_#, Traction_components\n');
fprintf(fid, '%d %d %f %f\n', load_element');
fclose(fid);