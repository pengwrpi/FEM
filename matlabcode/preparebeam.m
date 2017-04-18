clear all;

num_nodes = 200;
num_element = 99;

coord = zeros(num_nodes, 2);
index_element = zeros(num_element, 6);
load_element = zeros(num_element, 4);
for i = 1 : num_nodes
    coord(i,1) = num_nodes/2 - ceil(i/2);
    coord(i,2) = rem(i,2);
end

for i = 1 : num_element
    max_index = num_nodes - (i - 1) * 2;
    index_element(i, 1) = 1;
    index_element(i, 2) = 4;
    index_element(i, 3) = max_index;
    index_element(i, 4) = max_index - 2;
    index_element(i, 5) = max_index - 3;
    index_element(i, 6) = max_index - 1;
    load_element(i, 1) = i;
    load_element(i, 2) = 3; 
    load_element(i, 3) = 0.0;
    load_element(i, 4) = -1.0;
end

fid = fopen('inputparameter.txt','w');
fprintf(fid, '%f %f\n', coord');
fprintf(fid, '\n');
fprintf(fid, '%d %d %d %d %d %d\n', index_element');
fprintf(fid, '%d %d %f %f\n', load_element');
fclose(fid);