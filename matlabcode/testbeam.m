clear all;
num_node = 200;
fid = fopen('beam_displacement.txt','r');
coord = fscanf(fid, '%d %f %f', [3, num_node]);
fclose(fid);

plot((num_node - 2)/2:-1:0, coord(3, 1:2:num_node), 'o')