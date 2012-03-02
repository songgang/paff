function mask = scatter_binary_trajectory_nearest_neighbor(gbox, cps)
% cps : nb_time_points * nb_dim * nb_trajectories


% reshape cps to a 2D array: nb_points * nb_dim
dim = size(cps, 2);
cps1 = permute(cps, [1,3,2]);
cps1 = reshape(cps1, [prod(size(cps1))/dim, dim]);


% test scattering
% gbox is defined in g.box
leftB = gbox(1,1);
rightB = gbox(1,2);
bottomB = gbox(2,1);
topB = gbox(2,2);


x1 = round(cps1(:, 1)) - leftB + 1;
y1 = round(cps1(:, 2)) - bottomB + 1;

mask = sparse(y1, x1, 1, topB-bottomB+1, rightB-leftB+1);
mask = full(mask);
mask(mask>0) = 1;
