function mask_label = scatter_multiple_label_trajectory_after_distance_transform(g, cpslist)

% reshape cps to a 2D array: nb_points * nb_dim
dim = size(cpslist, 2);
nb_mask = length(g.aff);
nb_cps = size(cpslist, 3);

% test scattering
% gbox is defined in g.box
gbox = g.boundary.box;
leftB = gbox(1,1);
rightB = gbox(1,2);
bottomB = gbox(2,1);
topB = gbox(2,2);


mask_label = zeros(topB-bottomB+1, rightB-leftB+1);

for ii = 1:nb_mask
    
    idx_cps = (g.ind == ii);
    
    % reshape cps to a 2D array: nb_points * nb_dim
    cps0 = cpslist(:, :, idx_cps);
    cps1 = permute(cps0, [1,3,2]);
    cps1 = reshape(cps1, [prod(size(cps1))/dim, dim]);
    
    x1 = round(cps1(:, 1)) - leftB + 1;
    y1 = round(cps1(:, 2)) - bottomB + 1;

    mask = sparse(y1, x1, 1, topB-bottomB+1, rightB-leftB+1);
    mask = full(mask);
    mask(mask>0) = 1;
    
    maskfill = (bwdist(mask) <= g.h);
    
    mask_label(maskfill > 0) = ii;

end;

    


