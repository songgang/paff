function vfield = get_stationary_vield_from_labeled_mask(g, xfield, mask_label)

X = xfield(:, :, 1);
Y = xfield(:, :, 2);
y = transpose([X(:), Y(:)]);

% rasterize the trajectory using the neearest neighbor for each mask
vfield_1 = zeros(size(X));
vfield_2 = zeros(size(X));

nb_mask = length(g.aff);

for ii = 1:nb_mask

    idx_mask = (mask_label == ii);
    
    affL = g.aff{ii}.L;
    affv = g.aff{ii}.v;
    
    viiy = affL * y(:, idx_mask) + affv*ones(1, nnz(idx_mask(:)));
    
    vfield_1(idx_mask) = viiy(1, :);
    vfield_2(idx_mask) = viiy(2, :);
end;


vfield = cat(3, vfield_1, vfield_2);