function [vfield, ind_backbone] = scatter_interploate_trajectory_nearest_neighbor(bbox, cpslist, vcpslist)

% test scattering
leftB = bbox(1);
rightB = bbox(2);
bottomB = bbox(3);
topB = bbox(4);

[X, Y] = meshgrid(leftB:1:rightB, bottomB:1:topB);
vfield = zeros([size(X), 2]);

ind_backbone = [];

nb_cps = size(cpslist, 3);

for ii = 1:nb_cps
    x1 = round(cpslist(:, 1, ii)) - leftB + 1;
    y1 = round(cpslist(:, 2, ii)) - bottomB + 1;
    
    ind = sub2ind(size(vfield), y1, x1, ones(size(y1)));
    vfield(ind) = vcpslist(:, 1, ii);
    ind_backbone = [ind_backbone; ind];

    ind = sub2ind(size(vfield), y1, x1, 2 * ones(size(y1)));
    vfield(ind) = vcpslist(:, 2, ii);
    ind_backbone = [ind_backbone; ind];
end;
