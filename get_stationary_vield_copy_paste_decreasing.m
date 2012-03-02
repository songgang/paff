% create velocity field using RBF decreasing function
function vfield = get_stationary_vield_copy_paste_decreasing(g, xfield, cpslist, alpha)


nb_yaxis = size(xfield, 1);
nb_xaxis = size(xfield, 2);
X = xfield(:, :, 1);
Y = xfield(:, :, 2);
y = transpose([X(:), Y(:)]);


% rasterize the trajectory using the neearest neighbor for each mask


    
% compute the log of weights w
nb_y = prod(size(X));
nb_mask = length(g.aff);

logwlist = zeros(nb_mask, nb_y);
for ii = 1:nb_mask;
    idx = find(g.ind == ii);
    maskii = scatter_binary_trajectory_nearest_neighbor(g.boundary.box, cpslist(:, :, idx));
    % a1 = get_log_weight_using_distance_transform_tablegaussian(maskii, g.h, g.sigma1);
    a1 = get_log_weight_using_distance_transform_tablelaplacian(maskii, g.h, g.sigma1);
    
    logwlist(ii, :) = a1(:);
end;

v = zeros(g.dim, nb_y);
for ii = 1:nb_mask
    sw = zeros(1, nb_y);
    for jj = 1:nb_mask
        sw = sw + exp(logwlist(jj, :) - logwlist(ii, :));
    end;
    wii = 1 ./ sw;
    affL = g.aff{ii}.L;
    affv = g.aff{ii}.v;
    viiy = affL * y + affv*ones(1, nb_y);
    v = v + (ones(g.dim, 1)*wii).*viiy;
    
    % figure; imagesc(-150:150, -150:150, reshape(wii, [301, 301])); axis image; axis xy;
end;


% decreasing velocity from trajectory
for ii = 1:g.dim
    v(ii, :) = v(ii, :) .* alpha(:)';
end;


% fix boundaries

if 1 % boundary condition
    logwboundlist = zeros(g.dim*2, nb_y);
    for ii = 1:g.dim
        logwboundlist(ii*2-1, :) = loggpaff(y(ii, :) - g.boundary.box(ii,1), g.boundary.s);
        logwboundlist(ii*2, :) = loggpaff(y(ii, :) - g.boundary.box(ii,2), g.boundary.s);
    end;

    % velocity is scaled by distance to the fixed boundary
    swratio = -1 * inf(1, nb_y);
    for jj = 1:g.dim*2
        swratio = max(swratio, logwboundlist(jj, :));
    end;
    swratio = 1 - exp(swratio);
    v = v.* (ones(g.dim, 1) * swratio);
end;


if isnan(v)
    disp 'haha';
end;






vfield(:, :, 1) = reshape(v(1, :), [nb_yaxis, nb_xaxis]);
vfield(:, :, 2) = reshape(v(2, :), [nb_yaxis, nb_xaxis]);




