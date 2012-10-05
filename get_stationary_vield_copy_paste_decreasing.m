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


if 0
% plain way to compute weight using expoential of the distance to each mask
for ii = 1:nb_mask;
    idx = find(g.ind == ii);
    % maskii = scatter_binary_trajectory_nearest_neighbor(g.boundary.box, cpslist(:, :, idx));
    maskii = scatter_single_label_trajectory_after_distance_transform(g, cpslist(:, :, idx));
    % a1 = get_log_weight_using_distance_transform_tablegaussian(maskii, g.h, g.sigma1);
    a1 = get_log_weight_using_distance_transform_tablelaplacian(maskii, 0, g.sigma1);
    % a1 = get_log_weight_using_distance_transform_tablelaplacian(maskii, 0, g.sigma1);
    
    logwlist(ii, :) = a1(:);
end;

else

% test a complicated way of determine weight should be determined by which
% mask 

dist_per_mask = zeros(nb_yaxis, nb_xaxis, nb_mask);
for ii = 1:nb_mask
    idx = find(g.ind == ii);
    maskii = scatter_single_label_trajectory_after_distance_transform(g, cpslist(:, :, idx));
    d = bwdist(maskii);
    d = max(0, d-g.h);
    dist_per_mask(:, :, ii) = d;
end;
% divide influence region by finding the minimum distance
[void, influence_region] = min(dist_per_mask, [], 3);    
% for each influence region, computing SIGNED distance transform, such
% distance is used to compute the weight: w
% insdie mask and far from boundary: w = 1
% insdie mask and close to boundary: W = 1 to 0.5
% outside mask and close to boundary: w = 0.5 to 0
% outside mask and far from boundary: w = 0
for ii = 1:nb_mask
    current_region_mask = (influence_region == ii); 
    dout = bwdist(current_region_mask);
    din = bwdist(1 - current_region_mask);
    
    a1(din > sqrt(g.sigma1)) = 1;
    a1(din <= sqrt(g.sigma1) & din > 0) =  din(din <= sqrt(g.sigma1)  & din > 0) / sqrt(g.sigma1) * 0.5 + 0.5;
    a1(dout <= sqrt(g.sigma1) & dout > 0) =  0.5 - dout(dout <= sqrt(g.sigma1)  & dout > 0) / sqrt(g.sigma1) * 0.5;
    a1(dout > sqrt(g.sigma1)) =  0;
    
    % make sure inside each trajectory the weight is always 1
    a1(dist_per_mask(:, :, ii) == 0) = 1;
    
    logwlist(ii, :) = log(a1(:));
end;

end;



v = zeros(g.dim, nb_y);
for ii = 1:nb_mask
    sw = zeros(1, nb_y);
    for jj = 1:nb_mask
        sw = sw + exp(logwlist(jj, :) - logwlist(ii, :));
        
    end;
    
    for jj = 1:nb_mask
        % when logwlist(jj, :) == 0 means the points are inside other masks, 
        % sw has to be equal to 0 = 0/(0+1) for this case
        if (ii~=jj)
            sw(logwlist(jj, :) == 0) = Inf;
        else
            sw(logwlist(jj, :) == 0) = 1;
        end;
    end;
    
    wii = 1 ./ sw;
    
    
    affL = g.aff{ii}.L;
    affv = g.aff{ii}.v;
    viiy = affL * y + affv*ones(1, nb_y);
    v = v + (ones(g.dim, 1)*wii).*viiy;
    
%     figure; imagesc(-150:150, -150:150, reshape(wii, [301, 301])); axis image; axis xy;
end;


% decreasing velocity from trajectory
for ii = 1:g.dim
    v(ii, :) = v(ii, :) .* alpha(:)';
end;


% fix boundaries

if g.boundary.s > 0 % boundary condition
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




