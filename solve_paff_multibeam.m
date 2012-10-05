% simulate a piece of mask by multiple beams of particle trajectory


figno = 11;

leftB = -150;
rightB = 150;
bottomB = -150;
topB = 150;


% define the boundary boxes
% [dim1a, dim1b; dim2a, dim2b];
g.boundary.box = [leftB, rightB; bottomB, topB];
g.boundary.s = 2;

[X, Y] = meshgrid(leftB:1:rightB, bottomB:1:topB);

% test on the Gaussian smoothing for velocity field
% this is used for decreasing from boundaries
g.vfield_smooth_sigma = 0.1;

% g.s2 is to combine weight
g.s2 = 1;

g = get_control_point_and_transform_multibeam_Lshape(g);
nb_cps = g.nb_cps;
cps = g.cps;

dim = g.dim;

% two collision boundaries
% r1 = intra-beam collision, r1 \approx 1
% r2 = inter-beam collision, r2 can be fairly big, like 10^2 or 20^2
% r is the collision boundary
r1 = 1;
g.r1 = 1;

r2 = 2250; %225;
g.r2 = r2;

% dv is the difference of the velocity when conflicted
g.dv = 1;


% precompute the decreasing ratio
% compute the alpha decresing ratio from trajectory
g.h = 3; % h is the sampling radius of the points in the mask
g.sigma2 = 25; % local decreasing sigma outside of all trajectory
g.sigma1 = 100; % the sigma to combine affine fields, need to be smaller than collision radius


% solve iniitial ode on p and q
% pqhat0 = [p; q];
cpshat0 = reshape(cps, [nb_cps*dim, 1]);
tlist = 0:0.01:1;


if 0
    % use ode to solve trajectory and define velocity
    [T, cpshat] = ode45(@(t, y) vhat_paff_multibeam(g, t, y), tlist, cpshat0);
    nb_T = length(T);
    cpslist = reshape(cpshat, [nb_T, dim, nb_cps]); % [pqhat(:, 3), pqhat(:, 4)];
    g.nb_T = nb_T;
    
    
    % recompute the velocity at each point
    % cpslist: nb_pts * nb_dim * nb_cps (idx_point, u/v, idx_cps)
    cpslist1 = permute(cpslist, [2, 1, 3]);
    cpslist2 = repmat(cpslist1, [1, nb_cps, 1]);
    vcpslist1 = v_paff_ex_pqvec_multibeam(g, nan, reshape(cpslist1, [dim, nb_T*nb_cps]), cpslist2);
    vcpslist = reshape(vcpslist1, [dim, nb_T, nb_cps]);
    vcpslist = permute(vcpslist, [2, 1, 3]);
    
    
else
    % use predefined trjactories, assuming they won't conflict
    
    tlist = 0:0.05:1;
    [cpslist, vcpslist] = predefine_trajectory_multibeam(g, tlist);
    g.nb_T = length(tlist);

    
    % create a cpslist of extended time duration, with t > 1, to "dilate" the
    % trajector mask in the time 
    tlist_ex = [tlist, 1.05:0.05:1.1];
    [cpslist_ex, vcpslist_ex] = predefine_trajectory_multibeam(g, tlist_ex);
    % create the long-time mask (from t=0:1, or even extended, t=-0.2:1.2)
    mask_longtime_label_ex = scatter_multiple_label_trajectory_after_distance_transform(g, cpslist_ex(:, :, :));
    mask_longtime_ex = (mask_longtime_label_ex > 0);
    
end;

clist = collision_detection_conflicted_velocity_multibeam(g, cpslist, vcpslist);
clist_ex = collision_detection_conflicted_velocity_multibeam(g, cpslist_ex, vcpslist_ex);
disp 'collision time points'
clist

disp 'collision extended time points'
clist_ex

xfield_0 = cat(3, X, Y);
yfield_current = xfield_0;

cpspiecelist = zeros(length(clist), 2, nb_cps);
% qslist = zeros(length(clist)-1, 2);

for ii = 1:nb_cps
    cpspiecelist(1, :, ii) = cps(:, ii)';
end;




% cpslist1 = permute(cpslist, [1,3,2]);
% cpslist1 = reshape(cpslist1, [prod(size(cpslist1))/g.dim, g.dim]);
% maskall = scatter_binary_trajectory_nearest_neighbor(g.boundary.box, cpslist1);


% maskall = scatter_binary_trajectory_nearest_neighbor(g.boundary.box, cpslist);
maskall = scatter_multiple_label_trajectory_after_distance_transform(g, cpslist(:, :, :));
maskall = maskall > 0;
% maskall = mask_longtime_ex;
logalpha = get_log_weight_using_distance_transform_tablegaussian(maskall, g.h, g.sigma2);
alpha = exp(logalpha);


for ii = 1: length(clist)-1
    
    ind_t1 = clist(ii);
    ind_t2 = clist(ii+1);    

    
    % create velocity field using RBF decreasing function
    vfield = get_stationary_vield_copy_paste_decreasing(g, xfield_0, cpslist(ind_t1:ind_t2, :, :), alpha);

    
    if 0
    
    % smoothing the velocity field in time by averaging the future velocity
    % field, aiming to removing folding patten. The smoothing region is the
    % non-critical region 

    % fixed value masks
    mask_label = scatter_multiple_label_trajectory_after_distance_transform(g, cpslist(ind_t1:ind_t2, :, :));
    el = strel('disk', g.h * 5);
    mask_label_dilate = imdilate(mask_label, el);
%     [gx, gy] = gradient(mask_label);
%     mask_label_bounds = (gx.*gx + gy.*gy > 0);
    mask_noncritical = (mask_label_dilate == 0);
%    mask_noncritical = mask_noncritical & (1- (mask_label > 0));
    
    if ii < length(clist)-1
        ind_t3 = clist(ii+2);
        vfield_forward = get_stationary_vield_copy_paste_decreasing(g, xfield_0, cpslist(ind_t2:ind_t3, :, :), alpha);        
    else % smooth backwards
        ind_t3 = clist_ex(min(find(clist_ex > clist(end) )));
        vfield_forward = get_stationary_vield_copy_paste_decreasing(g, xfield_0, cpslist_ex(ind_t2:ind_t3, :, :), alpha);
    end;
    
    
    
    q=10;
    figure(ii+2000); clf;
    quiver(X(1:q:end,1:q:end), ...
    Y(1:q:end,1:q:end), ...
        vfield_forward(1:q:end,1:q:end, 1) , ...
        vfield_forward(1:q:end,1:q:end, 2) , 1);
    hold on;
    for kk = 1:nb_cps;
    plot(cpslist_ex(ind_t2:ind_t3, 1, kk), cpslist_ex(ind_t2:ind_t3, 2, kk), 'r.');
    end;
    hold off;
    %
    axis image;

    
    
    
    vfield_1 = vfield(:, :, 1);
    vfield_2 = vfield(:, :, 2);
    vfield_forward_1 = vfield_forward(:, :, 1);
    vfield_forward_2 = vfield_forward(:, :, 2);

    tmpw = 1;
    vfield_1(mask_noncritical) = tmpw * vfield_1(mask_noncritical) + (1-tmpw) * vfield_forward_1(mask_noncritical);
    vfield_2(mask_noncritical) = tmpw * vfield_2(mask_noncritical) + (1-tmpw) * vfield_forward_2(mask_noncritical);
    vfield(:, :, 1) = vfield_1;
    vfield(:, :, 2) = vfield_2;

    end;
    
    
    % vfield = smooth_field(vfield, 1, 'Gaussian');

    % another way to create vector field:
    % fix the values inside the mask as the accurate values
    % and use smoothing / interpolation for the rest of the field
    
    
%     mask_label = scatter_multiple_label_trajectory_after_distance_transform(g, cpslist(ind_t1:ind_t2, :, :));
%     vfield_inside_mask_accurate = get_stationary_vfield_from_labeled_mask(g, xfield_0, mask_label);
%     vfield = smooth_field(vfield_inside_mask_accurate, 20, 'PrecondtionVariationalWithBoundary', mask_label, vfield);
      
    
    % vfield = get_stationary_vield_copy_paste(g, xfield_0, cpslist(ind_t1:ind_t2, :, :), tlist(ind_t1:ind_t2));
    
    % [vfield_backbone, ind_backbone] = scatter_interploate_trajectory_nearest_neighbor([leftB, rightB, bottomB, topB], cpslist(ind_t1:ind_t2, :, :), vcpslist(ind_t1:ind_t2, :, :));
    % vfield = smooth_field(vfield_backbone, g.vfield_smooth_sigma, 'PrecondtionVariational');
    
    
    
    
    % replace v on trajectories with exame numbers
    % vfield(ind_backbone) = vfield_backbone(ind_backbone);
    % vfield = vfield_backbone;
    
    
    q=10;
    figure(ii+1000); clf;
%     quiver(X(1:q:end,1:q:end), ...
%         Y(1:q:end,1:q:end), ...
%         vfield(1:q:end,1:q:end, 1) / max(vfield(:)) * 2, ...
%         vfield(1:q:end,1:q:end, 2)  / max(vfield(:)) * 2, 0);

    quiver(X(1:q:end,1:q:end), ...
        Y(1:q:end,1:q:end), ...
        vfield(1:q:end,1:q:end, 1) , ...
        vfield(1:q:end,1:q:end, 2) , 1);

    hold on;
    for kk = 1:nb_cps;
        plot(cpslist(ind_t1:ind_t2, 1, kk), cpslist(ind_t1:ind_t2, 2, kk), 'r.');
        % quiver(cpslist(ind_t1:ind_t2, 1, kk), cpslist(ind_t1:ind_t2, 2, kk), vcpslist(ind_t1:ind_t2, 1 ,kk) * 0.02, vcpslist(ind_t1:ind_t2, 2 ,kk) * 0.02, 1, 'r');
    end;
    hold off;
    %
    axis image;
    
    yfield_delta = exp_mapping(vfield, X, Y, tlist(ind_t2)-tlist(ind_t1), 10);
    yfield_current = compose_phi(yfield_current, yfield_delta, X, Y);
    
    
    for k = 1:nb_cps
        idx_p = find(X==cps(1,k) & Y==cps(2,k));
        y1 = yfield_current(:,:,1);
        y2 = yfield_current(:,:,2);
        cpspiecelist(ii+1, :, k) = [y1(idx_p), y2(idx_p)];
    end;
    
end;


pad=0;
fil=5;

figure(109); clf;
% plot trajectory of control points (cps) using ode solution
hold on;
clrs='gbr';
markers='o*v';
for jj = 1:length(clist)
    for ii = 1:nb_cps
    %     plot(squeeze(cpslist(:, 1, ii)), squeeze(cpslist(:, 2, ii)), ['-', clrs(mod(ii, length(clrs))+1), '*']);
        plot(squeeze(cpslist(clist(jj), 1, ii)), squeeze(cpslist(clist(jj), 2, ii)), ...
            [clrs(mod(ii, length(clrs))+1), markers(mod(jj, length(markers))+1)],  'MarkerSize', 10);
    end;
end;
hold off;

% meshplot(X(pad*fil:fil:end-pad*fil, pad*fil:fil:end-pad*fil), Y(pad*fil:fil:end-pad*fil, pad*fil:fil:end-pad*fil), 'Color', 'g');
% meshplot(yfield_current(pad*fil:fil:end-pad*fil, pad*fil:fil:end-pad*fil, 1), yfield_current(pad*fil:fil:end-pad*fil, pad*fil:fil:end-pad*fil, 2), 'Color', 'b');

meshplot(X(pad*fil+1:fil:end-pad*fil, pad*fil+1:fil:end-pad*fil), Y(pad*fil+1:fil:end-pad*fil, pad*fil+1:fil:end-pad*fil), 'Color', 'g');
meshplot(yfield_current(pad*fil+1:fil:end-pad*fil, pad*fil+1:fil:end-pad*fil, 1), yfield_current(pad*fil+1:fil:end-pad*fil, pad*fil+1:fil:end-pad*fil, 2), 'Color', 'b');

% hold on;
% for jj = 1:nb_cps
%     for ii = size(cpspiecelist, 1):size(cpspiecelist, 1)
%         quiver(cpspiecelist(1,1,jj), ...
%             cpspiecelist(1,2,jj), ...
%             cpspiecelist(ii,1,jj)-cpspiecelist(1,1,jj), ...
%             cpspiecelist(ii,2,jj)-cpspiecelist(1,2,jj), ...
%             1,  clrs(mod(jj, length(clrs))+1), 'LineWidth', 1);
%     end;
% end;
% hold off;




% plot desired position
% cps_desired_target = get_affine_on_cps(g);
% hold on;
% for jj = 1:nb_cps
%         quiver(cps(1,jj), ...
%                g.cps(2,jj), ...
%                cps_desired_target(1,jj)-cps(1,jj), ...
%                cps_desired_target(2,jj)-cps(2,jj), ...
%                0, 'c', 'LineWidth', 2);
% end;
% hold off;




axis equal;

%
% figure(3); clf;
% plot(clist, tlist(clist), 'b*');
% title('clist');

