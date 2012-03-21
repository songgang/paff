function demo_velocity_field_step2_without_decreasing

% simulate a piece of mask by multiple beams of particle trajectory

g.non = 0;
g = get_control_point_and_transform_multibeam_two_rotations(g);
figure(12); clf;
disp_merging_velocity_field(g);
% print -dpng -f12 'data/trajectory_mask.png';




function disp_merging_velocity_field(g)

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

nb_cps = g.nb_cps;
cps = g.cps;

dim = g.dim;

% two collision boundaries
% r1 = intra-beam collision, r1 \approx 1
% r2 = inter-beam collision, r2 can be fairly big, like 10^2 or 20^2
% r is the collision boundary
r1 = 1;
g.r1 = 1;

r2 = 225;
g.r2 = r2;

% dv is the difference of the velocity when conflicted
g.dv = 1;

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
    
    
end;

clist = collision_detection_conflicted_velocity_multibeam(g, cpslist, vcpslist);
clist

xfield_0 = cat(3, X, Y);
yfield_current = xfield_0;

cpspiecelist = zeros(length(clist), 2, nb_cps);
% qslist = zeros(length(clist)-1, 2);

for ii = 1:nb_cps
    cpspiecelist(1, :, ii) = cps(:, ii)';
end;


% precompute the decreasing ratio
% compute the alpha decresing ratio from trajectory
g.h = 2; % h is the sampling radius of the points in the mask
g.sigma2 = Inf; % local decreasing sigma
g.sigma1 = 10; % the sigma to combine affine fields, need to be smaller than collision radius
g.boundary.s = 10;

% cpslist1 = permute(cpslist, [1,3,2]);
% cpslist1 = reshape(cpslist1, [prod(size(cpslist1))/g.dim, g.dim]);
% maskall = scatter_binary_trajectory_nearest_neighbor(g.boundary.box, cpslist1);

maskall = scatter_binary_trajectory_nearest_neighbor(g.boundary.box, cpslist);
logalpha = get_log_weight_using_distance_transform_tablegaussian(maskall, g.h, g.sigma2);
alpha = exp(logalpha);


for ii = 1: length(clist)-1
    
    ind_t1 = clist(ii);
    ind_t2 = clist(ii+1);
    
    
    % create velocity field using RBF decreasing function
    vfield = get_stationary_vield_copy_paste_decreasing(g, xfield_0, cpslist(ind_t1:ind_t2, :, :), alpha);
    % vfield = smooth_field(vfield, 1, 'Gaussian');
    
    mask_label = scatter_multiple_label_trajectory_after_distance_transform(g, cpslist(ind_t1:ind_t2, :, :));
    vfield_inside_mask_accurate = get_stationary_vfield_from_labeled_mask(g, xfield_0, mask_label);
    vfield = smooth_field(vfield_inside_mask_accurate, 20, 'PrecondtionVariationalWithBoundary', mask_label, vfield);
      
    
%     figure(15); clf;
%     imagesc(leftB:rightB, bottomB:topB, mask_label+1);colormap([1 1 1;0 1 0; 1 0 0]);
%     axis image;
    % print -dpng -f15 'data/trjactory_mask.png'
    
    
    figure(16); clf
    qq=20;
    q=10;
    image(leftB:rightB, bottomB:topB, mask_label+1); colormap([1 1 1;0 1 0; 1 0 0]);
    hold on;
    quiver(X(qq:q:end-qq+1,qq:q:end-qq+1), ...
        Y(qq:q:end-qq+1,qq:q:end-qq+1), ...
        vfield(qq:q:end-qq+1,qq:q:end-qq+1, 1) / max(vfield(:)) * 10, ...
        vfield(qq:q:end-qq+1,qq:q:end-qq+1, 2)  / max(vfield(:)) * 10, 1, 'LineWidth', 2);
    hold on;
    plot(g.aff{1}.polygon.x, g.aff{1}.polygon.y, 'g-', 'LineWidth', 2);
    plot(g.aff{2}.polygon.x, g.aff{2}.polygon.y, 'r-', 'LineWidth', 2);
    hold off;
    axis equal;
    axis([leftB, rightB, bottomB, topB]);
    print -dpng -f16 'data/merged_velocity_with_trajectory_without_decreasing.png'
    
    
    return;

    figure(18); clf
    qq=20;
    q=10;
    hold on;
    quiver(X(qq:q:end-qq+1,qq:q:end-qq+1), ...
        Y(qq:q:end-qq+1,qq:q:end-qq+1), ...
        vfield(qq:q:end-qq+1,qq:q:end-qq+1, 1) / max(vfield(:)) * 10, ...
        vfield(qq:q:end-qq+1,qq:q:end-qq+1, 2)  / max(vfield(:)) * 10, 1, 'LineWidth', 2);
    hold on;
    axis equal;
    axis([leftB, rightB, bottomB, topB]);
    print -dpng -f18 'data/merged_velocity.png'
    
    
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
fil=10;

figure(17); clf;
% % meshplot(X(pad*fil+1:fil:end-pad*fil, pad*fil+1:fil:end-pad*fil), Y(pad*fil+1:fil:end-pad*fil, pad*fil+1:fil:end-pad*fil), 'Color', 'g');
meshplot(yfield_current(pad*fil+1:fil:end-pad*fil, pad*fil+1:fil:end-pad*fil, 1), yfield_current(pad*fil+1:fil:end-pad*fil, pad*fil+1:fil:end-pad*fil, 2), 'Color', 'b', 'LineWidth', 2);
axis equal;
axis([leftB, rightB, bottomB, topB]);
print -dpng -f17 'data/merged_transform.png';



function g = get_control_point_and_transform_multibeam_two_rotations(g)


% control point
% column vector
% cps = [-10 10 0   0;
%         0  0  10 -10 ];
%      cps = [-10 0;
%               0 30];

% s is the sigma to computer affine velocity decreasing
% s can be removed in future implementation
s = 1;
dim = 2;

% simulation of two masks
% mask1 

[X, Y] = meshgrid(-80:10:-20, -20:10:20);
cps1 = [X(:)'; Y(:)']; % dim * nb_pts
ind1 = ones(1, size(cps1, 2));
cps = cps1;
ind = ind1;

[X, Y] = meshgrid(20:10:80, -40:10:40);
cps1 = [X(:)'; Y(:)']; % dim * nb_pts
ind1 = 2 * ones(1, size(cps1, 2));
cps = [cps, cps1];
ind = [ind, ind1];



g.cps = cps;
g.ind = ind;

nb_cps = size(cps, 2);


% define the affine transform for two masks
g.aff = cell(0);

[A, t] = get_A_and_t(40, [-50;0], [0; 0], [1, 1]);
[L, v] = get_Lv_from_At(A, t);
g.aff{end+1}.A = A;
g.aff{end}.t = t;
g.aff{end}.L = L;
g.aff{end}.v = v;
g.aff{end}.s = s;
g.aff{end}.polygon.x = [-80 -80 -20 -20 -80];
g.aff{end}.polygon.y = [20 -20 -20 20 20];


[A, t] = get_A_and_t(-40, [50; 0], [0; 0], [1,1]);
[L, v] = get_Lv_from_At(A, t);
g.aff{end+1}.A = A;
g.aff{end}.t = t;
g.aff{end}.L = L;
g.aff{end}.v = v;
g.aff{end}.s = s;
g.aff{end}.polygon.x = [20 20 80 80 20];
g.aff{end}.polygon.y = [40 -40 -40 40 40];


g.dim = dim;
g.nb_cps = nb_cps;
g.cps = cps;


