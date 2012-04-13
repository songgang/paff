function demo_velocity_field_step2

% simulate a piece of mask by multiple beams of particle trajectory

g.non = 0;
% print -dpng -f12 'data/trajectory_mask.png';

g = get_control_point_and_transform_multibeam_comp_arsigny(g);

mycase = 'arsigny';
disp_merging_velocity_field(g, mycase);

mycase = 'ours';
disp_merging_velocity_field(g, mycase);


function disp_merging_velocity_field(g, mycase)

leftB = -150;
rightB = 150;
bottomB = -150;
topB = 150;

switch mycase
    case 'arsigny',
        fignobase = 100;
    case 'ours',
        fignobase = 200;
end;


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
r1 = 0.25;
g.r1 = r1;

r2 = 100;
g.r2 = r2;

% dv is the difference of the velocity when conflicted
g.dv = 1;

% solve iniitial ode on p and q
% pqhat0 = [p; q];
cpshat0 = reshape(cps, [nb_cps*dim, 1]);
tlist = 0:0.01:1;


% use predefined trjactories, assuming they won't conflict

tlist = 0:0.02:1;
[cpslist_groundtruth, vcpslist] = predefine_trajectory_multibeam(g, tlist);


if strcmp(mycase, 'arsigny')
    % to fake Arsigny's results
    tlist = [0, 0];
    [cpslist, vcpslist] = predefine_trajectory_multibeam(g, tlist);
    g.nb_T = length(tlist);
    tlist = [0, 1];
elseif strcmp(mycase, 'ours')
    [cpslist, vcpslist] = predefine_trajectory_multibeam(g, tlist);
    g.nb_T = length(tlist);
end;
    

clist = collision_detection_conflicted_velocity_multibeam(g, cpslist, vcpslist);
%clist = [1, 2];

xfield_0 = cat(3, X, Y);
yfield_current = xfield_0;

cpspiecelist = zeros(length(clist), 2, nb_cps);
% qslist = zeros(length(clist)-1, 2);

for ii = 1:nb_cps
    cpspiecelist(1, :, ii) = cps(:, ii)';
end;


% precompute the decreasing ratio
% compute the alpha decresing ratio from trajectory
g.h = 20; % h is the sampling radius of the points in the mask
g.sigma2 = 100000; % local decreasing sigma
g.sigma1 = 100; % the sigma to combine affine fields, need to be smaller than collision radius
g.boundary.s = 2;

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
    % vfield = smooth_field(vfield_inside_mask_accurate, 20, 'PrecondtionVariationalWithBoundary', mask_label, vfield);
     vfield = smooth_field(vfield, 5, 'Gaussian');
      
    
    figure(fignobase+15); clf;
    qq=20;
    q=20;
    imagesc(leftB:rightB, bottomB:topB, mask_label+1);colormap([1 1 1;0 1 0; 1 0 0]);
    hold on;
    vfield_seed = vfield .* repmat(mask_label > 0, [1 1 2]);
    quiver(X(qq:q:end-qq+1,qq:q:end-qq+1), ...
        Y(qq:q:end-qq+1,qq:q:end-qq+1), ...
        vfield_seed(qq:q:end-qq+1,qq:q:end-qq+1, 1) / max(vfield(:)) * 10, ...
        vfield_seed(qq:q:end-qq+1,qq:q:end-qq+1, 2)  / max(vfield(:)) * 10, 1, 'LineWidth', 2);     
    hold off;
    axis xy;
    axis equal;
    axis([leftB, rightB, bottomB, topB]);
    switch mycase
    case 'arsigny',
        print -dpng -f115 'data/arsigny_velocity_seed.png'
    case 'ours',
        print -dpng -f215 'data/ours_velocity_seed.png'
    end;
    
    
    
    figure(fignobase+16); clf
    qq=20;
    q=20;
    image(leftB:rightB, bottomB:topB, mask_label+1); colormap([1 1 1;0 1 0; 1 0 0]);
    hold on;
    quiver(X(qq:q:end-qq+1,qq:q:end-qq+1), ...
        Y(qq:q:end-qq+1,qq:q:end-qq+1), ...
        vfield(qq:q:end-qq+1,qq:q:end-qq+1, 1) / max(vfield(:)) * 10, ...
        vfield(qq:q:end-qq+1,qq:q:end-qq+1, 2)  / max(vfield(:)) * 10, 1, 'LineWidth', 2);
    hold on;
    hold off;

    axis xy;
    axis equal;
    axis([leftB, rightB, bottomB, topB]);
    switch mycase
    case 'arsigny',
        print -dpng -f116 'data/arsigny_velocity_final.png'
    case 'ours',
        print -dpng -f216 'data/ours_velocity_final.png'
    end;

    
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

figure(fignobase+17); clf;
hold on;
q=1;
meshplot(yfield_current(pad*fil+1+q:fil:end-pad*fil-q, pad*fil+1+q:fil:end-pad*fil-q, 1), ...
    yfield_current(pad*fil+1+q:fil:end-pad*fil-q, pad*fil+1+q:fil:end-pad*fil-q, 2), 'Color', 'b', 'LineWidth', 2);
% draw the ground truth
hold on;
clrlist= 'gr';
for kk = 1:length(g.ind_disp)
    ii = g.ind_disp(kk);
    label = g.ind(ii);
    plot(squeeze(cpslist_groundtruth(:, 1, ii)), squeeze(cpslist_groundtruth(:, 2, ii)), [clrlist(label), '--'], 'LineWidth', 2);
    
end;
for kk = 1:2;
    ii = g.ind_boundary(find(g.ind(g.ind_boundary) == kk));
    label = kk;
    plot(g.cps(1, ii), g.cps(2, ii), [clrlist(label), '-'], 'LineWidth', 2);
    plot(squeeze(cpslist_groundtruth(end, 1, ii)), squeeze(cpslist_groundtruth(end, 2, ii)), [clrlist(label), '--'], 'LineWidth', 2);
    
end;
hold off;

% draw the final destination of those polygons
hold on;
clrlist= 'gr';
for kk = 1:length(g.ind_disp)
    ii = g.ind_disp(kk);
    label = g.ind(ii);
    phix = interp2(X, Y, yfield_current(:, :, 1), g.cps(1,ii), g.cps(2,ii));
    phiy = interp2(X, Y, yfield_current(:, :, 2), g.cps(1,ii), g.cps(2,ii));
    plot([g.cps(1,ii), phix], [g.cps(2,ii), phiy], [clrlist(label), '-'], 'LineWidth', 2);
end;

for kk = 1:2;
    ii = g.ind_boundary(find(g.ind(g.ind_boundary) == kk));
    label = kk;
    
    phix = interp2(X, Y, yfield_current(:, :, 1), g.cps(1,ii), g.cps(2,ii));
    phiy = interp2(X, Y, yfield_current(:, :, 2), g.cps(1,ii), g.cps(2,ii));
    plot([g.cps(1,ii), phix], [g.cps(2,ii), phiy], [clrlist(label), '*'], 'LineWidth', 2);
    
    plot(g.cps(1, ii), g.cps(2, ii), [clrlist(label), '-.'], 'LineWidth', 2);
    plot(phix, phiy, [clrlist(label), '-'], 'LineWidth', 2);
    
end;
hold off;
axis equal;
axis([leftB, rightB, bottomB, topB]);
switch mycase
case 'arsigny',
    print -dpng -f117 'data/arsigny_deformation.png'
case 'ours',
    print -dpng -f217 'data/ours_deformation.png'
end;








figure(fignobase+18); clf;
hold on;
% draw the ground truth
hold on;
clrlist= 'gr';
for kk = 1:length(g.ind_disp)
    ii = g.ind_disp(kk);
    label = g.ind(ii);
    plot(squeeze(cpslist_groundtruth(:, 1, ii)), squeeze(cpslist_groundtruth(:, 2, ii)), [clrlist(label), '--'], 'LineWidth', 2);
    
end;
for kk = 1:2;
    ii = g.ind_boundary(find(g.ind(g.ind_boundary) == kk));
    label = kk;
    plot(g.cps(1, ii), g.cps(2, ii), [clrlist(label), '-'], 'LineWidth', 2);
    plot(squeeze(cpslist_groundtruth(end, 1, ii)), squeeze(cpslist_groundtruth(end, 2, ii)), [clrlist(label), '--'], 'LineWidth', 2);
    
end;
hold off;

axis equal;
axis([leftB, rightB, bottomB, topB]);
switch mycase
case 'arsigny',
    print -dpng -f118 'data/arsigny_two_affine_input.png'
case 'ours',
    print -dpng -f218 'data/ours_two_affine_input.png'
end;








return;



function g = get_control_point_and_transform_multibeam_comp_arsigny(g)


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

[X, Y] = meshgrid(-80:20:-40, 60:20:100);
cps1 = [X(:)'; Y(:)']; % dim * nb_pts
nb_mask1 = size(cps1, 2);
ind1 = ones(1, size(cps1, 2));
cps = cps1;
ind = ind1;
ind_boundary = [find(X==-80 & Y==60), find(X==-80 & Y==100), find(X==-40 & Y==100), find(X==-40 & Y==60), find(X==-80 & Y==60)];
g.ind_boundary = ind_boundary;
ind_disp = round(nb_mask1 / 2);
g.ind_disp = ind_disp;

[X, Y] = meshgrid(-80:20:-40, -100:20:-60);
cps1 = [X(:)'; Y(:)']; % dim * nb_pts
ind1 = 2 * ones(1, size(cps1, 2));
cps = [cps, cps1];
ind = [ind, ind1];

g.cps = cps;
g.ind = ind;
ind_boundary = nb_mask1+ [find(X==-80 & Y==-100), find(X==-80 & Y==-60), find(X==-40 & Y==-60), find(X==-40 & Y==-100), find(X==-80 & Y==-100)];
g.ind_boundary =[g.ind_boundary, ind_boundary];
% =
ind_disp = round(size(cps1, 2) / 2) + nb_mask1;
g.ind_disp = [g.ind_disp, ind_disp];


nb_cps = size(cps, 2);


% define the affine transform for two masks
g.aff = cell(0);

[A, t] = get_A_and_t(0, [-60;60], [160; 0], [1, 1]);
[L, v] = get_Lv_from_At(A, t);
g.aff{end+1}.A = A;
g.aff{end}.t = t;
g.aff{end}.L = L;
g.aff{end}.v = v;
g.aff{end}.s = s;
g.aff{end}.polygon.x = [-80 -80 -20 -20 -80];
g.aff{end}.polygon.y = [20 -20 -20 20 20];

[A, t] = get_A_and_t(0, [0; 60], [120; 70], [1,1]);
% [A, t] = get_A_and_t(60, [0; 60], [0; 0], [1,1]);
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


