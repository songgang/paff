function demo_transform_objects

% warping the lung lobes
% translation

I = double(imread('data/01_Fixed_y161ud.png'));
Ilabel = double(imread('data/01_Fixed_y161ud_label.png'));
Ifinal = double(imread('data/01_Fixed_y161ud_final_warped.png'));

figure(16); clf;
imagesc(Ilabel); axis image;

cnt = 0;
for manual_ratio = 0:0.1:1
    cnt = disp_transform_lung(I, Ilabel, Ifinal, manual_ratio, 1, cnt);
end;

for manual_ratio = 0.1:0.1:1
    cnt = disp_transform_lung(I, Ilabel, Ifinal, manual_ratio, 2, cnt);
end;


for manual_ratio = 0.25:0.25:1
    cnt = disp_transform_lung(I, Ilabel, Ifinal, manual_ratio, 3, cnt);
end;


function cnt = disp_transform_lung(I, Ilabel, Ifinal, manual_ratio, step_id, cnt)


% print -dpng -f16 'data/01_Fixed_y161ud_lable_colored.png';


leftB = 1;
rightB = size(I, 2);
bottomB = 1;
topB = size(I, 1);


% define the boundary boxes
% [dim1a, dim1b; dim2a, dim2b];
g.boundary.box = [leftB, rightB; bottomB, topB];
g.boundary.s = 10;
g.Ilabel = Ilabel;

[X, Y] = meshgrid(leftB:1:rightB, bottomB:1:topB);

% test on the Gaussian smoothing for velocity field
% this is used for decreasing from boundaries
g.vfield_smooth_sigma = 0.1;

% g.s2 is to combine weight
% g.s2 = 1;

switch step_id
    case 1,
        g = get_control_point_and_transform_multibeam_lobes(g, manual_ratio);
    case 2,
        g = get_control_point_and_transform_multibeam_lobes_step2(g, manual_ratio);
    case 3,
        g = get_control_point_and_transform_multibeam_lobes_step3(g, manual_ratio);

end;



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

tlist = 0:0.05:1;
[cpslist, vcpslist] = predefine_trajectory_multibeam(g, tlist);
g.nb_T = length(tlist);

clist = collision_detection_conflicted_velocity_multibeam(g, cpslist, vcpslist);
clist

xfield_0 = cat(3, X, Y);
yfield_current = xfield_0;

cpspiecelist = zeros(length(clist), 2, nb_cps);

for ii = 1:nb_cps
    cpspiecelist(1, :, ii) = cps(:, ii)';
end;



% precompute the decreasing ratio
% compute the alpha decresing ratio from trajectory
g.h = 20; % h is the sampling radius of the points in the mask
g.sigma2 = 20; % local decreasing sigma
g.sigma1 = 1; % the sigma to combine affine fields, need to be smaller than collision radius

% descreasing weight map from trajectory
maskall = scatter_binary_trajectory_nearest_neighbor(g.boundary.box, cpslist);
logalpha = get_log_weight_using_distance_transform_tablegaussian(maskall, g.h, g.sigma2);
alpha = exp(logalpha);



figure(20); clf;

subplot(2,2,3);


% figure(13); clf;
imagesc(I); axis image; colormap gray;

for ii = 1: length(clist)-1    
    ind_t1 = clist(ii);
    ind_t2 = clist(ii+1);
    % create velocity field using RBF decreasing function
    vfield = get_stationary_vield_copy_paste_decreasing(g, xfield_0, cpslist(ind_t1:ind_t2, :, :), alpha);
    
    q=2;
    
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
fil=15;

% plot trajectory of control points (cps) using ode solution
hold on;
clrs='gbr';
for ii = 1:nb_cps
    %     plot(squeeze(cpslist(:, 1, ii)), squeeze(cpslist(:, 2, ii)), ['-', clrs(mod(ii, length(clrs))+1), '*']);
    %    plot(squeeze(cpslist([1, end], 1, ii)), squeeze(cpslist([1, end], 2, ii)), [clrs(mod(ii, length(clrs))+1), '*'],  'MarkerSize', 10);
end;
hold off;

% meshplot(X(pad*fil:fil:end-pad*fil, pad*fil:fil:end-pad*fil), Y(pad*fil:fil:end-pad*fil, pad*fil:fil:end-pad*fil), 'Color', 'g');
% meshplot(yfield_current(pad*fil:fil:end-pad*fil, pad*fil:fil:end-pad*fil, 1), yfield_current(pad*fil:fil:end-pad*fil, pad*fil:fil:end-pad*fil, 2), 'Color', 'b');

% meshplot(X(pad*fil+1:fil:end-pad*fil, pad*fil+1:fil:end-pad*fil), Y(pad*fil+1:fil:end-pad*fil, pad*fil+1:fil:end-pad*fil), 'Color', 'g');
meshplot(yfield_current(pad*fil+1:fil:end-pad*fil, pad*fil+1:fil:end-pad*fil, 1), yfield_current(pad*fil+1:fil:end-pad*fil, pad*fil+1:fil:end-pad*fil, 2), 'Color', 'g', 'LineWidth', 1);



axis([leftB, rightB, bottomB, topB]);
axis image;
axis ij;
axis off;

title('Deformation Field on Moving Image');


Ilabelwarped = griddata(yfield_current(:, :, 1), yfield_current(:, :, 2), Ilabel, X, Y, 'nearest');

subplot(2,2,1);

% figure(12); clf;
IlabelRGB = uint8(ind2rgb(Ilabelwarped+1, [0 0 0; 1 0 0; 0 1 0]) * 255);
h =imshow(IlabelRGB);
hold on;
h = imshow(Ifinal); 
hold off; 
set(h, 'AlphaData', 255-Ifinal);
title('Interactive Mask on Fix Image');

subplot(2,2,2);

% figure(14); clf;
imagesc(Ifinal); axis image; colormap gray;
title('Desired Fixed Image');
axis off;

Iwarped = griddata(yfield_current(:, :, 1), yfield_current(:, :, 2), I, X, Y);

subplot(2,2,4);

%figure(15); clf;
imagesc(Iwarped); axis image; colormap gray;
title('Interactive Warped Moving Image');
axis off;
% imwrite(uint8(Iwarped), 'data/01_Fixed_y161ud_final_warped.png');

cnt = cnt + 1;
imgname = sprintf('data/interactive_%05d.png', cnt); 
print('-dpng', '-f20', imgname);

function g = get_control_point_and_transform_multibeam_lobes(g, manual_ratio)
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


Ilabel = double(imread('data/01_Fixed_y161ud_label.png'));

% find control points inside each label
imgw = size(Ilabel, 2);
imgh = size(Ilabel, 1);

[X, Y] = meshgrid(1:20:imgw, 1:20:imgh);
ind_all = sub2ind(size(Ilabel), Y(:), X(:));

X1 = X(Ilabel(ind_all)==1);
Y1 = Y(Ilabel(ind_all)==1);
cps1 = [X1(:)'; Y1(:)'];
ind1 = ones(1, size(cps1, 2));
cps = cps1;
ind = ind1;

X1 = X(Ilabel(ind_all)==2);
Y1 = Y(Ilabel(ind_all)==2);
cps1 = [X1(:)'; Y1(:)']; % dim * nb_pts
ind1 = 2 * ones(1, size(cps1, 2));
cps = [cps, cps1];
ind = [ind, ind1];


g.cps = cps;
g.ind = ind;

nb_cps = size(cps, 2);

% define the affine transform for two masks
g.aff = cell(0);


% [A, t] = get_A_and_t(-15, [94; 141], [10; -10], [0.8,0.8])
[A, t] = get_A_and_t(0, [94; 141], [0; 0], [1, 1])
[L, v] = get_Lv_from_At(A, t);
g.aff{end+1}.A = A;
g.aff{end}.t = t;
g.aff{end}.L = L;
g.aff{end}.v = v;
g.aff{end}.s = s;


[A, t] = get_A_and_t(0, [332;142], [-60 * manual_ratio; -20], [1, 1])
[L, v] = get_Lv_from_At(A, t);
g.aff{end+1}.A = A;
g.aff{end}.t = t;
g.aff{end}.L = L;
g.aff{end}.v = v;
g.aff{end}.s = s;


g.dim = dim;
g.nb_cps = nb_cps;
g.cps = cps;
g.Ilabel = Ilabel;





function g = get_control_point_and_transform_multibeam_lobes_step2(g, manual_ratio)
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


Ilabel = double(imread('data/01_Fixed_y161ud_label.png'));

% find control points inside each label
imgw = size(Ilabel, 2);
imgh = size(Ilabel, 1);

[X, Y] = meshgrid(1:20:imgw, 1:20:imgh);
ind_all = sub2ind(size(Ilabel), Y(:), X(:));

X1 = X(Ilabel(ind_all)==1);
Y1 = Y(Ilabel(ind_all)==1);
cps1 = [X1(:)'; Y1(:)'];
ind1 = ones(1, size(cps1, 2));
cps = cps1;
ind = ind1;

X1 = X(Ilabel(ind_all)==2);
Y1 = Y(Ilabel(ind_all)==2);
cps1 = [X1(:)'; Y1(:)']; % dim * nb_pts
ind1 = 2 * ones(1, size(cps1, 2));
cps = [cps, cps1];
ind = [ind, ind1];


g.cps = cps;
g.ind = ind;

nb_cps = size(cps, 2);

% define the affine transform for two masks
g.aff = cell(0);


% [A, t] = get_A_and_t(-15, [94; 141], [10; -10], [0.8,0.8])
[A, t] = get_A_and_t(-30 * manual_ratio, [94; 141], [0; 0], [1, 1])
[L, v] = get_Lv_from_At(A, t);
g.aff{end+1}.A = A;
g.aff{end}.t = t;
g.aff{end}.L = L;
g.aff{end}.v = v;
g.aff{end}.s = s;


[A, t] = get_A_and_t(0, [332;142], [-60 * 1; -20], [1, 1])
[L, v] = get_Lv_from_At(A, t);
g.aff{end+1}.A = A;
g.aff{end}.t = t;
g.aff{end}.L = L;
g.aff{end}.v = v;
g.aff{end}.s = s;


g.dim = dim;
g.nb_cps = nb_cps;
g.cps = cps;
g.Ilabel = Ilabel;



function g = get_control_point_and_transform_multibeam_lobes_step3(g, manual_ratio)
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


Ilabel = g.Ilabel; % double(imread('data/01_Fixed_y161ud_label.png'));

% find control points inside each label
imgw = size(Ilabel, 2);
imgh = size(Ilabel, 1);

[X, Y] = meshgrid(1:20:imgw, 1:20:imgh);
ind_all = sub2ind(size(Ilabel), Y(:), X(:));

X1 = X(Ilabel(ind_all)==1);
Y1 = Y(Ilabel(ind_all)==1);
cps1 = [X1(:)'; Y1(:)'];
ind1 = ones(1, size(cps1, 2));
cps = cps1;
ind = ind1;

X1 = X(Ilabel(ind_all)==2);
Y1 = Y(Ilabel(ind_all)==2);
cps1 = [X1(:)'; Y1(:)']; % dim * nb_pts
ind1 = 2 * ones(1, size(cps1, 2));
cps = [cps, cps1];
ind = [ind, ind1];


g.cps = cps;
g.ind = ind;

nb_cps = size(cps, 2);

% define the affine transform for two masks
g.aff = cell(0);


% [A, t] = get_A_and_t(-15, [94; 141], [10; -10], [0.8,0.8])
[A, t] = get_A_and_t(-30, [94; 141], [0; 0], [1 - 0.2*manual_ratio, 1 - 0.2*manual_ratio])
[L, v] = get_Lv_from_At(A, t);
g.aff{end+1}.A = A;
g.aff{end}.t = t;
g.aff{end}.L = L;
g.aff{end}.v = v;
g.aff{end}.s = s;


[A, t] = get_A_and_t(0, [332;142], [-60 * 1; -20], [1, 1])
[L, v] = get_Lv_from_At(A, t);
g.aff{end+1}.A = A;
g.aff{end}.t = t;
g.aff{end}.L = L;
g.aff{end}.v = v;
g.aff{end}.s = s;


g.dim = dim;
g.nb_cps = nb_cps;
g.cps = cps;
g.Ilabel = Ilabel;

