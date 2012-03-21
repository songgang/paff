function demo_transform_objects

% warping the lung lobes

I = double(imread('data/01_Fixed_y161ud.png'));
Ilabel = double(imread('data/01_Fixed_y161ud_label.png'));

figure(16); clf;
imagesc(Ilabel); axis image;
print -dpng -f16 'data/01_Fixed_y161ud_lable_colored.png';


leftB = 1;
rightB = size(I, 2);
bottomB = 1;
topB = size(I, 1);


% define the boundary boxes
% [dim1a, dim1b; dim2a, dim2b];
g.boundary.box = [leftB, rightB; bottomB, topB];
g.boundary.s = 10;

[X, Y] = meshgrid(leftB:1:rightB, bottomB:1:topB);

% test on the Gaussian smoothing for velocity field
% this is used for decreasing from boundaries
g.vfield_smooth_sigma = 0.1;

% g.s2 is to combine weight
% g.s2 = 1;

g = get_control_point_and_transform_multibeam_lobes(g);


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

figure(13); clf;
imagesc(I); axis image; colormap gray;


for ii = 1: length(clist)-1
    
    ind_t1 = clist(ii);
    ind_t2 = clist(ii+1);
    % create velocity field using RBF decreasing function
    vfield = get_stationary_vield_copy_paste_decreasing(g, xfield_0, cpslist(ind_t1:ind_t2, :, :), alpha);
    % vfield = smooth_field(vfield, 1, 'Gaussian');
    
    %     mask_label = scatter_multiple_label_trajectory_after_distance_transform(g, cpslist(ind_t1:ind_t2, :, :));
    %     vfield_inside_mask_accurate = get_stationary_vfield_from_labeled_mask(g, xfield_0, mask_label);
    %     vfield = smooth_field(vfield_inside_mask_accurate, 20, 'PrecondtionVariationalWithBoundary', mask_label, vfield);
    
    
    
    
    
    q=2;
%     figure(13);
% %     hold on;
% %     quiver(X(1:q:end,1:q:end), ...
% %         Y(1:q:end,1:q:end), ...
% %         vfield(1:q:end,1:q:end, 1) / max(vfield(:)) * 4, ...
% %         vfield(1:q:end,1:q:end, 2)  / max(vfield(:)) * 4, 1);
% %     hold off;
% %     
%     hold on;
%     for kk = 1:nb_cps;
%         plot(cpslist(ind_t1:ind_t2, 1, kk), cpslist(ind_t1:ind_t2, 2, kk), 'r.-', 'MarkerSize', 1);
%         % quiver(cpslist(ind_t1:ind_t2, 1, kk), cpslist(ind_t1:ind_t2, 2, kk), vcpslist(ind_t1:ind_t2, 1 ,kk) * 0.02, vcpslist(ind_t1:ind_t2, 2 ,kk) * 0.02, 1, 'r');
%     end;
%     hold off;
    %
    
    
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


    figure(13);
%     hold on;
%     quiver(X(1:q:end,1:q:end), ...
%         Y(1:q:end,1:q:end), ...
%         vfield(1:q:end,1:q:end, 1) / max(vfield(:)) * 4, ...
%         vfield(1:q:end,1:q:end, 2)  / max(vfield(:)) * 4, 1);
%     hold off;
%     
    hold on;
    for kk = 1:nb_cps;
        plot(cpslist(:, 1, kk), cpslist(:, 2, kk), 'r.-', 'MarkerSize', 5);
        % quiver(cpslist(ind_t1:ind_t2, 1, kk), cpslist(ind_t1:ind_t2, 2, kk), vcpslist(ind_t1:ind_t2, 1 ,kk) * 0.02, vcpslist(ind_t1:ind_t2, 2 ,kk) * 0.02, 1, 'r');
    end;
    hold off;

    print -dpng -f13 'data/01_Fixed_y161ud_warp_field.png';

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





axis([leftB, rightB, bottomB, topB]);
axis image;
axis ij;

%
% figure(3); clf;
% plot(clist, tlist(clist), 'b*');
% title('clist');







Iwarped = griddata(yfield_current(:, :, 1), yfield_current(:, :, 2), I, X, Y);
figure(14); clf;
imagesc(Iwarped); axis image; colormap gray;
imwrite(uint8(Iwarped), 'data/01_Fixed_y161ud_warped.png');



function g = get_control_point_and_transform_multibeam_lobes(g)


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


[A, t] = get_A_and_t(-15, [94; 141], [10; -10], [0.8,0.8])
[L, v] = get_Lv_from_At(A, t);
g.aff{end+1}.A = A;
g.aff{end}.t = t;
g.aff{end}.L = L;
g.aff{end}.v = v;
g.aff{end}.s = s;


[A, t] = get_A_and_t(30, [332;142], [-20; -20], [1, 0.8])
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

