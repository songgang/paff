% test_scattered_interpolation

% a1 = load('cpslist_for_interp.mat');
% cpslist = a1.cpslist; % number_points X dim=2 X number_of_trajectories
% 
% [X,Y] = meshgrid(-150:150, -150:150);
% 
% ZI = griddata(,y,z,XI,YI)
leftB = -150;
rightB = 150;
bottomB = -150;
topB = 150;

dim = 2;

g.boundary.box = [leftB, rightB; bottomB, topB];
g.boundary.s = 2;

[X, Y] = meshgrid(leftB:1:rightB, bottomB:1:topB);

% test on the Gaussian smoothing for velocity field
g.vfield_smooth_sigma = 2;

% control point
% column vector
% cps = [-10 10 0   0;
%         0  0  10 -10 ];

%      cps = [-10 0;
%               0 30];
cps = [-50 0;
    0 0];


nb_cps = size(cps, 2);
g.aff = cell(0);

% s is the sigma to computer affine velocity decreasing
s = 2;
% g.s2 is to combine weight
g.s2 = 10;

[A, t] = get_A_and_t(30, [0;0], [50; 0], [1, 1]);
[L, v] = get_Lv_from_At(A, t);
g.aff{end+1}.A = A;
g.aff{end}.t = t;
g.aff{end}.L = L;
g.aff{end}.v = v;
g.aff{end}.s = s;

[A, t] = get_A_and_t(30, [0; 0], [0; 50], [1,1]);
[L, v] = get_Lv_from_At(A, t);
g.aff{end+1}.A = A;
g.aff{end}.t = t;
g.aff{end}.L = L;
g.aff{end}.v = v;
g.aff{end}.s = s;


% [A, t] = get_A_and_t(0, [0; 0], [0;0], [1,1]);
% [L, v] = get_Lv_from_At(A, t);
% g.aff{end+1}.A = A;
% g.aff{end}.t = t;
% g.aff{end}.L = L;
% g.aff{end}.v = v;
% g.aff{end}.s = s;

% [A, t] = get_A_and_t(0, [0; -10], [0;20], [1,1]);
% [L, v] = get_Lv_from_At(A, t);
% g.aff{4}.A = A;
% g.aff{4}.t = t;
% g.aff{4}.L = L;
% g.aff{4}.v = v;
% g.aff{4}.s = s;

% [A, t] = get_A_and_t(0, [-10;0], [40;0], [2, 1]);
% [L, v] = get_Lv_from_At(A, t);
% g.aff{1}.A = A;
% g.aff{1}.t = t;
% g.aff{1}.L = L;
% g.aff{1}.v = v;
% g.aff{1}.s = s;
%
% [A, t] = get_A_and_t(0, [0; 10], [0;-40], [1,1]);
% [L, v] = get_Lv_from_At(A, t);
% g.aff{2}.A = A;
% g.aff{2}.t = t;
% g.aff{2}.L = L;
% g.aff{2}.v = v;
% g.aff{2}.s = s;
%

g.dim = dim;
g.nb_cps = nb_cps;
g.cps = cps;

r = 100;

xfield_0 = cat(3, X, Y);

% solve iniitial ode on p and q
% pqhat0 = [p; q];
cpshat0 = reshape(cps, [nb_cps*dim, 1]);
tlist = 0:0.01:1;

[T, cpshat] = ode45(@(t, y) vhat_paff(g, t, y), tlist, cpshat0);
nb_T = length(T);
cpslist = reshape(cpshat, [nb_T, dim, nb_cps]); % [pqhat(:, 3), pqhat(:, 4)];


g.nb_T = nb_T;


% recompute the velocity at each point
% cpslist: nb_pts * nb_dim * nb_cps (idx_point, u/v, idx_cps)
cpslist1 = permute(cpslist, [2, 1, 3]);
cpslist2 = repmat(cpslist1, [1, nb_cps, 1]);
vcpslist1 = v_paff_ex_pqvec(g, nan, reshape(cpslist1, [dim, nb_T*nb_cps]), cpslist2);
vcpslist = reshape(vcpslist1, [dim, nb_T, nb_cps]);
vcpslist = permute(vcpslist, [2, 1, 3]);

vfield = scatter_interploate_trajectory_nearest_neighbor([leftB, rightB, bottomB, topB], cpslist, vcpslist);

figure(12); clf;
% plot trajectory of control points (cps) using ode solution
hold on;
clrs='gbr';
for ii = 1:nb_cps
    x1 = squeeze(cpslist(:, 1, ii));
    y1 = squeeze(cpslist(:, 2, ii));
    clrr = clrs(mod(ii, length(clrs))+1);
    plot(squeeze(cpslist(:, 1, ii)), squeeze(cpslist(:, 2, ii)), [clrr, '.']);
    u1 = squeeze(vcpslist(:, 1, ii)) * 0.01;
    v1 = squeeze(vcpslist(:, 2, ii)) *0.01;
    quiver(x1, y1, u1, v1, 0, clrr);
end;
hold off;

hold on;
plot(X, Y, '.g');
quiver(X, Y, vfield(:, :, 1) * 0.01, vfield(:, :, 2) * 0.01, 0, 'g');
hold off;

axis equal;
axis([leftB, rightB, bottomB, topB]);