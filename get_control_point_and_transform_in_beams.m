function g = get_control_point_and_transform_multibeam(g)


% control point
% column vector
% cps = [-10 10 0   0;
%         0  0  10 -10 ];
%      cps = [-10 0;
%               0 30];

% s is the sigma to computer affine velocity decreasing
% s can be removed in future implementation
s = 2;
dim = 2;

% simulation of two masks
% mask1 
[X, Y] = meshgrid(-60:10:-40, -10:10:10);
cps1 = [X(:)'; Y(:)']; % dim * nb_pts
ind1 = ones(1, size(cps1, 1));
cps = cps1;
ind = ind1;

% mask2 
[X, Y] = meshgrid(-5:10:5, -10:10:20);
cps1 = [X(:)'; Y(:)']; % dim * nb_pts
ind1 = 2 * ones(1, size(cps1, 1));
cps = [cps, cps1];
ind = [ind, ind1];


nb_cps = size(cps, 2);


% define the affine transform for two masks
g.aff = cell(0);

[A, t] = get_A_and_t(0, [0;0], [50; 0], [1, 1]);
[L, v] = get_Lv_from_At(A, t);
g.aff{end+1}.A = A;
g.aff{end}.t = t;
g.aff{end}.L = L;
g.aff{end}.v = v;
g.aff{end}.s = s;

[A, t] = get_A_and_t(0, [0; 0], [0; 50], [1,1]);
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
