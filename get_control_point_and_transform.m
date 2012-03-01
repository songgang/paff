function g = get_control_point_and_transform(g)


% control point
% column vector
% cps = [-10 10 0   0;
%         0  0  10 -10 ];

%      cps = [-10 0;
%               0 30];

% s is the sigma to computer affine velocity decreasing
s = 2;
dim = 2;

cps = [-50 0;
    0 0];

nb_cps = size(cps, 2);
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
