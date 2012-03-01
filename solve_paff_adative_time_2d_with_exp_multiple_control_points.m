function solve_paff_adative_time_2d_with_exp

figno = 11;

leftB = -150;
rightB = 150;
bottomB = -150;
topB = 150;

dim = 2;

% define the boundary boxes
% [dim1a, dim1b; dim2a, dim2b];
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

% use matlab ode to solve the field
% use compose diff to solve the field

clist(1) = 1;
for ii = 2:nb_T;
    is_collision = 0;
    for jj = 1:nb_cps
        for kk = jj+1:nb_cps
            
            d2 = dist2(cpslist(clist(end):ii, :, jj), cpslist(clist(end):ii, :, kk));
            if (min(d2(:)) < r)
                clist(end+1) = ii;
                is_collision = 1;
                break;
            end;
            
        end;
        
        if (is_collision)
            break;
        end;
    end;
    
    if (clist(end) == nb_T)
        break;
    end;
    
end;

if (clist(end) ~= nb_T)
    clist(end+1) = nb_T;
end;

% clist = clist(1:20:end);



yfield_current = xfield_0;

cpspiecelist = zeros(length(clist), 2, nb_cps);
% qslist = zeros(length(clist)-1, 2);

for ii = 1:nb_cps
    cpspiecelist(1, :, ii) = cps(:, ii)';
end;

for ii = 1:length(clist)-1
    
    ind_t1 = clist(ii);
    ind_t2 = clist(ii+1);
    
    vfield = get_stationary_vield_copy_paste(g, xfield_0, cpslist(ind_t1:ind_t2, :, :), tlist(ind_t1:ind_t2));
    vfield = smooth_field(vfield, g.vfield_smooth_sigma);
    
    % vfield = get_stationary_vield_copy_paste_average_multiple_trajectory(g, xfield_0, cpslist(ind_t1:ind_t2, :, :), tlist(ind_t1:ind_t2));
    % vfield = get_stationary_vield_first_trajectory_point(g, xfield_0);
    %
    %     q=5;
    %     figure;
    %     quiver(X(1:5:end,1:5:end), ...
    %         Y(1:5:end,1:5:end), ...
    %         vfield(1:5:end,1:5:end, 1) / max(vfield(:)) * 4, ...
    %     vfield(1:5:end,1:5:end, 2)  / max(vfield(:)) * 4, 0);
    %
    %     hold on;
    %     for kk = 1:nb_cps;
    %         plot(cpslist(ind_t1:ind_t2, 1, kk), cpslist(ind_t1:ind_t2, 2, kk), 'r*-');;
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
fil=5;

figure(figno); clf;
% plot trajectory of control points (cps) using ode solution
hold on;
clrs='gbr';
for ii = 1:nb_cps
    plot(squeeze(cpslist(:, 1, ii)), squeeze(cpslist(:, 2, ii)), ['-', clrs(mod(ii, length(clrs))+1), '*']);
end;
hold off;

% meshplot(X(pad*fil:fil:end-pad*fil, pad*fil:fil:end-pad*fil), Y(pad*fil:fil:end-pad*fil, pad*fil:fil:end-pad*fil), 'Color', 'g');
% meshplot(yfield_current(pad*fil:fil:end-pad*fil, pad*fil:fil:end-pad*fil, 1), yfield_current(pad*fil:fil:end-pad*fil, pad*fil:fil:end-pad*fil, 2), 'Color', 'b');

meshplot(X(pad*fil+1:fil:end-pad*fil, pad*fil+1:fil:end-pad*fil), Y(pad*fil+1:fil:end-pad*fil, pad*fil+1:fil:end-pad*fil), 'Color', 'g');
meshplot(yfield_current(pad*fil+1:fil:end-pad*fil, pad*fil+1:fil:end-pad*fil, 1), yfield_current(pad*fil+1:fil:end-pad*fil, pad*fil+1:fil:end-pad*fil, 2), 'Color', 'b');

hold on;
for jj = 1:nb_cps
    for ii = size(cpspiecelist, 1):size(cpspiecelist, 1)
        quiver(cpspiecelist(1,1,jj), ...
            cpspiecelist(1,2,jj), ...
            cpspiecelist(ii,1,jj)-cpspiecelist(1,1,jj), ...
            cpspiecelist(ii,2,jj)-cpspiecelist(1,2,jj), ...
            0,  clrs(mod(jj, length(clrs))+1), 'LineWidth', 2);
    end;
end;
hold off;




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




figure(3); clf;
plot(clist, tlist(clist), 'b*');
title('clist');
















function cps_desired_target = get_affine_on_cps(g)

cps_desired_target = zeros(g.dim, g.nb_cps);
for ii = 1:g.nb_cps
    cps_desired_target(:, ii) = g.aff{ii}.A * g.cps(:, ii) + g.aff{ii}.t;
end;




function vfield = get_stationary_vield_first_trajectory_point(g, xfield)
% function vfield = get_stationary_vield_copy_paste(xfield, plist, qlist, tlist)

nb_yaxis = size(xfield, 1);
nb_xaxis = size(xfield, 2);
X = xfield(:, :, 1);
Y = xfield(:, :, 2);

% points are row vectors
cps = g.cps;
cps_fox_X = cps(:,:,ones(1, length(X(:))));
cps_for_X = permute(cps_fox_X, [1,3,2]);
% cps(nb_dim, nb_y, id_cps)


% v is indep of t, set t = nan
% use column vector for control points positions at time t(idx_min)


v = v_paff_ex_pqvec(g, nan, [X(:)'; Y(:)'], cps_for_X); %    plist(idx_min, :)', qlist(idx_min, :)');
vfield(:, :, 1) = reshape(v(1, :), [nb_yaxis, nb_xaxis]);
vfield(:, :, 2) = reshape(v(2, :), [nb_yaxis, nb_xaxis]);

