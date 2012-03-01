% function solve_paff_adative_time_2d_with_exp

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
g.vfield_smooth_sigma = 20;

% g.s2 is to combine weight
g.s2 = 20;

g = get_control_point_and_transform(g);
nb_cps = g.nb_cps;
cps = g.cps;

dim = g.dim;

% r is the collision boundary
r = 400;
g.r = r;

% dv is the difference of the velocity when conflicted
g.dv = 1;


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

clist = collision_detection_conflicted_velocity(g, cpslist, vcpslist);



xfield_0 = cat(3, X, Y);
yfield_current = xfield_0;

cpspiecelist = zeros(length(clist), 2, nb_cps);
% qslist = zeros(length(clist)-1, 2);

for ii = 1:nb_cps
    cpspiecelist(1, :, ii) = cps(:, ii)';
end;

for ii = 1: length(clist)-1
    
    ind_t1 = clist(ii);
    ind_t2 = clist(ii+1);
    
    % vfield = get_stationary_vield_copy_paste(g, xfield_0, cpslist(ind_t1:ind_t2, :, :), tlist(ind_t1:ind_t2));
    % vfield = smooth_field(vfield, g.vfield_smooth_sigma);
    
    vfield = scatter_interploate_trajectory_nearest_neighbor([leftB, rightB, bottomB, topB], cpslist(ind_t1:ind_t2, :, :), vcpslist(ind_t1:ind_t2, :, :));
    vfield = smooth_field(vfield, g.vfield_smooth_sigma, 'Variational');
    
        q=10;
        figure;
        quiver(X(1:q:end,1:q:end), ...
            Y(1:q:end,1:q:end), ...
            vfield(1:q:end,1:q:end, 1) / max(vfield(:)) * 4, ...
        vfield(1:q:end,1:q:end, 2)  / max(vfield(:)) * 4, 1);
    
        hold on;
        for kk = 1:nb_cps;
%             plot(cpslist(ind_t1:ind_t2, 1, kk), cpslist(ind_t1:ind_t2, 2, kk), 'r.');
            quiver(cpslist(ind_t1:ind_t2, 1, kk), cpslist(ind_t1:ind_t2, 2, kk), vcpslist(ind_t1:ind_t2, 1 ,kk) * 0.02, vcpslist(ind_t1:ind_t2, 2 ,kk) * 0.02, 1, 'r');
        end;
        hold off;
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
            1,  clrs(mod(jj, length(clrs))+1), 'LineWidth', 2);
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




















