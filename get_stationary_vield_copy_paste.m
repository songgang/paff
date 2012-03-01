function vfield = get_stationary_vield_copy_paste(g, xfield, cpslist, tlist)
% function vfield = get_stationary_vield_copy_paste(xfield, plist, qlist, tlist)


nb_yaxis = size(xfield, 1);
nb_xaxis = size(xfield, 2);
X = xfield(:, :, 1);
Y = xfield(:, :, 2);

% points are row vectors
nb_t = length(tlist);
D = dist2([X(:), Y(:)], reshape(permute(cpslist, [1,3,2]), [nb_t * g.nb_cps, g.dim, ]) );
[non, idx_minA] = min(D, [], 2);
nb_t = length(tlist);
idx_trajectory = floor((idx_minA-1) /  nb_t) + 1;
idx_min = mod(idx_minA-1, nb_t) + 1;

% figure; clf;
% imagesc(reshape(idx_minA, size(X)));
% axis xy;
%
% figure; clf;
% imagesc(reshape(idx_min, size(X)));
% axis xy;
%
% figure; clf;
% imagesc(reshape(idx_trajectory, size(X)));
% axis xy;
%
% figure; imagesc([-150, 150], [-150, 150], reshape(idx_minA, size(X)));
% hold on; plot(cpslist(:, 1, 2), cpslist(:, 2, 2), 'y*'); hold off;
% hold on; plot(cpslist(:, 1, 1), cpslist(:, 2, 1), 'y*'); hold off;
% axis xy;

% v is indep of t, set t = nan
% use column vector for control points positions at time t(idx_min)

cps_for_X = permute(cpslist(idx_min, :, :), [2,1,3]);
v = v_paff_ex_pqvec(g, nan, [X(:)'; Y(:)'], cps_for_X); %    plist(idx_min, :)', qlist(idx_min, :)');
vfield(:, :, 1) = reshape(v(1, :), [nb_yaxis, nb_xaxis]);
vfield(:, :, 2) = reshape(v(2, :), [nb_yaxis, nb_xaxis]);
