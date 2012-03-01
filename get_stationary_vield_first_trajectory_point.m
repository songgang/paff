
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

