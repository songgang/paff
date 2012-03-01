function vfield = get_stationary_vield_copy_paste_average_multiple_trajectory(g, xfield, cpslist, tlist)

nb_yaxis = size(xfield, 1);
nb_xaxis = size(xfield, 2);
X = xfield(:, :, 1);
Y = xfield(:, :, 2);

% points are row vectors
% for each trajectory find the closet time
% then sum up (not average) all trajectory

nb_t = length(tlist);
v = zeros(g.dim, length(X(:)));
normv = sum(v.*v, 1);
for kk = 1:g.nb_cps
    D = dist2([X(:), Y(:)], cpslist(:, :, kk));
    [non, idx_minA] = min(D, [], 2);
    nb_t = length(tlist);
    idx_trajectory = floor((idx_minA-1) /  nb_t) + 1;
    idx_min = mod(idx_minA-1, nb_t) + 1;
    
    cps_for_X = permute(cpslist(idx_min, :, :), [2,1,3]);
    v1 = v_paff_ex_pqvec(g, nan, [X(:)'; Y(:)'], cps_for_X);
    normv1 = sum(v1.*v1, 1);
    idx1= find(normv1 > normv);
    v(:, idx1) = v1(:, idx1);
    normv = sum(v.*v, 1);
    
    
end;

vfield(:, :, 1) = reshape(v(1, :), [nb_yaxis, nb_xaxis]);
vfield(:, :, 2) = reshape(v(2, :), [nb_yaxis, nb_xaxis]);


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
