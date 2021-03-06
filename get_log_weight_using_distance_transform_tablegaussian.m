function logw = get_log_weight_using_distance_transform_tablegaussian(mask, h, sigma)

d = bwdist(mask);
d = max(0, d-h);

logw = - 1/(sigma*sigma) .* (d.*d);
% logw = -1 * log( 2*pi*sigma^2) + logw;