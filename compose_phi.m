% compute phi = \phi2 \circ \phi1 (X)
function phi = compose_phi(phi1, phi2, X, Y)

phi(:, :, 1) = interp2(X, Y, phi2(:,:,1), phi1(:,:,1), phi1(:,:,2), 'linear', 1000); 
phi(:, :, 2) = interp2(X, Y, phi2(:,:,2), phi1(:,:,1), phi1(:,:,2), 'linear', 1000); 

phi = remove_out_of_bound_with_orig_value(phi, phi1, X, Y);

function phi = remove_out_of_bound_with_orig_value(phi, phi1, X, Y)

phix = phi(:, :, 1);
phi1x = phi1(:, :, 1);
idx = (phix > max(X(:))) | (phix < min(X(:))) ;
phix(idx) = phi1x(idx); 
phi(:, :, 1) = phix;


phiy = phi(:, :, 2);
phi1y = phi1(:, :, 2);
idx = (phiy > max(Y(:))) | (phiy < min(Y(:))) ;
phiy(idx) = phi1y(idx); 
phi(:, :, 2) = phiy;
