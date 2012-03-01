% test diffusion style variational smooth on velocity field
% try to minimize the functional in the Soblev space
% E(h) = \int |g| (g-h).^2 dx
% s.t. h \in Soblev space
% 
% updating scheme:
%  h = h - Soblev gradient of E(h)
%    = h - Gaussian ( -2* |g| * (g-h))

v = zeros(100, 100 ,2);
v(50, 10:55, 1) = 1;
v(50, 10:55, 2) = 0;

v(10:50, 60, 1) = 0;
v(10:50, 60, 2) = -1;

% v=vfield;

sigma = 10;

% maxv = max(v(:));
% vsmooth = smooth_field(v / maxv, sigma, 'Variational');
% vsmooth = vsmooth * maxv;

vsmooth = smooth_field(v, sigma, 'Variational');

figure(1); clf;

hold on;
[X, Y] = meshgrid(1:5:size(v, 2), 1:5:size(v, 1));
meshplot(X, Y, 'Color', 'g');
hold off;

hold on;
quiver(v(:,:,1), v(:, :, 2), 0, 'b');
quiver(vsmooth(:, :, 1), vsmooth(:, :, 2), 0, 'r'); 

hold off;

axis image;
axis equal;

figure; imagesc(sqrt(v(:, :, 1).^2+v(:,:,2).^2)); colorbar;
figure; imagesc(sqrt(vsmooth(:, :, 1).^2+vsmooth(:,:,2).^2)); colorbar;