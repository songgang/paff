function vsmooth = smooth_field(v, sigma, varargin)

if nargin <= 2
    method = 'Gaussian';
else
    method = varargin{1};
end;


if sigma == 0
    vsmooth = v;
    return;
end;

switch method
    case 'Gaussian'
        dim = size(v, 3);
        vsmooth =zeros(size(v));
        for ii = 1:dim
            vsmooth(:, :, ii) = myGaussianLPF2D(v(:, :, ii), sigma);
        end;
    case 'Variational'
        dim = size(v, 3);
        vsmooth =zeros(size(v));
        vnorm = sqrt(v(:,:,1).^2+v(:,:,2).^2);
        for ii = 1:dim
            
            v1 = v(:, :, ii);
            mv1 = max((vnorm(:)));
            vsmooth(:, :, ii) = myVariational2D(vnorm/mv1, v1 / mv1, sigma);
            vsmooth(:, :, ii) = vsmooth(:, :, ii)  * mv1;
%             vsmooth(:, :, ii) = myVariational2D(vnorm, v(:, :, ii), sigma);
        end;
        
    case 'PrecondtionVariational'
        dim = size(v, 3);
        vsmooth =zeros(size(v));
        vnorm = sqrt(v(:,:,1).^2+v(:,:,2).^2);
        vnorm(vnorm > 0) = 1;
        for ii = 1:dim
            vsmooth(:, :, ii) = myPreconditionVariational2D(vnorm, v(:, :, ii), sigma);
        end;
        
    case 'PrecondtionVariationalWithBoundary'
        boundary = varargin{2};
        v0 = varargin{3};
        dim = size(v, 3);
        vsmooth =zeros(size(v));
        for ii = 1:dim
            vsmooth(:, :, ii) = myPreconditionVariationalWithBoundary2D(boundary, v(:, :, ii), sigma, v0(:, :, ii));
        end;
       
end;


function v1 = myVariational2D(w, v, sigma, ddim)

v1 = zeros(size(v));
% v1 = myGaussianLPF2D(v, sigma);
% v1 = v*0.1 + v1*0.9;

maxw = max(w(:));

for ii = 1:100
     a = -2 * w .* (v - v1) * 1 * 0.5; %  / maxw;
% a = -2 *  (v - v1) / 1;  % maxw;
%     a1 = myGaussianLPF2D(a, sigma);

% a1 = myGaussianLPF2D(a, sigma);
a1 = myGaussianTableLPF2D(a, sigma);


% a1 = mySimpleLPF2D(a, sigma);
      v1 = v1 - a1;

%     b = myGaussianLPF2D(v1, sigma);
%     g = a + 2 * b;
%     v1 = v1 - g;


%     if ddim == 1
%     
%     figure; imagesc(v1); colorbar;
%     end;
end;

for ii = 1:0
    a = -2 * w .* (v - v1);
      a1 = myGaussianLPF2D(a, 1);
% a1 = mySimpleLPF2D(a, sigma);
     v1 = v1 - a1;

%     b = myGaussianLPF2D(v1, sigma);
%     g = a + 2 * b;
%     v1 = v1 - g;


%     if ddim == 1
%     
%     figure; imagesc(v1); colorbar;
%     end;
end;




function filtered = mySimpleLPF2D(image, sigma)

% a = ceil(3*sigma+1);
% [X, Y] = ndgrid(-a:a, -a:a);
% 
% G = 1./((sqrt(2*pi)*sigma).^2) .* exp( -1 *(X.*X+Y.*Y) / (2*sigma*sigma));

G = [0 1 0;1 4 1; 0 1 0];
G = G / sum(G(:));

filtered = imfilter2d(image, G);


function v1 = myPreconditionVariational2D(w, v, sigma, ddim)

v1 = zeros(size(v));
maxw = max(w(:));

idx_tr = (w ~= 0);
se = strel('disk', ceil(sigma+1) );



for ii = 1:5
     g1 = -2 * w .* (v - v1); %  / maxw;
     g = myGaussianLPF2D(g1, sigma);
     
     idx_gnonzero = (g~=0);
     if isempty(find(idx_gnonzero))
         disp 'converge';
         break; % converge
     end;
     
     
     idx_tr1 = idx_tr & idx_gnonzero;
     l_tr =  ( v1(idx_tr1) - v(idx_tr1) ) ./ (g(idx_tr1));
 
     % use rbf
     [y, x] = ind2sub(size(v), find(idx_tr1));
     [X, Y] = meshgrid(1:size(v, 2), 1:size(v, 1));
     lrbf = rbf_interpolation(l_tr', x, y, 1, X, Y);
     % lrbf(idx_tr1) = l_tr * 0.9 + lrbf(idx_tr1) * 0.1;
     v1 = v1 - lrbf .* g;
     
%      % use dilation
%      l = zeros(size(v));
%      l(idx_tr) = l_tr;
%      ldilated = imdilate(l, se);
%     v1 = v1 - ldilated .* g;

     v1(141,124)
     
     figure(100); clf;imagesc(v1);
end;

return;











% x,y are row vectors
% X, Y are from meshgrid
% g is a row vector

function f = rbf_interpolation(g, x, y, sigma, X, Y)

nb_cps = size(x, 1);
nb_y = prod(size(X));

logwlist = zeros(nb_cps, nb_y);
for ii = 1:nb_cps
    logwlist(ii,:) = loggpaff(transpose([X(:), Y(:)] - ones(nb_y, 1) * [x(ii), y(ii)]), sigma);
end;

wlist = exp(logwlist);
sw = sum(wlist, 1);
f = g * wlist ./ sw;
f = reshape(f ,size(X));

