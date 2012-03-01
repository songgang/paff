function filtered = myGaussianTableLPF2D(image, sigma)

s1 = sigma;
r0 = sigma;

a = ceil(3*(s1+r0)+1);
[X, Y] = ndgrid(-a:a, -a:a);


r2 = r0*r0;

d2 = X.*X +Y.*Y;
G1 = 1./(2*pi*s1.^2) .* exp( -1 *(d2 - r0.*r0) / (2*s1*s1));
GH=G1;
GH(d2<=r2) = 1./(2*pi*s1.^2);





% G = 1./(2*pi*s1.^2) .* exp( -1 *(X.*X+Y.*Y) / (2*s1*s1));
% G = fspecial('gaussian', 2*a+1, sigma);

% b = ceil(12*sigma + 1);
% H = ones(b, b);
% H = H / sum(H(:));
% 
% GH = conv2(G, H, 'full');










GH = GH / sum(GH(:));

filtered = imfilter2d(image, GH);
