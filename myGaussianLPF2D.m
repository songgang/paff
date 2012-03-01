function filtered = myGaussianLPF2D(image, sigma)

a = ceil(3*sigma+1);
[X, Y] = ndgrid(-a:a, -a:a);

G = 1./(2*pi*sigma.^2) .* exp( -1 *(X.*X+Y.*Y) / (2*sigma*sigma));
% G = fspecial('gaussian', 2*a+1, sigma);

filtered = imfilter2d(image, G);
% filtered = imfilter(image, G);