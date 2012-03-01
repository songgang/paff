function filtered = myGaussianLPF(image, sigma)

a = ceil(3*sigma+1);
[X, Y, Z] = ndgrid(-a:a, -a:a, -a:a);

G = 1./((sqrt(2*pi)*sigma).^3) .* exp( -1 *(X.*X+Y.*Y+Z.*Z) / (2*sigma*sigma));

filtered = imfilter3d(image, G);