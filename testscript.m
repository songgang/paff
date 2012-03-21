function testscript

% f = myGaussianList2D(NaN, [10 20 30 40]);
% figure; 
% imagesc(f);
% axis image;
% colorbar;
% 
% figure(100);hold on;
% plot(f((end+1)/2, :), 'b*');

myGaussianMod2D(1/20);
myGaussianMod2D(1/4);
myGaussianMod2D(1/2);
myGaussianMod2D(1/1);

function G = myGaussianList2D(image, sigma_list)

maxsigma = min(sigma_list);

a = ceil(3*maxsigma+1);
[X, Y] = ndgrid(-a:a, -a:a);

G = zeros(size(X));

for ii = 1:length(sigma_list)
    sigma = sigma_list(ii);

    G = G + exp( -1 *(X.*X+Y.*Y) / (2*sigma*sigma));
end;

minsigma = min(sigma_list);
G = G * 1/(2*pi*minsigma)/ max(G(:));



function G = myGaussianMod2D(sigma)


a1 = load('durer', 'X');
I = a1.X;


a = ceil(3*1/sigma)+1;
szK = [a, a];

% [X, Y] = ndgrid(-a:a, -a:a);
% 
% G = zeros(size(X));
% G =  exp( -1 *(X.*X+Y.*Y) / (2*sigma*sigma));
% % G((end+1)/2-1:(end+1)/2+1, (end+1)/2-1:(end+1)/2+1) = [ 0 -1 0; -1 4 -1; 0 -1 0];

% fG(isnan(fG)) = 0;

% ifG = real(ifftn(fG));
% ifG(abs(ifG)>0.01) = 0;



szf = size(I) + szK - 1;

% Compute FFT of the image
fft_I = fftn(I, szf);

% Compute FFT of the kernel
% fft_K = fftn(K, szf);

cX = (szf(2)+1)/2;
cY = (szf(1)+1)/2;
lX = 1 - cX;
rX = szf(2) - cX;
lY = 1 - cY;
rY = szf(1) - cY;

[X,Y] = meshgrid(lX:rX, lY:rY);


% fG1 = (fftn(G, szf));
fG1 =  exp( -1 *(X.*X+Y.*Y) / (2/sigma/sigma));
fG = fG1 ./ (1- 0.95 * fG1 );
% fG = 1e-8./(fG1+1e-10);
fft_K1 = fG;

fft_K = fftshift(fft_K1);

fft_K = ones(size(fft_K));

% figure; imagesc(real(ifftn(fft_K)));
% figure; imagesc((fG)); axis image;

% Multiply
fft_I = fft_I .* fft_K;

% Save memory
clear('fft_K');

% Compute inverse
F = real(ifftn(fft_I));

% Take the central part
%shift = floor(szK/2);
%F = F(shift(1)+1:end-shift(1),shift(2)+1:end-shift(2));

F = F(1:size(I, 1), 1:size(I, 2));



% figure; clf; imagesc(I); colormap gray;
figure; clf; subplot(1,2,1);imagesc(F); colormap gray; colorbar;
% subplot(1,2,2); imagesc(I-F); colormap gray;







% 
% 
% 
% 
% figure; imagesc(G);
% figure; imagesc(abs(fG));
% figure; imagesc(ifG);