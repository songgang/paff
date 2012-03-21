I = double(imread('data/r16slice.png'));

theta = -30 * pi/ 180;
A = [1.2*cos(theta), -sin(theta); sin(theta), 0.8*cos(theta)];
c = sz.' / 2;
t0 = [0; 0];
t = -1 * A * c + c+ t0;

Ainv = inv(A);
tinv = -Ainv * t;




sz = size(I);

[X0, Y0] = meshgrid(1:sz(2), 1:sz(1));

x0 = [X0(:), Y0(:)].';
oneb = ones(1,size(x,2));
x1 = A * x0 + t * oneb;
X1 = reshape(x1(1, :), size(X0));
Y1 = reshape(x1(2, :), size(X1));

x0 = [X0(:), Y0(:)].';
oneb = ones(1,size(x,2));
x1inv = Ainv * x0 + tinv * oneb;
X1inv = reshape(x1inv(1, :), size(X0));
Y1inv = reshape(x1inv(2, :), size(X1));


I1 = interp2(X, Y, I, x1(1, :), x1(2, :), '*Linear', 0);
I1 = reshape(I1, sz);

imwrite(uint8(I1), 'data/r16slice_affined.png');

figure(1); clf;
imagesc(I); axis image; colormap gray;

figure(2); clf;
imagesc(I1); axis image; colormap gray;


fil=10;
figure(3); clf;
meshplot(X1inv(1:fil:end, 1:fil:end), Y1inv(1:fil:end, 1:fil:end), 'Color', 'b'); 
axis image; axis ij;