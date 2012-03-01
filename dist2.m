
% A = [p1; p2; p3], each row is a vector
function d2 = dist2(A, B)

nb_A = size(A, 1);
nb_B = size(B, 1);
d2 = zeros(nb_A, nb_B);

dim = size(A, 2);

for d = 1:dim
    tmp = A(:, d) * ones(1, nb_B) - ones(nb_A, 1) * B(:, d)';
    d2 = d2 + tmp.*tmp;
end;
