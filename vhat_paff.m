function vhat=vhat_paff(g, t, yhat)

% use column vector for point coordinates vector to favor affine matrix operation
% vhat must be a column vector

cps = reshape(yhat, [g.dim, 1, g.nb_cps]);
vhat = zeros(length(yhat), 1);

for ii = 1:g.nb_cps
    vhat((ii-1)*g.dim+1 : ii*g.dim) = v_paff_ex_pqvec(g, t, yhat((ii-1)*g.dim+1 : ii*g.dim), cps);
end;

return;
