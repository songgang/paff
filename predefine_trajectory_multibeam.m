function [cpslist, vcpslist] = predefine_trajectory_multibeam(g, tlist)
% tlist must start from zero

cps = g.cps;
nb_cps = g.nb_cps;
dim = g.dim;
nb_T = length(tlist);

cpslist = zeros(nb_T, dim, nb_cps);
vcpslist = zeros(nb_T, dim, nb_cps);

cpslist(1, :, :) = reshape(cps, [1, dim, nb_cps]);
for ii = 1:nb_T
    for jj = 1:nb_cps
        affL = g.aff{g.ind(jj)}.L;
        affv = g.aff{g.ind(jj)}.v;
        y = squeeze(cpslist(ii, :, jj));
        vy = affL * y' + affv;
        vcpslist(ii, :, jj) = reshape(vy, [1, dim, 1]);
        if ii < nb_T
            cpslist(ii+1, :, jj) = cpslist(ii, :, jj) + vcpslist(ii, :, jj) * (tlist(ii+1)-tlist(ii));
        end;
    end;
end;