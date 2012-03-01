function cps_desired_target = get_affine_on_cps(g)

cps_desired_target = zeros(g.dim, g.nb_cps);
for ii = 1:g.nb_cps
    cps_desired_target(:, ii) = g.aff{ii}.A * g.cps(:, ii) + g.aff{ii}.t;
end;