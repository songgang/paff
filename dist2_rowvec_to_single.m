function d2 = dist2_rowvec_to_single(plist, q)

nb_p = size(plist, 1);
d = plist - ones(nb_p, 1) * q;
d2 = sum(d.*d, 2);
