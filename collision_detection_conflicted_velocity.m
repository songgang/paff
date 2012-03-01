function clist = collision_detection_conflicted_velocity(g, cpslist, vcpslist)

% if two trajectories are very close AND the velocities on close points are
% not equal (within a threshold), then it is conflict

nb_T = g.nb_T;
nb_cps = g.nb_cps;
r = g.r;
dv = g.dv; % difference of velocity

clist(1) = 1;
for ii = 2:nb_T;
    is_collision = 0;
    for jj = 1:nb_cps
        for kk = jj+1:nb_cps
            
            offset = clist(end);
            d2 = dist2(cpslist(clist(end):ii, :, jj), cpslist(clist(end):ii, :, kk));
            
            [mind, mind_idx] = min(d2(:));
            
            if (mind < r)
                
                % check if velocities are different
                
                [J, K] = ind2sub(size(d2), mind_idx);
                J = J + offset - 1;
                K = K + offset - 1;
                
                cpslist(J, :, jj)
                cpslist(K, :, kk)
                vcpslist(K, :, kk)
                vcpslist(J, :, jj)

                
                vd2 = dist2(vcpslist(J, :, jj), vcpslist(K, :, kk))
                if vd2 > dv
                    clist(end+1) = ii;
                    is_collision = 1;
                    break;
                end;
            end;
            
        end;
        
        if (is_collision)
            break;
        end;
    end;
    
    if (clist(end) == nb_T)
        break;
    end;
    
end;

if (clist(end) ~= nb_T)
    clist(end+1) = nb_T;
end;

