% detection collision of multiple trajectory
% by if trajectories are overlapped in the spatial domain
%
% input: cpslist, nb_time_points * spatial_dim * nb_trajectory
%
% output: clist, a 1*K vector, of index in cpslist first dimension(time
% points), clist(ii) to clist(ii+1) is continuous 

function clist = collision_detection(g, cpslist)

% nb_T = size(cpslist, 1);
% nb_cps = size(cpslist, 3);

nb_T = g.nb_T;
nb_cps = g.nb_cps;
r = g.r;

clist(1) = 1;
for ii = 2:nb_T;
    is_collision = 0;
    for jj = 1:nb_cps
        for kk = jj+1:nb_cps
            
            d2 = dist2(cpslist(clist(end):ii, :, jj), cpslist(clist(end):ii, :, kk));
            if (min(d2(:)) < r)
                clist(end+1) = ii;
                is_collision = 1;
                break;
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

% clist = clist(1:20:end);
