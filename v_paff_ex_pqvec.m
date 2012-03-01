% y, p, q are column vectors
% cps is a 3D array, each slice in z is for one control point
% cps(nb_dim, nb_y, id_cps)
% return v as column vector, each column is a velocity vector for the point
%  at same position in y 
function v = v_paff_ex_pqvec(g, t, y, cps)

nb_y = size(y, 2);

% add support from boundary lines
% each row is all y for one control point
% each column is one y for all control points
logwlist = zeros(g.nb_cps, nb_y);
for ii = 1:g.nb_cps
    logwlist(ii,:) = loggpaff(y - squeeze(cps(:, :, ii)), g.aff{ii}.s);
end;

logwboundlist = zeros(g.dim*2, nb_y);
for ii = 1:g.dim
    logwboundlist(ii*2-1, :) = loggpaff(y(ii, :) - g.boundary.box(ii,1), g.boundary.s);
    logwboundlist(ii*2, :) = loggpaff(y(ii, :) - g.boundary.box(ii,2), g.boundary.s);
end;

if 1
    
    % weighted by gaussian distance
    v = zeros(g.dim, nb_y);
    for ii = 1:g.nb_cps
        sw = zeros(1, nb_y);
        for jj = 1:g.nb_cps
            sw = sw + exp(logwlist(jj, :) - logwlist(ii, :));
        end;
        wii = 1 ./ sw;
        viiy = g.aff{ii}.L * y + g.aff{ii}.v*ones(1, nb_y);
        v = v + (ones(g.dim, 1)*wii).*viiy;
    end;
    % scaled by minimum of nearest distances to each trajectory
    if 1
        logw2list = zeros(g.nb_cps, nb_y);
        for ii = 1:g.nb_cps
            logw2list(ii,:) = loggpaff(y - squeeze(cps(:, :, ii)), g.s2);
        end;
        
        [maxlogwlist, idx_cps_idx] = max(logw2list, [], 1);
        v = v.* (ones(g.dim, 1) * exp(maxlogwlist));
    end;
    
else
    
    v = zeros(g.dim, nb_y);
    [maxlogwlist, idx_cps_idx] = max(logwlist, [], 1);
    for ii = 1:g.nb_cps
        sw = zeros(1, nb_y);
        sw = exp(logwlist(ii, :));
        wii = sw .* (idx_cps_idx == ii);
        viiy = g.aff{ii}.L * y + g.aff{ii}.v*ones(1, nb_y);
        v = v + (ones(g.dim, 1)*wii).*viiy;
    end;
    
end;


if 1 % boundary condition
    % velocity is scaled by distance to the fixed boundary
    swratio = -1 * inf(1, nb_y);
    for jj = 1:g.dim*2
        %         sw = sw + exp(logwboundlist(jj, :) - logwlist(ii, :));
        swratio = max(swratio, logwboundlist(jj, :));
    end;
    swratio = 1 - exp(swratio);
    v = v.* (ones(g.dim, 1) * swratio);
end;

% logw1 = loggpaff(y - p );
% logw2 = loggpaff(y - q );
%
% w1 = 1 ./ (1 + exp(logw2 - logw1));
% w2 = 1 ./ (1 + exp(logw1 - logw2));
%
%
% v1 = L1 * y + v1*ones(1, nb);
% v2 = L2 * y + v2*ones(1, nb);
%
% v = ([1;1]*w1).*v1 + ([1;1]*w2).*v2;

if isnan(v)
    disp 'haha';
end;

