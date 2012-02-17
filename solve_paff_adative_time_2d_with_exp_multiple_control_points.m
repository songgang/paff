function solve_paff_adative_time_2d_with_exp

figno = 8;

leftB = -40;
rightB = 40;
bottomB = -40;
topB = 40;

dim = 2;

% define the boundary boxes
% [dim1a, dim1b; dim2a, dim2b];
g.boundary.box = [leftB, rightB; bottomB, topB];
g.boundary.s = 3;

[X, Y] = meshgrid(leftB:1:rightB, bottomB:1:topB);

% control point
% column vector
% cps = [-10 10 0   0;  
%         0  0  10 -10 ];

%      cps = [-10 0;  
%               0 30];

cps = [-10 0;  
        0  10];


nb_cps = size(cps, 2);
g.aff = cell(0);

s = 10;

[A, t] = get_A_and_t(0, [0;0], [20; 0], [1, 1]);
[L, v] = get_Lv_from_At(A, t);
g.aff{end+1}.A = A;
g.aff{end}.t = t;
g.aff{end}.L = L;
g.aff{end}.v = v;
g.aff{end}.s = s;

% [A, t] = get_A_and_t(0, [0; 0], [20;0], [1,1]);
% [L, v] = get_Lv_from_At(A, t);
% g.aff{end+1}.A = A;
% g.aff{end}.t = t;
% g.aff{end}.L = L;
% g.aff{end}.v = v;
% g.aff{end}.s = s;


[A, t] = get_A_and_t(45, [0; 20], [0;-20], [1,1]);
[L, v] = get_Lv_from_At(A, t);
g.aff{end+1}.A = A;
g.aff{end}.t = t;
g.aff{end}.L = L;
g.aff{end}.v = v;
g.aff{end}.s = s;

% [A, t] = get_A_and_t(0, [0; -10], [0;20], [1,1]);
% [L, v] = get_Lv_from_At(A, t);
% g.aff{4}.A = A;
% g.aff{4}.t = t;
% g.aff{4}.L = L;
% g.aff{4}.v = v;
% g.aff{4}.s = s;

% [A, t] = get_A_and_t(0, [-10;0], [40;0], [2, 1]);
% [L, v] = get_Lv_from_At(A, t);
% g.aff{1}.A = A;
% g.aff{1}.t = t;
% g.aff{1}.L = L;
% g.aff{1}.v = v;
% g.aff{1}.s = s;
% 
% [A, t] = get_A_and_t(0, [0; 10], [0;-40], [1,1]);
% [L, v] = get_Lv_from_At(A, t);
% g.aff{2}.A = A;
% g.aff{2}.t = t;
% g.aff{2}.L = L;
% g.aff{2}.v = v;
% g.aff{2}.s = s;
% 

g.dim = dim;
g.nb_cps = nb_cps;
g.cps = cps;

r = 2;




xfield_0 = cat(3, X, Y);

% solve iniitial ode on p and q
% pqhat0 = [p; q];
cpshat0 = reshape(cps, [nb_cps*dim, 1]);
tlist = 0:0.01:1;

[T, cpshat] = ode45(@(t, y) vhat_paff(g, t, y), tlist, cpshat0);
nb_T = length(T);
cpslist = reshape(cpshat, [nb_T, dim, nb_cps]); % [pqhat(:, 3), pqhat(:, 4)];


g.nb_T = nb_T;

% use matlab ode to solve the field
% use compose diff to solve the field

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





yfield_current = xfield_0;

cpspiecelist = zeros(length(clist), 2, nb_cps);
% qslist = zeros(length(clist)-1, 2);

for ii = 1:nb_cps
    cpspiecelist(1, :, ii) = cps(:, ii)';
end;

for ii = 1:length(clist)-1
    
    ind_t1 = clist(ii);
    ind_t2 = clist(ii+1);
    
    vfield = get_stationary_vield_copy_paste(g, xfield_0, cpslist(ind_t1:ind_t2, :, :), tlist(ind_t1:ind_t2));
    yfield_delta = exp_mapping(vfield, X, Y, tlist(ind_t2)-tlist(ind_t1), 10);
    yfield_current = compose_phi(yfield_current, yfield_delta, X, Y);
    
    
    for k = 1:nb_cps
        idx_p = find(X==cps(1,k) & Y==cps(2,k));
        y1 = yfield_current(:,:,1);
        y2 = yfield_current(:,:,2);
        cpspiecelist(ii+1, :, k) = [y1(idx_p), y2(idx_p)];
    end;
    
end;
    

pad=0;
fil=2;

figure(figno); clf;
% plot trajectory of control points (cps) using ode solution
hold on;
clrs='gbr';
for ii = 1:nb_cps
    plot(squeeze(cpslist(:, 1, ii)), squeeze(cpslist(:, 2, ii)), ['-', clrs(mod(ii, length(clrs))+1), '*']);
end;
hold off;

% meshplot(X(pad*fil:fil:end-pad*fil, pad*fil:fil:end-pad*fil), Y(pad*fil:fil:end-pad*fil, pad*fil:fil:end-pad*fil), 'Color', 'g');
% meshplot(yfield_current(pad*fil:fil:end-pad*fil, pad*fil:fil:end-pad*fil, 1), yfield_current(pad*fil:fil:end-pad*fil, pad*fil:fil:end-pad*fil, 2), 'Color', 'b');

meshplot(X(pad*fil+1:fil:end-pad*fil, pad*fil+1:fil:end-pad*fil), Y(pad*fil+1:fil:end-pad*fil, pad*fil+1:fil:end-pad*fil), 'Color', 'g');
meshplot(yfield_current(pad*fil+1:fil:end-pad*fil, pad*fil+1:fil:end-pad*fil, 1), yfield_current(pad*fil+1:fil:end-pad*fil, pad*fil+1:fil:end-pad*fil, 2), 'Color', 'b');

hold on;
for jj = 1:nb_cps
    for ii = size(cpspiecelist, 1):size(cpspiecelist, 1)
        quiver(cpspiecelist(1,1,jj), ...
               cpspiecelist(1,2,jj), ...
               cpspiecelist(ii,1,jj)-cpspiecelist(1,1,jj), ...
               cpspiecelist(ii,2,jj)-cpspiecelist(1,2,jj), ...
               0,  clrs(mod(jj, length(clrs))+1), 'LineWidth', 2);
    end;
end;
hold off;




% plot desired position
cps_desired_target = get_affine_on_cps(g);
hold on;
for jj = 1:nb_cps
        quiver(cps(1,jj), ...
               g.cps(2,jj), ...
               cps_desired_target(1,jj)-cps(1,jj), ...
               cps_desired_target(2,jj)-cps(2,jj), ...
               0, 'c', 'LineWidth', 2);
end;
hold off;
axis equal;



axis equal;




figure(3); clf;
plot(clist, 'b*');
title('clist');



function vhat=vhat_paff(g, t, yhat)

% use column vector for point coordinates vector to favor affine matrix operation
% vhat must be a column vector

cps = reshape(yhat, [g.dim, 1, g.nb_cps]); 
vhat = zeros(length(yhat), 1);

for ii = 1:g.nb_cps
    vhat((ii-1)*g.dim+1 : ii*g.dim) = v_paff_ex_pqvec(g, t, yhat((ii-1)*g.dim+1 : ii*g.dim), cps);
end;

return;


function w = loggpaff(d, s)

% s = 80;
% w = - log(sqrt(2*pi) * s) + (- sum(d.*d, 1) ./ (s*s*2));
w = (- sum(d.*d, 1) ./ (s*s*2));

    
function d2 = dist2_rowvec_to_single(plist, q)

nb_p = size(plist, 1);
d = plist - ones(nb_p, 1) * q;
d2 = sum(d.*d, 2);


function vfield = get_stationary_vield_copy_paste(g, xfield, cpslist, tlist)
% function vfield = get_stationary_vield_copy_paste(xfield, plist, qlist, tlist)


nb_yaxis = size(xfield, 1);
nb_xaxis = size(xfield, 2);
X = xfield(:, :, 1);
Y = xfield(:, :, 2);

% points are row vectors
nb_t = length(tlist);
D = dist2([X(:), Y(:)], reshape(permute(cpslist, [1,3,2]), [nb_t * g.nb_cps, g.dim, ]) );
[non, idx_min] = min(D, [], 2);
nb_t = length(tlist);
idx_min = mod(idx_min, nb_t) + 1;


% v is indep of t, set t = nan
% use column vector for control points positions at time t(idx_min)

cps_for_X = permute(cpslist(idx_min, :, :), [2,1,3]);
v = v_paff_ex_pqvec(g, nan, [X(:)'; Y(:)'], cps_for_X); %    plist(idx_min, :)', qlist(idx_min, :)');
vfield(:, :, 1) = reshape(v(1, :), [nb_yaxis, nb_xaxis]);
vfield(:, :, 2) = reshape(v(2, :), [nb_yaxis, nb_xaxis]);

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





function [L, v] = get_Lv_from_At(A, t)
Lv = logm([A, t; 0 0 1]);
L = Lv(1:2, 1:2);
v = Lv(1:2, 3);

function [A, t] = get_A_and_t(theta, c, t0, s)

A = [cos(theta/180*pi), -sin(theta/180*pi); sin(theta/180*pi), cos(theta/180*pi)];
A = A * diag(s);
t = (eye(2) - A) * c + t0;


% y, p, q are column vectors
% cps is a 3D array, each slice in z is for one control point
% cps(nb_dim, nb_y, id_cps)  
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


v = zeros(g.dim, nb_y);

% weighted by gaussian distance
if 0
for ii = 1:g.nb_cps
    sw = zeros(1, nb_y);
    for jj = 1:g.nb_cps
        sw = sw + exp(logwlist(jj, :) - logwlist(ii, :));
    end;
    wii = 1 ./ sw;
    viiy = g.aff{ii}.L * y + g.aff{ii}.v*ones(1, nb_y);
    v = v + (ones(g.dim, 1)*wii).*viiy;
end;
else
% scaled by minimum of nearest distances to each trajectory
[non, idx_cps_idx] = max(logwlist, [], 1);
for ii = 1:g.nb_cps
    sw = zeros(1, nb_y);
    sw = exp(logwlist(ii, :));
    wii = sw .* (idx_cps_idx == ii);
    viiy = g.aff{ii}.L * y + g.aff{ii}.v*ones(1, nb_y);
    v = v + (ones(g.dim, 1)*wii).*viiy;
end;

end;


% velocity is scaled by distance to the fixed boundary
swratio = -1 * inf(1, nb_y);
for jj = 1:g.nb_cps
%         sw = sw + exp(logwboundlist(jj, :) - logwlist(ii, :));
    swratio = max(swratio, logwboundlist(jj, :));
end;
swratio = 1 - exp(swratio);
v = v.* (ones(g.dim, 1) * swratio);

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

function cps_desired_target = get_affine_on_cps(g)

cps_desired_target = zeros(g.dim, g.nb_cps);
for ii = 1:g.nb_cps
    cps_desired_target(:, ii) = g.aff{ii}.A * g.cps(:, ii) + g.aff{ii}.t;
end;
