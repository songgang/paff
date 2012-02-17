function solve_paff_adative_time_2d_with_exp

leftB = -20;
rightB = 20;
bottomB = -20;
topB = 20;

[X, Y] = meshgrid(-40:1:40, -40:1:40);
xfield_0 = cat(3, X, Y);

% control point
p = [-5; 0];
q = [5; 0];
r = 1;

% solve iniitial ode on p and q
pqhat0 = [p; p; q];
tlist = 0:0.01:1;

% 
% 
% xlist = 0:0.1:6;
% p = 1,
% q = 5;
% tlist = 0:0.01:1;
% pqhat0 = [p, p, q];
% r = 0.5;

[T, pqhat] = ode45(@vhat_paff, tlist, pqhat0);


plist = [pqhat(:, 3), pqhat(:, 4)];
qlist = [pqhat(:, 5), pqhat(:, 6)];



% use matlab ode to solve the field

% use compose diff to solve the field

ind_t1 = 1;
clist = [];
clist(end+1) = ind_t1;
while (ind_t1 < length(tlist))
   ind_t2_a = find(dist2_rowvec_to_single(plist(ind_t1:end, :), qlist(ind_t1, :)) < r, 1, 'first');
   ind_t2_b = find(dist2_rowvec_to_single(qlist(ind_t1:end, :), plist(ind_t1, :)) < r, 1, 'first');
   
   if (isempty(ind_t2_a)) 
       ind_t2_a = length(tlist);
   else
       ind_t2_a = ind_t2_a + ind_t1-1;
   end;
   
   if (isempty(ind_t2_b))
       ind_t2_b = length(tlist);
   else
       ind_t2_b = ind_t2_b + ind_t1-1;

   end;
   
 
   ind_t2 = min(ind_t2_a, ind_t2_b);
   clist(end+1) = ind_t2;
   ind_t1 = ind_t2 + 1;
   ind_t1
end;


% y = phi(x) = x + offset(x)

yfield_current = xfield_0;


pslist = zeros(length(clist)-1, 2);
qslist = zeros(length(clist)-1, 2);

for ii = 1:length(clist)-1
    
    
    ind_t1 = clist(ii);
    ind_t2 = clist(ii+1);
    
    pa = plist(ind_t1, :);
    qa = qlist(ind_t1, :);
    
    % vfield = get_stationary_vield(xfield_0, pa', qa');
    vfield = get_stationary_vield_copy_paste(xfield_0, plist(ind_t1:ind_t2, :), qlist(ind_t1:ind_t2, :), tlist(ind_t1:ind_t2));
    yfield_delta = exp_mapping(vfield, X, Y, tlist(ind_t2)-tlist(ind_t1), 5);
    yfield_current = compose_phi(yfield_current, yfield_delta, X, Y);
    
    idx_p = find(X==p(1) & Y==p(2));
    idx_q = find(X==q(1) & Y==q(2));
    y1 = yfield_current(:,:,1);
    y2 = yfield_current(:,:,2);

    pslist(ii, :) = [y1(idx_p), y2(idx_p)];
    qslist(ii, :) = [y1(idx_q), y2(idx_q)];
    
end;
    


% figure(1); clf;
% hold on;
% plot(plist(:, 1), plist(:, 2), '-r*');
% plot(qlist(:, 1), qlist(:, 2), '-g*');
% % plot(pqhat(:, 1), tlist, '-b');
% 
% 
% plot(xfield_0(:, :, 1), xfield_0(:, :, 2), 'b.'); 
% plot(yfield_current(:, :, 1), yfield_current(:, :, 2), 'k.');
% 
% 
% quiver(xfield_0(:, :, 1), xfield_0(:, :, 2), yfield_current(:, :, 1)-xfield_0(:, :, 1), yfield_current(:, :, 2)-xfield_0(:, :, 2), 0);
% 
% 
% hold off;
% axis([leftB, rightB, bottomB, topB]);



figure(2); clf;
hold on;

plot(plist(:, 1), plist(:, 2), '-r*');
plot(qlist(:, 1), qlist(:, 2), '-g*');
% plot(xfield_0(:, :, 1), xfield_0(:, :, 2), 'b.'); 
% plot(yfield_current(:, :, 1), yfield_current(:, :, 2), 'k.');
hold off;

pad=20;
fil=1;

meshplot(X(pad*fil:fil:end-pad*fil, pad*fil:fil:end-pad*fil), Y(pad*fil:fil:end-pad*fil, pad*fil:fil:end-pad*fil), 'Color', 'r');
meshplot(yfield_current(pad*fil:fil:end-pad*fil, pad*fil:fil:end-pad*fil, 1), yfield_current(pad*fil:fil:end-pad*fil, pad*fil:fil:end-pad*fil, 2), 'Color', 'b');

idx_p = find(X==p(1) & Y==p(2));
idx_q = find(X==q(1) & Y==q(2));
y1 = yfield_current(:,:,1);
y2 = yfield_current(:,:,2);

hold on;
quiver(p(1), p(2), y1(idx_p)-p(1), y2(idx_p)-p(2), 0, 'r', 'LineWidth', 2);
quiver(q(1), q(2), y1(idx_q)-q(1), y2(idx_q)-q(2), 0, 'g', 'LineWidth', 2);

for ii = 1:size(pslist, 1)

    quiver(p(1), p(2), pslist(ii,1)-p(1), pslist(ii,2)-p(2), 0, 'r', 'LineWidth', 2);
    quiver(q(1), q(2), qslist(ii,1)-q(1), qslist(ii,2)-q(2), 0, 'g', 'LineWidth', 2);


end;



hold off;

axis equal;


figure(3); clf;
plot(clist, 'b*');
title('clist');

% y, p, q are column vectors
function v = v_paff_ex(t, y, p, q)
     
scale = 100;

[A1, t1] = get_A_and_t(0, [-10;0], [20;0]);
[A2, t2] = get_A_and_t(15, [10; 0], [0;0]);

[L1, v1] = get_Lv_from_At(A1, t1);
[L2, v2] = get_Lv_from_At(A2, t2);

% L1 = [0 -0.02; 0.02 0] * 100; scale;
% v1 = [1; 0];
% L2 = [0 0.015; -0.015 0] * 0; scale;
% v2 = [3; 0];

nb = size(y, 2);


logw1 = loggpaff(y-p*ones(1, nb));
logw2 = loggpaff(y-q*ones(1, nb));

w1 = 1 ./ (1 + exp(logw2 - logw1));
w2 = 1 ./ (1 + exp(logw1 - logw2));

% w1 = gpaff() + 0.01;
% w2 = gpaff(y-q*ones(1, nb)) + 0.01;
% wsum = w1+w2;
% w1 = w1 ./ wsum;
% w2 = w2 ./ wsum;



v1 = L1 * y + v1*ones(1, nb);
v2 = L2 * y + v2*ones(1, nb);

%   w2 = 0.9 * ones(1,nb);
%   w1 = 0.1 * ones(1, nb);


v = ([1;1]*w1).*v1 + ([1;1]*w2).*v2;

if isnan(v)
    disp 'haha';
end;
 




function solve_paff

xlist = 1:0.1:4;
% xlist = [0 1];
ylist = {};
tlist = {};
for x = xlist
    [T, Y] = ode45(@vhat_paff, [0, 1], [x 2 3]);   
    ylist{end+1} = Y(:, 1);
    tlist{end+1} = T(:);
    
end;

figure(2); clf;
for ii = 1:size(ylist, 2);
    hold on;
    plot(ylist{ii}, tlist{ii}, '-k');
    hold off;
    
    
end;


hold on;
plot(Y(:, end-1), T(:), '-c');
plot(Y(:, end), T(:), '-m');
hold off;
grid on;
    




function vhat=vhat_paff(t, yhat)

% use column vector for point coordinates vector to favor affine matrix operation
y = [yhat(1); yhat(2)];
p = [yhat(3); yhat(4)];
q = [yhat(5); yhat(6)];

% vhat must be a column vector
vhat = zeros(6, 1);

vhat(1:2) = v_paff_ex_pqvec(t, y, p, q);
vhat(3:4) = v_paff_ex_pqvec(t, p, p, q);
vhat(5:6) = v_paff_ex_pqvec(t, q, p, q);


function w = loggpaff(d)

s = 2;
w = - log(sqrt(2*pi) * s) + (- sum(d.*d, 1) ./ (s*s*2));


function w = gpaff(d)

s = 0.1;
w = 1/(sqrt(2*pi) * s) .* exp(- sum(d.*d, 1) ./ (s*s*2));

% ad = abs(d);
% 
% if (abs(d) <= 0.1)
%     w = 1;
% elseif (abs(d) <= 0.2)
%     w = (1 + 0.1*4.5 - ad * 4.5);
% else
%     w = 0.1;
% end;
% 
% w = sqrt(w);
% w = sqrt(w);
    
function d2 = dist2_rowvec_to_single(plist, q)

nb_p = size(plist, 1);
d = plist - ones(nb_p, 1) * q;
d2 = sum(d.*d, 2);


function vfield = get_stationary_vield(xfield, p, q)

nb_yaxis = size(xfield, 1);
nb_xaxis = size(xfield, 2);
X = xfield(:, :, 1);
Y = xfield(:, :, 2);

% v is indep of t, set t = nan
v = v_paff_ex(nan, [X(:)'; Y(:)'], p, q);
vfield(:, :, 1) = reshape(v(1, :), [nb_yaxis, nb_xaxis]);
vfield(:, :, 2) = reshape(v(2, :), [nb_yaxis, nb_xaxis]);


function vfield = get_stationary_vield_copy_paste(xfield, plist, qlist, tlist)


nb_yaxis = size(xfield, 1);
nb_xaxis = size(xfield, 2);
X = xfield(:, :, 1);
Y = xfield(:, :, 2);

% points are row vectors
D = dist2([X(:), Y(:)], [plist; qlist]);
[non, idx_min] = min(D, [], 2);
nb_t = length(tlist);
idx_min(idx_min > nb_t) = idx_min(idx_min > nb_t) - nb_t;



% v is indep of t, set t = nan
v = v_paff_ex_pqvec(nan, [X(:)'; Y(:)'], plist(idx_min, :)', qlist(idx_min, :)');
vfield(:, :, 1) = reshape(v(1, :), [nb_yaxis, nb_xaxis]);
vfield(:, :, 2) = reshape(v(2, :), [nb_yaxis, nb_xaxis]);


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

function [A, t] = get_A_and_t(theta, c, t0)

A = [cos(theta/180*pi), -sin(theta/180*pi); sin(theta/180*pi), cos(theta/180*pi)];
t = (eye(2) - A) * c + t0;


% y, p, q are column vectors
function v = v_paff_ex_pqvec(t, y, p, q)
     
scale = 100;
[A1, t1] = get_A_and_t(0, [-5;0], [20;0]);
[A2, t2] = get_A_and_t(40, [10;0], [0;0]);
[L1, v1] = get_Lv_from_At(A1, t1);
[L2, v2] = get_Lv_from_At(A2, t2);

% L1 = [0 -0.02; 0.02 0] * 100; scale;
% v1 = [1; 0];
% L2 = [0 0.015; -0.015 0] * 0; scale;
% v2 = [3; 0];
nb = size(y, 2);

logw1 = loggpaff(y - p );
logw2 = loggpaff(y - q );

w1 = 1 ./ (1 + exp(logw2 - logw1));
w2 = 1 ./ (1 + exp(logw1 - logw2));


v1y = L1 * y + v1*ones(1, nb);
v2y = L2 * y + v2*ones(1, nb);

v = ([1;1]*w1).*v1y + ([1;1]*w2).*v2y;

if isnan(v)
    disp 'haha';
end;

