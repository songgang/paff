function solve_paff_adative_time


xlist = 0:0.1:6;
p = 1,
q = 5;
tlist = 0:0.01:1;
pqhat0 = [p, p, q];
r = 0.5;

[T, pqhat] = ode45(@vhat_paff, tlist, pqhat0);


plist = pqhat(:, 2);
qlist = pqhat(:, 3);


figure(1); clf;
hold on;
plot(plist, tlist, '-r*');
plot(qlist, tlist, '-g*');
% plot(pqhat(:, 1), tlist, '-b');
hold off;

ind_t1 = 1;
clist = [];
clist(end+1) = ind_t1;
while (ind_t1 < length(tlist))
   ind_t2_a = find(abs(plist(ind_t1+1:end)  - qlist(ind_t1)) < r, 1, 'first');
   ind_t2_b = find(abs(qlist(ind_t1+1:end)  - plist(ind_t1)) < r, 1, 'first');
   
   if (isempty(ind_t2_a)) 
       ind_t2_a = length(tlist);
   else
       ind_t2_a = ind_t2_a + ind_t1;
   end;
   
   if (isempty(ind_t2_b))
       ind_t2_b = length(tlist);
   else
       ind_t2_b = ind_t2_b + ind_t1;

   end;
   
   
 
    ind_t2 = min(ind_t2_a, ind_t2_b);
   clist(end+1) = ind_t2;
   ind_t1 = ind_t2;
   ind_t1
end;

disp clist;


xlist_last = xlist;
ylist = zeros(length(plist), length(xlist_last));

for ii = 1:length(clist)-1
    xlist_current = xlist_last;
    
    ind_t1 = clist(ii);
    ind_t2 = clist(ii+1);
    
    p = plist(ind_t1);
    q = qlist(ind_t1);
    
    for ind_x = 1:length(xlist)
        x = xlist_current(ind_x);
        
        % tlist1 = zeros(1, (ind_t2 - ind_t1 + 1) * 5 + 1);
        % tlist1 = tlist(ind_t1) : (tlist(ind_t1+1) - tlist(ind_t1)) / 5 : tlist(ind_t2);
        
        % [T, Y1] = ode45(@(t, y) v_paff_ex(t, y, p, q), tlist1, x);   
        
        
        % Y = Y1(1:5:end);
        
        [T, Y] = ode45(@(t, y) v_paff_ex(t, y, p, q), tlist(ind_t1:ind_t2), x);   
        
        if length(T) ~= ind_t2 - ind_t1 + 1
            disp 'haha';
        end;
        
        if (ind_t2 > ind_t1 + 1)
            ylist(ind_t1:ind_t2, ind_x) = Y;
        else
            ylist(ind_t1, ind_x) = Y(1);
            ylist(ind_t2, ind_x) = Y(end);
        end;
    end;
    xlist_last = ylist(ind_t2, :);
end;


figure(1); 
for ind_x = 1:length(xlist)
    hold on;
    plot(ylist(:, ind_x), tlist, 'b-');
    hold off;
end;







function v = v_paff_ex(t, y, p, q)
        
w1 = gpaff(y-p) + 0.01;
w2 = gpaff(y-q) + 0.01;
wsum = w1+w2;
w1 = w1 / wsum;
w2 = w2 / wsum;
v1 = log(3) * y;
v2 = log(2) * y;

v = w1*v1+w2*v2;

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

y = yhat(1);
p = yhat(2);
q = yhat(3);

vhat = zeros(3, 1);

vhat(1) = v_paff_ex(t, y, p, q);
vhat(2) = v_paff_ex(t, p, p, q);
vhat(3) = v_paff_ex(t, q, p, q);

% 
% w1 = gpaff(y-p);
% w2 = gpaff(y-q);
% wsum = w1+w2;
% w1 = w1 / wsum;
% w2 = w2 / wsum;
% v1 = r1 * y;
% v2 = r2 * y;
% vhat = zeros(3, 1);
% vhat(1) = w1*v1+w2*v2;
% 
% w1 = gpaff(p-p);
% w2 = gpaff(p-q);
% wsum = w1+w2;
% w1 = w1 / wsum;
% w2 = w2 / wsum;
% v1 = r1 * p;
% v2 = r2 * p;
% vhat(2) = w1*v1+w2*v2;
% 
% w1 = gpaff(q-p);
% w2 = gpaff(q-q);
% wsum = w1+w2;
% w1 = w1 / wsum;
% w2 = w2 / wsum;
% v1 = r1 * q;
% v2 = r2 * q;
% vhat(3) = w1*v1+w2*v2;


% vhat(2) = v1;
% vhat(3) = v2;


function w = gpaff(d)

s = 10;
w = 1/(sqrt(2*pi) * s) * exp(- d*d / (s*s*2));

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
    
