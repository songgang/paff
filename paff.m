function solve_paff

xlist = 0:0.1:4;
% xlist = [0 1];
ylist = {};
tlist = {};
for x = xlist
    [T, Y] = ode45(@vhat_paff, [0, 1], [x 1 2]);   
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

r1 = log(3);
r2 = log(2);

y = yhat(1);
p = yhat(2);
q = yhat(3);

w1 = gpaff(y-p);
w2 = gpaff(y-q);
wsum = w1+w2;
w1 = w1 / wsum;
w2 = w2 / wsum;
v1 = r1 * y;
v2 = r2 * y;
vhat = zeros(3, 1);
vhat(1) = w1*v1+w2*v2;

w1 = gpaff(p-p);
w2 = gpaff(p-q);
wsum = w1+w2;
w1 = w1 / wsum;
w2 = w2 / wsum;
v1 = r1 * p;
v2 = r2 * p;
vhat(2) = w1*v1+w2*v2;

w1 = gpaff(q-p);
w2 = gpaff(q-q);
wsum = w1+w2;
w1 = w1 / wsum;
w2 = w2 / wsum;
v1 = r1 * q;
v2 = r2 * q;
vhat(3) = w1*v1+w2*v2;



% vhat(2) = v1;
% vhat(3) = v2;


function w = gpaff(d)

s = 0.2;
w = 1/(sqrt(2*pi) * s) * exp(- d*d / (s*s*2));

% ad = abs(d);
% 
% if (abs(d) <= 0.1)
%     w = 1;
% elseif (abs(d) <= 0.2)
%     w = (1 - ad * 4.5);
% else
%     w = 0.1;
% end;
% 
% w = sqrt(w);
% w = sqrt(w);
    
