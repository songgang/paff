function solve_paff_orig

xlist = -1:0.25:2;
ylist = {};
tlist = {};
for x = xlist
    [T, Y] = ode45(@v_paff, [0, 1], [x]);
    ylist{end+1} = Y(:, 1);
    tlist{end+1} = T;
end;

figure(3); clf;
for ii = 1:size(ylist, 2);
    hold on;
    plot(ylist{ii}, tlist{ii}, '-b');
    hold off;
    
    
end;


hold on;
idx = find(xlist==0);
plot(ylist{idx}, tlist{idx}, '-r');
idx = find(xlist==1);
plot(ylist{idx}, tlist{idx}, '-g');
hold off;




grid on;
    


function v=v_paff(t, y)

v1 = 1;
v2 = -1;

w1 = gpaff(y-(0+v1*t));
w2 = gpaff(y-(1+v2*t));

wsum = w1+w2;
w1 = w1 / wsum;
w2 = w2 / wsum;


v = zeros(1, 1);
v = w1*v1+w2*v2;


function w = gpaff(d)

s = 0.1;
w = 1/(sqrt(2*pi) * s) * exp(- d*d / (s*s*2));