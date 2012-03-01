function w = loggpaff(d, s)

% s = 80;
% w = - log(sqrt(2*pi) * s) + (- sum(d.*d, 1) ./ (s*s*2));
w = - (sum(d.*d, 1)./ (s*s*2)).^0.5;
