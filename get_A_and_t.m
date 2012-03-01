function [A, t] = get_A_and_t(theta, c, t0, s)

A = [cos(theta/180*pi), -sin(theta/180*pi); sin(theta/180*pi), cos(theta/180*pi)];
A = A * diag(s);
t = (eye(2) - A) * c + t0;
