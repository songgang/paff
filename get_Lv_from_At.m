function [L, v] = get_Lv_from_At(A, t)
Lv = logm([A, t; 0 0 1]);
L = Lv(1:2, 1:2);
v = Lv(1:2, 3);

