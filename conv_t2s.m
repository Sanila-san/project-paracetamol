% �������� T-���������� � S-���������

function S = conv_t2s(T)

DT = T(1, 1) * T(2, 2) - T(1, 2) * T(2, 1);
S = 1 / T(2, 2) * [T(1, 2) DT; 1 -T(2, 1)];

return