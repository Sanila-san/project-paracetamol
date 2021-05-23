% �������� �������� ���������� ����� � S-���������
% (L - ����� �����; Z0 - ������� ��������)

function S = conv_gz2s(gamma, Z, l, Z0)

G = (Z - Z0) / (Z + Z0);
S1 = [G 1 - G; 1 + G -G];
S2 = [0 exp(-gamma * l); exp(-gamma * l) 0];
S3 = [-G 1 + G; 1 - G G];
S = conv_t2s(conv_s2t(S1) * conv_s2t(S2) * conv_s2t(S3));

return