%% ‘ункци€ расчета параметров копланарного волновода:
% (the electromagnetic propagation properties of CPW)
% (по материалам статьи: Heinrich Quasi-TEM ... 1993)
% ¬ходные данные (частота, размеры и параметры материалов):
% f - частота; t - толщина металлического сло€, wg - ширина земли,
% s - зазор земл€-сигнал, w - ширина сигнальной линии;
% er - диэлектрическа€ проницаемость подложки,
% tan_d - тангенс угла потерь, k - проводимость металла.
% ¬ыходные данные:
% C, G, R, L - параметры линии на единицу длины

function [C, G, R, L] = CPW_EM_prop(f, t, wg, s, w, er, tan_d, k)

%  онстанты:
e0 = 8.85418781762039e-12; % электрическа€ посто€нна€
nu0 = 4 * pi * 1e-7; % магнитна€ посто€нна€


% Appendix A.1:
k0 = w / (w + 2 * s);
k1 = k0 * sqrt((1 - ((w + 2 * s) / (w + 2 * s + 2 * wg))^2) / ...
               (1 - (w / (w + 2 * s + 2 * wg))^2));
k2 = k0 * sqrt((1 - ((w + 2 * s) / (4 * w + 2 * s))^2) / ...
               (1 - (w / (4 * w + 2 * s))^2));

% Appendix A.2:
a = w / 2; b = w / 2 + s; tH = t / 2;
pc0 = b / (2 * a) * 1 / (CPW_spec_func('Ktouch', k0))^2;
pc1 = 1 + log(8 * pi * a / (a + b)) + a / (a + b) * log(b / a);
pc2 = pc1 - 2 * a / b * (CPW_spec_func('Ktouch', k0))^2;
pc3 = 2 * b^2 / (a * (b + a)) * CPW_spec_func('Etouch', k0) / CPW_spec_func('Ktouch', k0);
pc4 = (b - a) / (b + a) * (log(8 * pi * a / (a + b)) + a / b);
pc5 = (b - a) / (b + a) * log(3);
pc6 = (b - a) / (b + a) * log(24 * pi * b * (a + b) / (b - a)^2) - ...
       b / (b + a) * log(b / a);

% Appendix A.3:
arg = [tH s a b pc0 pc1 pc2 pc3 pc4 pc5 pc6];
FLc = CPW_spec_func('FLc', arg);
FLg = CPW_spec_func('FLg', arg);


%% II. Calculation of The Capacitance
arg = [t k1 pc0 pc1 pc2 s]; Fup = CPW_spec_func('F', arg);
Flow = CPW_spec_func('K', k1) / CPW_spec_func('Ktouch', k1);
C = 2 * e0 * (Fup + er * Flow); % формулы (2)
omega = 2 * pi * f; % G зависит от частоты:
G = 2 * omega * e0 * er * tan_d * Flow;

arg = [t / 2 k1 pc0 pc1 pc2 s]; F0 = CPW_spec_func('F', arg);
Le8 = nu0 / (4  * F0); % формула (3)


%% III. The Skin-Effect Values of R and Li
% формулы (4):
RcSE = sqrt(omega * nu0 / (2 * k)) * FLc / (4 * F0^2);
RgSE = sqrt(omega * nu0 / (2 * k)) * FLg / (4 * F0^2);
LiSE = sqrt(nu0 / (2 * omega * k)) * (FLc + FLg) / (4 * F0^2);


%% IV. The Line Resistance R
omegaC1 = sqrt(2) * 4 / (nu0 * k * t * w); % формулы (8)
omegaC2 = 8 / (nu0 * k) * ((w + t) / (w * t))^2;
omegaG1 = 2 / (nu0 * k * t * wg);
omegaG2 = 2 / (nu0 * k) * ((2 * wg + t) / (wg * t))^2;
Rc0 = 1 / (k * w * t);
Rc1 = sqrt(omegaC2 * nu0 / (2 * k)) * FLc / (4 * F0^2);
Rg0 = 1 / (2 * k * wg * t);
Rg1 = sqrt(omegaG2 * nu0 / (2 * k)) * FLg / (4 * F0^2);
vc = log(Rc0 / Rc1) / log(omegaC1 / omegaC2);
vg = log(Rg0 / Rg1) / log(omegaG1 / omegaG2);
% Appendix A.5:
gammaC = (omegaC1 / omegaC2)^2; % = (w * t / (sqrt(2) * (w + t)^2))^2
gammaG = (omegaG1 / omegaG2)^2; % = (wg * t / ((2 * wg + t)^2))^2

a4c = (gammaC * vc + 1 / 4 * (1 / 2 - vc) * (4 - vc * (1 - gammaC^2))) / ...
      (4 - vc - 1 / 4 * (1 / 2 - vc) * (4 - vc * (1 - gammaC^2)));
a3c = 1 / 4 * (1 / 2 - vc) * (1 + a4c);
a2c = 1 / gammaC * (a4c - a3c);
a1c = a2c + gammaC * a3c;
a4g = (gammaG * vg + 1 / 4 * (1 / 2 - vg) * (4 - vg * (1 - gammaG^2))) / ...
      (4 - vg - 1 / 4 * (1 / 2 - vg) * (4 - vg * (1 - gammaG^2)));
a3g = 1 / 4 * (1 / 2 - vg) * (1 + a4g);
a2g = 1 / gammaG * (a4g - a3g);
a1g = a2g + gammaG * a3g;

arg = [f Rc0 Rc1 Rg0 Rg1 omegaC1 omegaC2 omegaG1 omegaG2 ...
      vc vg a1c a2c a3c a4c a1g a2g a3g a4g nu0 k FLc FLg F0];
Rc = CPW_spec_func('Rc', arg); Rg = CPW_spec_func('Rg', arg);
R = Rc + Rg; % формула (5)


%% V. The Line Inductance L
LDC0 = CPW_spec_func('LDC', [w wg s t nu0]);
F1 = F0 + CPW_spec_func('K', k2) / CPW_spec_func('Ktouch', k2) - ...
          CPW_spec_func('K', k1) / CPW_spec_func('Ktouch', k1);
Lz1 = CPW_spec_func('LDC', [w 3 / 2 * w s t nu0]) - nu0 / 4 * 1 / F1; % формула (10)
omegaL0 = 4 / (nu0 * k * t * wg);
omegaL1 = 4 / (nu0 * k * t * w); % формулы (11)
omegaL2 = 18 / (nu0 * k * t^2);
Lz2 = sqrt(nu0 / (2 * omegaL2 * k)) * (FLc + FLg) / (4 * F0^2);
vz1 = log((LDC0 - Le8) / Lz1) / log(omegaL0 / omegaL1);
vz2 = log(Lz1 / Lz2) / log(omegaL1 / omegaL2);

% Appendix A.6:
nu1 = (w / wg)^4 * vz1 / (4 - vz1);
nu2 = (w / wg)^2 * vz1 / (4 - vz1);
nu3 = (2 * t / (9 * w))^3 * (vz2 - 1 / 2) / (vz2 + 5 / 2);
nu4 = (2 * t / (9 * w)) * (vz2 + 1 / 2) / (vz2 + 5 / 2);
a3L = ((vz2 - vz1) * (1 + nu1) * (1 - nu4) + 4 * nu2 + nu4 * (1 - 3 * nu1)) / ...
      ((vz1 - vz2) * (1 + nu1) * (1 - nu3) + 4 - nu3 * (1 - 3 * nu1));
a2L = 1 / (1 + nu1) * (a3L * (1 - nu3) - nu2 - nu4);
a4L = -9 / 2 * w / t * (nu4 + a3L * nu3);
a5L = (2 / 9 * t / w)^2 * a3L + a4L;
a1L = vz1 / (4 - vz1) + nu2 * a2L;
a0L = (1 - Le8 / LDC0) * (a1L + (w / wg)^2 * a2L);

arg = [f LDC0 Lz1 Lz2 Le8 omegaL0 omegaL1 omegaL2 ...
       vz1 vz2 a0L a1L a2L a3L a4L a5L nu0 k FLc FLg F0];
L = CPW_spec_func('L', arg); % формула (9)

return;
