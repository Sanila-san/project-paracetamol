%% Специальные функции для расчета параметров CPW:
% (type - тип функции, arg - входной аргумент)

function f = CPW_spec_func(type, arg)

switch type
    
    case 'F' % arg = [ttouch k1 pc0 pc1 pc2 s]
        tt = arg(1); k1 = arg(2); % расчеты по формуле (1):
        pc0 = arg(3); pc1 = arg(4); pc2 = arg(5); s = arg(6);
        if tt <= s / 2
            f = CPW_spec_func('K', k1) / CPW_spec_func('Ktouch', k1) + pc0 * ...
                (tt / s * (pc1 - log(2 * tt / s)) + (tt / s)^2 * ...
                 (1 - 3 / 2 * pc2 + pc2 * log(2 * tt / s)));
        else
            f = CPW_spec_func('K', k1) / CPW_spec_func('Ktouch', k1) + ...
                pc0 / 8 * (pc2 + 2) + tt / s;
        end;
        
    case 'K', f = ellipke(arg^2);
    case 'E', [~, f] = ellipke(arg^2); % M = k^2 
    case 'Ktouch', f = ellipke(sqrt(1 - arg^2)^2);
    case 'Etouch', [~, f] = ellipke(sqrt(1 - arg^2)^2);
        
    case 'FLc'
        tH = arg(1); s = arg(2); % расчеты по формуле (17):
        a = arg(3); b = arg(4); pc0 = arg(5); pc1 = arg(6); pc2 = arg(7);
        pc3 = arg(8); pc4 = arg(9); pc5 = arg(10); pc6 = arg(11);
        if tH <= s / 2
            f = pc0 / s * (1 / (a + b) * (pi * b + b * log(8 * pi * a / (a + b)) - ...
                (b - a) * log((b - a) / (b + a)) - b * log(2 * tH / s)) + ...
                tH / s * (pc1 * pc3 - pc2 - b / a * pc4 + pc5 + (pc2 - pc3 + b / a - 1 - pc5) * log(2 * tH / s)) + ...
                (tH / s)^2 * (pc3 * (1 - 3 / 2 * pc1) + 3 / 2 * pc1 - 2 * pc2 + 1 + 3 / 2 * b / a * pc4 - ...
                b / a * (b - a) / (b + a) + (2 * pc2 + pc1 * (pc3 - 1) - b / a * pc4) * log(2 * tH / s)));
        else
            f = 1 / (2 * s) + tH / s^2 + pc0 / s * (pi * b / (a + b) + 1 / 2 * pc6 + ...
                1 / 8 * (-pc1 + pc3 * (pc1 + 2) - b / a * pc4 - 2 * (a^2 + b^2) / (a * (a + b))));
        end;
        
    case 'FLg'
        tH = arg(1); s = arg(2); % расчеты по формуле (18):
        a = arg(3); b = arg(4); pc0 = arg(5); pc1 = arg(6); pc2 = arg(7);
        pc3 = arg(8); pc4 = arg(9); pc5 = arg(10); pc6 = arg(11);
        if tH <= s / 2
            f = pc0 / s * (1 / (a + b) * (pi * a + a * log(8 * pi * b / (b - a)) + ...
                b * log((b - a) / (b + a)) - a * log(2 * tH / s)) + ...
                tH / s * (a / b * pc1 * pc3 + (1 - a / b) * pc1 - pc2 - pc4 - pc5 + (-a / b * pc3 + pc2 + a / b - 1 + pc5) * log(2 * tH / s)) + ...
                (tH / s)^2 * (a / b * pc3 * (1 - 3 / 2 * pc1) + 3 / 2 * a / b * pc1 - 2 * pc2 + 2 - a / b + 3 / 2 * pc4 - ...
                (b - a) / (b + a) + (2 * pc2 + a / b * pc1 * (pc3 - 1) - pc4) * log(2 * tH / s)));
        else
            f = 1 / (2 * s) + tH / s^2 + pc0 / s * (pi * a / (a + b) - 1 / 2 * pc6 + ...
                1 / 8 * (-a / b * pc1 + a / b * pc3 * (pc1 + 2) - pc4 - 2 * (a^2 + b^2) / (b * (a + b))));
        end;
        
    case 'Rc'
        freq = arg(1); % расчеты по формуле (6)
        Rc0 = arg(2); Rc1 = arg(3); Rg0 = arg(4); Rg1 = arg(5);
        omegaC1 = arg(6); omegaC2 = arg(7); omegaG1 = arg(8); omegaG2 = arg(9);
        vc = arg(10); vg = arg(11);
        a1c = arg(12); a2c = arg(13); a3c = arg(14); a4c = arg(15);
        a1g = arg(16); a2g = arg(17); a3g = arg(18); a4g = arg(19);
        nu0 = arg(20); k = arg(21); FLc = arg(22); FLg = arg(23); F0 = arg(24);
        omega = 2 * pi * freq; % круговая частота
        if omega <= omegaC1
            f = Rc0 * (1 + a1c * (omega / omegaC1)^2);
        else if omega <= omegaC2
                f = Rc1 * (omega / omegaC2)^vc * (1 + a2c * (omegaC1 / omega)^2 + ...
                    a3c * (omega / omegaC2)^2);
            else
                f = sqrt(omega * nu0 / (2 * k)) * FLc / (4 * F0^2) * ...
                    (1 + a4c * (omegaC2 / omega)^2);
            end;
        end;
        
    case 'Rg'
        freq = arg(1); % расчеты по формуле (7)
        Rc0 = arg(2); Rc1 = arg(3); Rg0 = arg(4); Rg1 = arg(5);
        omegaC1 = arg(6); omegaC2 = arg(7); omegaG1 = arg(8); omegaG2 = arg(9);
        vc = arg(10); vg = arg(11);
        a1c = arg(12); a2c = arg(13); a3c = arg(14); a4c = arg(15);
        a1g = arg(16); a2g = arg(17); a3g = arg(18); a4g = arg(19);
        nu0 = arg(20); k = arg(21); FLc = arg(22); FLg = arg(23); F0 = arg(24);
        omega = 2 * pi * freq; % круговая частота
        if omega <= omegaG1
            f = Rg0 * (1 + a1g * (omega / omegaG1)^2);
        else if omega <= omegaG2
                f = Rg1 * (omega / omegaG2)^vg * (1 + a2g * (omegaG1 / omega)^2 + ...
                    a3g * (omega / omegaG2)^2);
            else
                f = sqrt(omega * nu0 / (2 * k)) * FLg / (4 * F0^2) * ...
                    (1 + a4g * (omegaG2 / omega)^2);
            end;
        end;
        
    case 'LDC' % Appendix A.4:
        w1 = arg(1); w2 = arg(2); s = arg(3); t = arg(4); nu0 = arg(5);
        f = nu0 / (8 * pi) * (4 / w1^2 * CPW_spec_func('gl', [w1 t]) + ...
            1 / w2^2 * (CPW_spec_func('gl', [w1 + 2 * s t]) + CPW_spec_func('gl', [w1 + 2 * w2 + 2 * s t]) + ...
            2 * CPW_spec_func('gl', [w2 t]) - 2 * CPW_spec_func('gl', [w1 + w2 + 2 * s t])) - ...
            4 / (w1 * w2) * (CPW_spec_func('gl', [w1 + w2 + s t]) - CPW_spec_func('gl', [w1 + s t]) + ...
            CPW_spec_func('gl', [s t]) - CPW_spec_func('gl', [w2 + s t]))); % (19)
    case 'gl'
        x = arg(1); t = arg(2);
        f = (1 / 12 * t^2 - 1 / 2 * x^2) * log(1 + (x / t)^2) + ...
            1 / 12 * x^4 / t^2 * log(1 + (t / x)^2) - ...
            2 / 3 * x * t * (atan(x / t) + (x / t)^2 * atan(t / x)); % (20)
        
    case 'L'
        freq = arg(1); % расчеты по формуле (9)
        LDC0 = arg(2); Lz1 = arg(3); Lz2 = arg(4); Le8 = arg(5);
        omegaL0 = arg(6); omegaL1 = arg(7); omegaL2 = arg(8);
        vz1 = arg(9); vz2 = arg(10); a0L = arg(11); a1L = arg(12);
        a2L = arg(13); a3L = arg(14); a4L = arg(15); a5L = arg(16);
        nu0 = arg(17); k = arg(18); FLc = arg(19); FLg = arg(20); F0 = arg(21);
        omega = 2 * pi * freq; % круговая частота
        if omega <= omegaL0
            f = LDC0 * (1 + a0L * (omega / omegaL0)^2);
        else if omega <= omegaL1
                f = Le8 + Lz1 * (omega / omegaL1)^vz1 * (1 + ...
                    a1L * (omegaL0 / omega)^2 + a2L * (omega / omegaL1)^2);
            else if omega <= omegaL2
                    f = Le8 + Lz2 * (omega / omegaL2)^vz2 * (1 + ...
                        a3L * (omegaL1 / omega)^2 + a4L * (omega / omegaL2));
                else % > omegaL2 (skin-effect range):
                    f = Le8 + sqrt(nu0 / (2 * omega * k)) * (FLc + FLg) / ...
                        (4 * F0^2) * (1 + a5L * (omegaL2 / omega));
                end;
            end;
        end;
        
end;

return;