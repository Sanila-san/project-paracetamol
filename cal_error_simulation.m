clc; clear all;
tic
maxSweep = 100;

logdate = datestr(now(), 30);

Var_Gap = [];
Var_Width = [];

Var_S11l1 = [];
Var_S21l1 = [];
Var_S12l1 = [];
Var_S22l1 = [];
Var_Zl1 = [];

Var_S11l2 = [];
Var_S21l2 = [];
Var_S12l2 = [];
Var_S22l2 = [];
Var_Zl2 = [];

Var_S11p = [];
Var_S21p = [];
Var_S12p = [];
Var_S22p = [];


Var_Elf = [];
Var_Elr = [];
Var_e11 = [];
Var_e22 = [];
Var_e10e01 = [];
Var_e32e23 = [];
Var_e10e32 = [];
Var_e01e23 = [];
Var_e00 = [];
Var_e33 = [];



Z0 = 50; % base impedance
Error = 0; % LEGACY unused binary coefficient for backward compatibility

% As we have 3 lines the length of both lines is defined below
% Correct lengths are for respective bands in GHz:
% 450 um: 16-125
% 900 um: 7,9-63
% 1800 um: 3,9-31
% 3500 um: 2,1-16
% 5250 um: 1,4-11
% e6 = MHz, e9 = GHz
Length = 5250e-6;           % Physical length (L) of the IDEAL line, m 
Lengthl1 = 5250e-6;           % Physical length of Line (L), m 
Lengthl2 = 50e-6;           % Physical length of Thru (L2), m 

F = (1400e6 : 200e6 : 11e9)'; N = length(F); % frequency range and step

% here we generate S-parameters of the ideal transmission line


c = 299792458; % скорость света в вакууме, м/с

%[F_1 S11_1 S21_1 S12_1 S22_1] = VreadS2P('_Metas_23um_120GHz.s2p');
%[F_1 S11_1 S21_1 S12_1 S22_1] = VreadS2P('_GGB_CS-5_102_Line (VNA tools).s2p');
%[F_2 S11_2 S21_2 S12_2 S22_2] = VreadS2P('_S_CPW (no radiation loss).s2p');


% Description of Line
DielConstant = 9.59;   % Dielectric constant (Er)
LossTangent = 0.001;       % Loss tangent (tan d)
Conductivity = 35e6;       % Conductivity (k), S/m 

Width = 50e-6;  % Width of signal conductor (w), m 
Gap = 23e-6;  % Gap width (s), m
WidthGround = 280e-6;      % Width of ground conductor (wg), m
Thickness = 4.2e-6;          % Thickness (t), m
      

% Расчет исходных параметров линий:
for f = 1 : N
    [C, G, R, L] = CPW_EM_prop(F(f, 1), Thickness, WidthGround, Gap, Width, DielConstant, LossTangent, Conductivity);
    gamma1(f, 1) = sqrt((1i * 2 * pi * F(f, 1) * L + R) * (1i * 2 * pi * F(f, 1) * C + G));
    Z1(f, 1) = sqrt((1i * 2 * pi * F(f, 1) * L + R) / (1i * 2 * pi * F(f, 1) * C + G));    
end;

% Расчет S-параметров линий и их соединения:
S11ideal = zeros(N, 1); S12ideal = zeros(N, 1);
S21ideal = zeros(N, 1); S22ideal = zeros(N, 1);

S11l1r = zeros(N, 1); S12l1r = zeros(N, 1);
S21l1r = zeros(N, 1); S22l1r = zeros(N, 1);

S11l2r = zeros(N, 1); S12l2r = zeros(N, 1);
S21l2r = zeros(N, 1); S22l2r = zeros(N, 1);

for f = 1 : N
    S1 = conv_gz2s(gamma1(f, 1), Z1(f, 1), Lengthl1, Z0);
    %S2 = conv_gz2s(gamma2(f, 1), Z2(f, 1), l2, Z0);
    %S = conv_t2s(conv_s2t(S1) * conv_s2t(S2)); % две линии
    %S = conv_t2s(conv_s2t(S1) * conv_s2t(S2) * conv_s2t(S1)); % три линии
    S11l1r(f, 1) = S1(1, 1); S12l1r(f, 1) = S1(1, 2);
    S21l1r(f, 1) = S1(2, 1); S22l1r(f, 1) = S1(2, 2);
end;



for f = 1 : N
    S1 = conv_gz2s(gamma1(f, 1), Z1(f, 1), Lengthl1, Z0);
    %S2 = conv_gz2s(gamma2(f, 1), Z2(f, 1), l2, Z0);
    %S = conv_t2s(conv_s2t(S1) * conv_s2t(S2)); % две линии
    %S = conv_t2s(conv_s2t(S1) * conv_s2t(S2) * conv_s2t(S1)); % три линии
    S11l2r(f, 1) = S1(1, 1); S12l2r(f, 1) = S1(1, 2);
    S21l2r(f, 1) = S1(2, 1); S22l2r(f, 1) = S1(2, 2);
end;

%%% DEBUG DEBUG DEBUG
S12l1r = S12l1r .* 0 .+ 1+0i;
S21l1r = S21l1r .* 0 .+ 1+0i;
S11l1r = S11l1r .* 0+0i;
S22l1r = S22l1r .* 0+0i;


toc
% this piece sweeps entire calculation over initial phase of ideal DUT defined
% below
    avgS11p = [];
    avgS21p = [];


tic
idealDUTs11 = 0.33333;
idealDUTs21 = 1.0000;
for initphase = [0:30:330]



% this code creates S-parameters for ideal DUT, e.g. attenuator or thru.
% in our work 4 cases are considered:
% Sxx = 0.005 Sxy = 0.0003 for matched 40 dB attenuator
% Sxx = 0.005 Sxy = 1.0000 for matched 0 dB attenuator
% Sxx = 0.33333 Sxy = 0.0003 for mismatched 40 dB attenuator, VSWR = 2.0
% Sxx = 0.33333 Sxy = 1.0000 for mismatched 0 dB attenuator, VSWR = 2.0

S11ideal = complex_circle(idealDUTs11, initphase, 30, N);
S21ideal = complex_circle(idealDUTs21, initphase, 30, N);
S22ideal = S11ideal;
S12ideal = S21ideal;

   S11ideal = S11ideal';
   S21ideal = S21ideal';
   S12ideal = S12ideal';
   S22ideal = S22ideal';

% after S-parameters are generated we can simulate multiple calibrations

for i = 1:maxSweep


  % Description of Line
  DielConstant = 9.59 + randn(1, 1)/4 * 1e-2;   % Dielectric constant (Er)
  LossTangent = 0.001;       % Loss tangent (tan d)
  Conductivity = 35e6;       % Conductivity (k), S/m 
  
% Central conductor width and gaps are randomly variable within some limits
  
  Width = 50e-6 + randn(1, 1)/4 * 1e-6;  % Width of signal conductor (w), m 
  Gap = 23e-6 - randn(1, 1)/4 * 1e-6;  % Gap width (s), m
  WidthGround = 280e-6 + randn(1, 1)/4 * 1e-6;      % Width of ground conductor (wg), m
  Thickness = 4.2e-6;          % Thickness (t), m
  
  Var_Gap = [Var_Gap; Gap * 1e6];
  Var_Width = [Var_Width, Width * 1e6];  

  % Расчет исходных параметров линий:
  for f = 1 : N
      [C, G, R, L] = CPW_EM_prop(F(f, 1), Thickness, WidthGround, Gap, Width, DielConstant, LossTangent, Conductivity);
      gamma1(f, 1) = sqrt((1i * 2 * pi * F(f, 1) * L + R) * (1i * 2 * pi * F(f, 1) * C + G));
      Z1(f, 1) = sqrt((1i * 2 * pi * F(f, 1) * L + R) / (1i * 2 * pi * F(f, 1) * C + G));    
  end;
  
  % Расчет S-параметров 1 линии и их соединения:
  S11 = zeros(N, 1); S12 = zeros(N, 1);
  S21 = zeros(N, 1); S22 = zeros(N, 1);
  
  for f = 1 : N
      S1 = conv_gz2s(gamma1(f, 1), Z1(f, 1), Lengthl1, Z0);
      %S2 = conv_gz2s(gamma2(f, 1), Z2(f, 1), l2, Z0);
      %S = conv_t2s(conv_s2t(S1) * conv_s2t(S2)); % две линии
      %S = conv_t2s(conv_s2t(S1) * conv_s2t(S2) * conv_s2t(S1)); % три линии
      S11(f, 1) = S1(1, 1); S12(f, 1) = S1(1, 2);
      S21(f, 1) = S1(2, 1); S22(f, 1) = S1(2, 2);
  
  end;
      Var_S11l1 = [Var_S11l1, S11];      
      Var_S21l1 = [Var_S21l1, S21];
      Var_S12l1 = [Var_S12l1, S12];
      Var_S22l1 = [Var_S22l1, S22];
      
      % debugging output to figure out if this stuff works
      cstrcat('Width: ',num2str(Var_Width(i)),', ', 'Gap: ',num2str(Var_Gap(i)),', ', 'initial phase ', num2str(initphase), ', i = ', num2str(i))
##      i
      
      % Расчет исходных параметров линий:
  for f = 1 : N
      [C, G, R, L] = CPW_EM_prop(F(f, 1), Thickness, WidthGround, Gap, Width, DielConstant, LossTangent, Conductivity);
      gamma1(f, 1) = sqrt((1i * 2 * pi * F(f, 1) * L + R) * (1i * 2 * pi * F(f, 1) * C + G));
      Z1(f, 1) = sqrt((1i * 2 * pi * F(f, 1) * L + R) / (1i * 2 * pi * F(f, 1) * C + G));    
  end;
  
  % Расчет S-параметров 2 линии и их соединения:
  S11 = zeros(N, 1); S12 = zeros(N, 1);
  S21 = zeros(N, 1); S22 = zeros(N, 1);

%%% DEBUG DEBUG DEBUG Lengthl2 превратилась в Lengthl1 !!!111
  
  for f = 1 : N
      S1 = conv_gz2s(gamma1(f, 1), Z1(f, 1), Lengthl1, Z0);
      %S2 = conv_gz2s(gamma2(f, 1), Z2(f, 1), l2, Z0);
      %S = conv_t2s(conv_s2t(S1) * conv_s2t(S2)); % две линии
      %S = conv_t2s(conv_s2t(S1) * conv_s2t(S2) * conv_s2t(S1)); % три линии
      S11(f, 1) = S1(1, 1); S12(f, 1) = S1(1, 2);
      S21(f, 1) = S1(2, 1); S22(f, 1) = S1(2, 2);  
  end;
  

  
      Var_S11l2 = [Var_S11l2, S11];
      Var_S21l2 = [Var_S21l2, S21];
      Var_S12l2 = [Var_S12l2, S12];
      Var_S22l2 = [Var_S22l2, S22];
##      i
 
    
   S11l1 = Var_S11l1(1:size(Var_S11l1, 1),i);   
   S11l2 = Var_S11l2(1:size(Var_S11l2, 1),i);
   S21l1 = Var_S21l1(1:size(Var_S21l1, 1),i);
   S12l1 = Var_S12l1(1:size(Var_S12l1, 1),i);
   S12l2 = Var_S12l2(1:size(Var_S12l2, 1),i);
   S21l2 = Var_S21l2(1:size(Var_S21l2, 1),i);
   S22l1 = Var_S22l1(1:size(Var_S22l1, 1),i);
   S22l2 = Var_S22l2(1:size(Var_S22l2, 1),i);
   
   
   %%% DEBUG DEBUG DEBUG
S12l1 = S12l1 .* 0 .+ 1+0i;
S21l1 = S21l1 .* 0 .+ 1+0i;
S11l1 = S11l1 .* 0+0i;
S22l1 = S22l1 .* 0+0i;
   
   
   
   Gr = -1 + 0i;
   ed = 1 + 0i;
   S11r = Gr;
   S22r = Gr;
##   S11l1r = S11l1; 
##   S11l2r = S11l2; 
##   S12l1r = S12l1;
##   S12l2r = S12l2; 
##   S21l1r = S21l1; 
##   S21l2r = S21l2; 
##   S22l1r = S22l1;
##   S22l2r = S22l2;
  
  % Calculating error values by MI-3411-2013 formulations
   
   A1 = S11l1 .* S22l1 .- S21l1 .* S12l1;
   A2 = S11l2 .* S22l2 .- S21l2 .* S12l2;
   B1 = (S11r .- S11l1r) ./ (S11r .- S11l2r);
   B2 = (S22r .- S22l1r) ./ (S22r .- S22l2r);
   C1 = S21l1r .* S21l2 ./ (S21l2r .* S21l1);
   C2 = S12l1r .* S12l2 ./ (S12l2r .* S12l1);
   Elf = (B1 .* (Gr .- S11l2) .- C1 .* (Gr .+ S11l1)) ./ (C1 .* (A1 .- Gr .* S22l1) .- B1 .* (A2 .- Gr .* S22l2));
   Elr = (B2 .* (Gr .- S22l2) .- C2 .* (Gr .+ S22l1)) ./ (C2 .* (A1 .- Gr .* S11l1) .- B2 .* (A2 .- Gr .* S11l2));
   e11 = (C1 .* (ed .- S22l1 .* Elf) .- (ed .- S22l2 .* Elf)) ./ (C1 .* (S11l1 .- Elf .* A1) .- (S11l2 .- Elf .* A2));
   e22 = (C2 .* (ed .- S11l1 .* Elr) .- (ed .- S11l2 .* Elr)) ./ (C2 .* (S22l1 .- Elr .* A1) .- (S22l2 .- Elr .* A2));
   e10e01 = (S11r .- S11l1r) .* (ed .- e11 .* Gr) .* (ed .- e11 .* S11l1 .- S22l1 .* Elf .+ e11 .* Elf .* A1) ./ (Gr .- Gr .* S22l1 .* Elf .- S11l1 .+ Elf .* A1);
   e32e23 = (S22r .- S22l1r) .* (ed .- e22 .* Gr) .* (ed .- e22 .* S22l1 .- S11l1 .* Elr .+ e22 .* Elr .* A1) ./ (Gr .- Gr .* S11l1 .* Elr .- S22l1 .+ Elr .* A1);
   e10e32 = S21l1r .* (ed .- e11 .* S11l1 .- S22l1 .* Elf .+ e11 .* Elf .* A1) ./ S21l1;
   e01e23 = S12l1r .* (ed .- e22 .* S22l1 .- S11l1 .* Elr .+ e22 .* Elr .* A1) ./ S12l1;
   e00 = S11r .- (e10e01 .* Gr ./ (ed .- e11 .* Gr));
   e33 = S22r .- (e32e23 .* Gr ./ (ed .- e22 .* Gr));
   
   % simulation of measuring ideal line with "calibrated" VNA
      

   
   a = (S11ideal .- e00) ./ e10e01;
   b = S21ideal ./ e10e32;
   c = S12ideal ./ e01e23;
   d = (S22ideal .- e33) ./ e32e23;
   S11p = ((ed .+ e22 .* d) .* a .- (Elf .* b .* c)) ./ ((ed .+ e11 .* a) .* (ed .+ e22 .* d) .- (Elf .* Elr .* c .* b));
   S21p = ((ed .+ (e22 .- Elf) .* d) .* b) ./ ((ed .+ e11 .* a) .* (ed .+ e22 .* d) .- (Elf .* Elr .* c .* b));
   S12p = ((ed .+ (e11 .- Elr) .* a) .* c) ./ ((ed .+ e11 .* a) .* (ed .+ e22 .* d) .- (Elf .* Elr .* c .* b));
   S22p = ((ed .+ e11 .* a) .* d .- (Elr .* b .* c)) ./ ((ed .+ e11 .* a) .* (ed .+ e22 .* d) .- (Elf .* Elr .* c .* b));

   
% this commented code is mostly for test purposes
##   plot(F, abs(S11p));
##   title("S11p over frequency")
##   hold on;
   
   
   % this is test monitoring output for checking if this code works at all
##   i
##   Elf(1)
##   Elr(1)
##   e10e01(1)
##   e00(1)


   
   % Compositing error vectors over frequency
   
   Var_Elf = [Var_Elf, Elf];
   Var_Elr = [Var_Elr, Elr];
   Var_e11 = [Var_e11, e11];
   Var_e22 = [Var_e22, e22];
   Var_e10e01 = [Var_e10e01, e10e01];
   Var_e32e23 = [Var_e32e23, e32e23];
   Var_e10e32 = [Var_e10e32, e10e32];
   Var_e01e23 = [Var_e01e23, e01e23];
   Var_e00 = [Var_e00, e00];
   Var_e33 = [Var_e33, e33];

    Var_S11p = [Var_S11p, S11p];
    Var_S21p = [Var_S21p, S21p];
    Var_S12p = [Var_S12p, S12p];
    Var_S22p = [Var_S22p, S22p];

    %%%%
    % abs(Var_S11p(1, 1:size(Var_S11p, 2)))
    
%%%%%%
##    avgS11p = [avgS11p, mean(abs(S11p))];
##    avgS21p = [avgS21p, mean(abs(S21p))];
    
endfor
% this piece draws plots for eXX and eXXeXX parameters defined by MI 3411-2013 document
##   figure
##   subplot(1,2,1)   
##   plot(F, abs(Var_e10e01))
##   hold on;
##   plot(F, abs(Var_e32e23))
##   hold on;
##   plot(F, abs(Var_e10e32))
##   hold on;
##   plot(F, abs(Var_e01e23))
##   xlabel('Frequency, Hz');
##   ylabel('ee-parameters');
##   
## d  subplot(1,2,2)
##   plot(F, abs(Var_Elf))
##   hold on;
##   plot(F, abs(Var_Elr))
##   hold on;
##   plot(F, abs(Var_e00))
##   hold on;
##   plot(F, abs(Var_e33))
##   xlabel('Frequency, Hz');
##   ylabel('e-parameters');

##    Var_S11p = [Var_S11p, (max(Var_S11p')'- min(Var_S11p')')./2];
##    Var_S21p = [Var_S21p, (max(Var_S21p')'- min(Var_S21p')')./2];



deltaS11p = [];
deltaS21p = [];

deltaS11p = [deltaS11p, max(Var_S11p')', min(Var_S11p')',...
    (Var_S11p(:,size(Var_S11p,2)-1).-Var_S11p(:,size(Var_S11p,2)))];
deltaS21p = [deltaS21p, max(Var_S21p')', min(Var_S21p')',...
    (Var_S21p(:,size(Var_S21p,2)-1).-Var_S21p(:,size(Var_S21p,2)))];
##
##    deltaS11p = [deltaS11p, max(deltaS11p')', min(deltaS11p')',...
##    (deltaS11p(:,size(deltaS11p,2)-1).-deltaS11p(:,size(deltaS11p,2)))];
##deltaS21p = [deltaS21p, max(deltaS21p')', min(deltaS21p')',...
##    (deltaS21p(:,size(deltaS21p,2)-1).-deltaS21p(:,size(deltaS21p,2)))];

       % this writes raw calculated data into files by one file per initial phase angle

##    logdate = datestr(now(), 30);
##    filenameS11 = cstrcat('Initial_phase_',num2str(initphase), '_Line_', num2str(Lengthl1*1e6),'um_', 'Sxxi_', num2str(idealDUTs11),...
##     '_Sxyi_', num2str(idealDUTs21) ,"_Var_S11p_", logdate, ".txt");
##    filenameS21 = cstrcat('Initial_phase_',num2str(initphase),'_Line_', num2str(Lengthl1*1e6),'um_', 'Sxxi_', num2str(idealDUTs11),...
##     '_Sxyi_', num2str(idealDUTs21), "_Var_S21p_", logdate, ".txt");
    filenameAvgS11 = cstrcat('Avg_Initial_phase_',num2str(initphase), '_Line_', num2str(Lengthl1*1e6),'um_', 'Sxxi_', num2str(idealDUTs11),...
     '_Sxyi_', num2str(idealDUTs21) ,"_Var_S11p_", logdate, ".txt");
    filenameAvgS21 = cstrcat('Avg_Initial_phase_',num2str(initphase), '_Line_', num2str(Lengthl1*1e6),'um_', 'Sxxi_', num2str(idealDUTs11),...
     '_Sxyi_', num2str(idealDUTs21) ,"_Var_S21p_", logdate, ".txt"); 
##    dlmwrite(filenameS11, [F, abs(S11ideal),abs(Var_S11p)]);
##    dlmwrite(filenameS21, [F, abs(S21ideal),abs(Var_S21p)]);
    dlmwrite(filenameAvgS11, [abs(deltaS11p)]);
    dlmwrite(filenameAvgS21, [abs(deltaS21p)]);
    dlmwrite(cstrcat("AVGERR_", filenameAvgS11), [mean(abs(deltaS11p(:,3)))]);
    dlmwrite(cstrcat("AVGERR_", filenameAvgS21), [mean(abs(deltaS21p(:,3)))]);
    
    
endfor   
toc   

##    Var_S11p = [Var_S11p, (max(Var_S11p')'- min(Var_S11p')')./2];
##    Var_S21p = [Var_S21p, (max(Var_S21p')'- min(Var_S21p')')./2];

## this writes max eror data into the end of the matrix with calculated
## S-parameters. Max error data actually means absolute measurement error.

    Var_S11p = [Var_S11p, max(Var_S11p')', min(Var_S11p')',...
    (Var_S11p(:,size(Var_S11p,2)-1).-Var_S11p(:,size(Var_S11p,2)))];
    Var_S21p = [Var_S21p, max(Var_S21p')', min(Var_S21p')',...
    (Var_S21p(:,size(Var_S21p,2)-1).-Var_S21p(:,size(Var_S21p,2)))./2];
%    this writes raw calculated data into files


    filenameTS11 = cstrcat('TOTAL_Line_', num2str(Lengthl1*1e6),'um_', 'Sxxi_', num2str(idealDUTs11),...
     '_Sxyi_', num2str(idealDUTs21) ,"_Var_S11p_", logdate, ".txt");
    filenameTS21 = cstrcat('TOTAL_Line_', num2str(Lengthl1*1e6),'um_', 'Sxxi_', num2str(idealDUTs11),...
     '_Sxyi_', num2str(idealDUTs21), "_Var_S21p_", logdate, ".txt");
    dlmwrite(filenameTS11, [F, abs(S11ideal),abs(Var_S11p)]);
    dlmwrite(filenameTS21, [F, abs(S21ideal),abs(Var_S21p)]);
   
   

   
##   figure
##   subplot(1,2,1)   
##   semilogy(F/1e9, 20*log(abs(Var_S11p)), 'Color', [0 0.4 1], 'LineWidth', 2)
##   grid on;
##   hold on;
##   semilogy(F/1e9, 20*log(abs(S11ideal)), '+', 'Color', [1 1 0], 'LineWidth', 1)
##   xlabel('Frequency, GHz');
##   ylabel('Measured S21, S12');
##   title("Measured S11, S22, linear");
##   
##   subplot(1,2,2)
##   semilogy(F/1e9, 20*log(abs(Var_S21p)), 'Color', [0 0.4 1], 'LineWidth', 2)
##   grid on;
##   hold on;
##   semilogy(F/1e9, 20*log10(abs(S21ideal)), '+', 'Color', [1 1 0], 'LineWidth', 1)
##   xlabel('Frequency, GHz');
##   ylabel('Measured S21, S12');
##   title("Measured S21, S12, linear"); 

   
      figure
   subplot(1,2,1)   
   semilogy(F/1e9, abs(Var_S11p), 'Color', [0 0.4 1], 'LineWidth', 2)
   grid on;
   hold on;
   semilogy(F/1e9, abs(S11ideal),  'Color', [1 1 0], 'LineWidth', 1)
   xlabel('Frequency, GHz');
   ylabel('Measured S21, S12');
   title("Measured S11, S22, linear");
   
   subplot(1,2,2)
   semilogy(F/1e9, abs(Var_S21p), 'Color', [0 0.4 1], 'LineWidth', 2)
   grid on;
   hold on;
   semilogy(F/1e9, abs(S21ideal),  'Color', [1 1 0], 'LineWidth', 1)
   xlabel('Frequency, GHz');
   ylabel('Measured S21, S12');
   title("Measured S21, S12, linear"); 
   print -djpg "-S1350,400"  scattering_S11_0333_S21_1_linmag_wide_log.jpg




% let's plot polar constellations
   axisS11 = abs(S11ideal(1))*1.25;
   axisS21 = abs(S21ideal(1))*1.25;
   figure
   subplot(1,2,1)
   plot(Var_S11p, '.','Color', [0 0.6 1], 'LineWidth', 6)
   grid on;
   hold on;
   plot(S11ideal, 'o', 'Color', [1 0 0], 'LineWidth', 3)
   axis([-axisS11 axisS11 -axisS11 axisS11]);
   title("Measured S11, S22, polar");
   
   subplot(1,2,2)
   plot(Var_S21p, '+','Color', [0 0.6 1], 'LineWidth', 4)
   grid on;
   hold on;
   plot(S21ideal, 'o', 'Color', [1 0 0], 'LineWidth', 3)
   axis([-axisS21 axisS21 -axisS21 axisS21]);
   title("Measured S21, S12, polar");
   print -djpg "-S1350,600"  scattering_S11_0333_S21_1_linmag_polar.jpg
      
##   figure
##   plot(Var_Gap, '+','Color', [0 0.6 1], 'LineWidth', 2)
##   hold on;
##   plot(Var_Width, 'o','Color', [1 0.3 0], 'LineWidth', 2)
##   xlabel('Iteration');
##   ylabel('Gap (+) and Width (+)');
##   title('Gap and Width over iterations')
   

   

