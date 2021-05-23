%% Функция чтения *.s2p файлов

function [F, S11, S21, S12, S22] = VreadS2P(file)
k = 0; f = fopen(file);
while ~feof(f)
line = fgetl(f);
     switch line(1) 
        case '#'
            Symbol = textscan(line, '%s');                               
         case '!'          
         otherwise
              k=k+1;
              NUMBER=str2num(line); 
              F(k,1)=NUMBER(1);
           if (Symbol{1}{4}=='dB')|(Symbol{1}{4}=='DB')
              % dB - модуль в dB и фаза в градусах
              S11(k, 1) = 10^(NUMBER(2) / 20) * cos(NUMBER(3) * pi / 180) +...
                   1i * 10^(NUMBER(2) / 20) * sin(NUMBER(3) * pi / 180);
              S21(k, 1) = 10^(NUMBER(4) / 20) * cos(NUMBER(5) * pi / 180)+...
                   1i * 10^(NUMBER(4) / 20) * sin(NUMBER(5) * pi / 180);
              S12(k, 1) = 10^(NUMBER(6) / 20) * cos(NUMBER(7) * pi / 180)+...
                   1i * 10^(NUMBER(6) / 20) * sin(NUMBER(7) * pi / 180);
              S22(k, 1) = 10^(NUMBER(8) / 20) * cos(NUMBER(9) * pi / 180)+...
                   1i * 10^(NUMBER(8) / 20) * sin(NUMBER(9) * pi / 180);            
           end;            
           if (Symbol{1}{4}=='RI')
            % RI % действительная и мнимая части
            S11(k, 1) = NUMBER(2) + 1i * NUMBER(3);
            S21(k, 1) = NUMBER(4) + 1i * NUMBER(5);
            S12(k, 1) = NUMBER(6) + 1i * NUMBER(7);
            S22(k, 1) = NUMBER(8) + 1i * NUMBER(9);
           end;  
           if (Symbol{1}{4}=='MA')
           % MA % модуль в линейном масштабе и фаза в градусах               
            S11(k, 1) = NUMBER(2) * cos(NUMBER(3) * pi / 180)+...
                   1i * NUMBER(2) * sin(NUMBER(3) * pi / 180);
            S21(k, 1) = NUMBER(4) * cos(NUMBER(5) * pi / 180)+...
                   1i * NUMBER(4) * sin(NUMBER(5) * pi / 180);
            S12(k, 1) = NUMBER(6) * cos(NUMBER(7) * pi / 180)+...
                   1i * NUMBER(6) * sin(NUMBER(7) * pi / 180);
            S22(k, 1) = NUMBER(8) * cos(NUMBER(9) * pi / 180)+...
                   1i * NUMBER(8) * sin(NUMBER(9) * pi / 180);    
           end;
     end;    
end;
fclose(f);
return

