## This function returns point on complex circle by given magnitude and angle

function [f] = complex_circle(mag, initial_phase, angle_step, steps)

if size(mag, 2) > 1
  mag = mag';
endif
  
if size(mag,2) > 1
  mag = mag(:,1);
endif

f = zeros(1,steps);

a = 0:steps-1;

f = mag*(cos(deg2rad(initial_phase + angle_step*a)) + 1i * sin(deg2rad(initial_phase + angle_step*a)));
##
##for i = 0:steps+1
##  f = [f, mag * cos(rad2deg(angle_step*i)) + 1i * mag * sin(rad2deg(angle_step*i))];
##  endfor

endfunction
##plot(complex_circle(0.5, 10, 0, 10), '+')
##hold on;
##plot(complex_circle(0.5, 20, 0, 10), '.')
##hold on;
##plot(complex_circle(0.5, 30, 0, 10), '*')
##hold on;
##plot(complex_circle(0.5, 40, 0, 10), '+')
##hold on;
##plot(complex_circle(0.5, 50, 0, 10), '+')
##hold on;
##plot(complex_circle(0.5, 0, 45, 8), '*')
