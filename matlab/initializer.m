clear all;
clc;
clf;
close all;



tic
[opt_A,fval,exitflag,output_fmin,lambda,gradient,constraints] = optimizer();
toc

Num_blades = opt_A(1);
rotor_radius = opt_A(2)/10;
omega = opt_A(3)/2/pi*60*100;
Num_rotors = opt_A(4);
Blade_pitch = opt_A(5)*180/pi/10;
fuel_rate = opt_A(6)/10^3;
engine_size = opt_A(7)*100;
Flight_time = -fval*1e4/3600;

