% Windowing functions

clear
clc
close all

%  Flat top windowing
L=1000;
n=0:L-1;
L1 = L-1;

a0=0.21557895;
a1=0.41663158;
a2=0.277263158;
a3=0.083578947;
a4=0.006947368;

w = a0 - a1*cos(2*pi*n/L1) + a2*cos(4*pi*n/L1) - a3*cos(6*pi*n/L1) + a4*cos(8*pi*n/L1);

plot(w, 'LineWidth', 2)
grid on
