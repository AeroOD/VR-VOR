%Originally From:
%Rey-Martinez, J., Batuecas-Caletrio, A., Matiño, E., Trinidad-Ruiz, 
%G., Altuna, X., & Perez-Fernandez, N. (2018). Mathematical methods for 
%measuring the visually enhanced vestibulo?ocular reflex and preliminary
%results from healthy subjects and patient groups. 
%Frontiers in neurology, 9, 69. https://doi.org/10.3389/fneur.2018.00069
%Original source files: https://github.com/bendermh/VVOR
%
%***Requires Signal Processing Toolbox To Function ***

function [f,P1] = fourier(data)
Fs = 220;                     % Sampling frequency
[L,~] = size(data);           % Length of signal
X = data;                     % Data
Y = fft(X);
P2 = abs(Y/L);
P1 = P2(1:fix(L/2)+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
