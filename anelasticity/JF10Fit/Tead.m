function [Tad] = Tead(delz,Tp);
% Calculates the adiabatic temperature gradient for a given Tp
cp = 1350; alv = 2.9E-5; grav = 9.98;
Tad = alv * grav * Tp * delz/cp + Tp;
