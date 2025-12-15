clc;
format short g;

% Problem data.

h = 0.02;
b = 0.20;
L = 2.0;
E = 25e9;
alpha = 1.0e-5;

% Geometric properties.

A = b*h;
I = b*h^3/12;
r = sqrt(I/A)
lbd = L/r

% Critical temperature.

dTcr = pi^2/(alpha*lbd^2)

% Temperature x displacement curve.

dT = linspace(dTcr, 2.5*dTcr, 31);
v = 2*r*sqrt(dT/dTcr - 1);
dT'
v'