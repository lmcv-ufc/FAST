clear;
clc;

% Geometria.

A = 1
L = 1

% Material.

nu = 0.3
E  = 200e9
fy = 400e6
K = E/(3*(1 - 2*nu))
G = E/(2*(1 + nu))
C = K + 4/3*G
Lambda = K - 2/3*G

% Carga

sig = 1e6;
P = sig*A;

% Carregamento elastico.

epsy = fy/(2*G)
sigy = C*epsy
uy = epsy*L;


