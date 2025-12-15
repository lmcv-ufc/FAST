% ==========================================================================
% PlateBuckling.m - Evaluate the buckling load of a symple supported 
% square plate subjected to biaxial loading.
% ==========================================================================
% Created:  20-Mai-2018   Evandro Parente Junior
%
% Modified: 14-Oct-2024   Evandro Parente Junior
%           Moved all functions to a single file.
% ========================================================================== 

function PlateBuckling()

  clc;
  clear;
  format short g;
  
  % Geometry and loading.

  a = 2.0;         % Length dimension (m)
  b = 2.0;         % Width dimension (m)
  h = 1.25e-3;     % Thickness (m)
  k = 0;           % k = Ny/Nx.
  
  % Material properties (Steel).
  
  E1 = 200e9;
  E2 = E1;
  n12 = 0.3;
  G12 = E1/(2*(1 + n12));
  
  % Material properties (AS4/3506-1).

  E1  = 147.0e9;   % Elastic modulus parallel to the fibers (GPa)
  E2  =  10.3e9;   % Elastic modulus perpendicular to the fibers (GPa)
  n12 =  0.27;     % Poisson' ratio
  G12 =  7.0e9;    % Transverse modulus (GPa)
    
  % Composite layup.
  
  nl = 8;
  Lamina = [h,  45, E1, E2, n12, G12;  % Symmetric angle-ply
            h, -45, E1, E2, n12, G12;
            h,  45, E1, E2, n12, G12;
            h, -45, E1, E2, n12, G12;
            h, -45, E1, E2, n12, G12;
            h,  45, E1, E2, n12, G12;
            h, -45, E1, E2, n12, G12;
            h,  45, E1, E2, n12, G12];
        
  Lamina = [h,  90, E1, E2, n12, G12;  % Symmetric cross-ply
            h,   0, E1, E2, n12, G12;
            h,  90, E1, E2, n12, G12;
            h,   0, E1, E2, n12, G12;
            h,   0, E1, E2, n12, G12;
            h,  90, E1, E2, n12, G12;
            h,   0, E1, E2, n12, G12;
            h,  90, E1, E2, n12, G12];
               
  [Ncr, mcr, ncr] = CalcBuckLoad(a, b, Lamina, k);
  fprintf('Ncr = %e\n', Ncr);
  fprintf('Mode = (%d, %d)\n', mcr, ncr);
end

function [Ncr, mcr, ncr] = CalcBuckLoad(a, b, Lamina, k)

  % Compute ABD matrices.

  [A, B, D] = CalcABD(Lamina)
  
  % Compute the buckling load.
  
  Ncr = 0;
  mcr = 0;
  ncr = 0;
  for m = 1:20
    for n = 1:20
      N = pi^2*(D(1,1)*(m/a)^4 + 2*(D(1,2)+ 2*D(3,3))*(m/a)^2*(n/b)^2 + D(2,2)*(n/b)^4)/((m/a)^2 + k*(n/b)^2);
      if (m == 1 && n == 1) 
        Ncr = N;
        mcr = m;
        ncr = n;
      elseif (N < Ncr)
        Ncr = N;
        mcr = m;
        ncr = n;
      end
    end
  end
end

% ================================ CalcABD =================================

function [A, B, D, h] = CalcABD(Lamina)

  % Evalute the laminate thickness (h).

  n = length(Lamina);
  h = 0.0;                
  for i = 1:n 
    h = h + Lamina(i, 1);  
  end

  % Evaluate the {z} vector (bottom/top coordinates of each ply).

  z = zeros(n);
  z(1) = -h/2.0;  % Origin at laminate center
  for i = 1:n
    z(i+1) = z(i) + Lamina(i,1);
  end

  % Evaluate [A], [B], and [D] matrices.

  A = zeros(3, 3);
  B = zeros(3, 3);
  D = zeros(3, 3);
  for i = 1:n 
    Q = MatrixQ(Lamina(i,:));
    T = MatrixT(Lamina(i,:));
    Qb = T'*Q*T;
    A = A + Qb*(z(i+1) - z(i));
    B = B + Qb*(z(i+1)^2 - z(i)^2)/2.0;
    D = D + Qb*(z(i+1)^3 - z(i)^3)/3.0;
  end
end

% ================================ MatrizQ =================================

function Q = MatrixQ(lamina)
 
  E1  = lamina(3);
  E2  = lamina(4);
  n12 = lamina(5);
  G12 = lamina(6);
  n21 = n12*E2/E1;
  
  Q11 = E1/(1.0 - n12*n21);
  Q22 = E2/(1.0 - n12*n21);
  Q12 = n12*Q22;
  Q66 = G12;
  Q = [Q11, Q12, 0.0;
       Q12, Q22, 0.0;
       0.0, 0.0, Q66];
end

% ================================ MatrizT =================================

function T = MatrixT(lamina)
 
  ang = lamina(2);
  s  = sin(deg2rad(ang));
  c  = cos(deg2rad(ang));
  cs = c*s;
  s2 = s*s;
  c2 = c*c;  
  T = [   c2,   s2,    cs;
          s2,   c2,   -cs;
       -2*cs, 2*cs, c2-s2];
end

% ======================================================= End of file ======