% Compute reference values by linear interpolation.

% Input data.
P = 600;
R = 100;
E = 10^7;
I = 1/12;

TabWR = [
8.37e-05	7.99e-04
7.58e-01	1.15e-01
9.73e-01	1.43e-01
1.07e+00	1.56e-01
1.16e+00	1.69e-01
1.28e+00	1.86e-01
1.40e+00	1.99e-01
1.51e+00	2.14e-01
1.64e+00	2.29e-01
1.80e+00	2.48e-01
1.92e+00	2.61e-01
2.06e+00	2.77e-01
2.21e+00	2.92e-01
2.41e+00	3.12e-01
2.70e+00	3.36e-01
2.94e+00	3.55e-01
3.17e+00	3.73e-01
3.38e+00	3.88e-01
3.68e+00	4.07e-01
3.89e+00	4.18e-01
4.15e+00	4.31e-01
4.38e+00	4.42e-01
4.57e+00	4.51e-01
4.83e+00	4.62e-01
5.10e+00	4.72e-01
5.42e+00	4.83e-01
5.72e+00	4.93e-01
5.98e+00	5.01e-01
6.22e+00	5.08e-01
6.50e+00	5.16e-01
6.83e+00	5.23e-01
7.08e+00	5.30e-01
];

TabVR = [
8.37e-05	7.99e-04
1.81e-01	7.16e-04
5.11e-01	5.36e-03
8.50e-01	1.08e-02
1.15e+00	1.71e-02
1.42e+00	2.49e-02
1.70e+00	3.60e-02
1.97e+00	4.62e-02
2.33e+00	6.29e-02
2.67e+00	7.95e-02
3.02e+00	9.53e-02
3.29e+00	1.06e-01
3.58e+00	1.18e-01
3.85e+00	1.28e-01
4.16e+00	1.40e-01
4.44e+00	1.51e-01
4.75e+00	1.62e-01
5.01e+00	1.73e-01
5.30e+00	1.81e-01
5.61e+00	1.92e-01
5.98e+00	2.03e-01
6.28e+00	2.12e-01
6.61e+00	2.21e-01
6.96e+00	2.29e-01
7.13e+00	2.32e-01
];

TabUR = [
7.94e-03	7.95e-04
4.56e-01	5.91e-04
9.91e-01	5.94e-03
1.48e+00	1.45e-02
1.98e+00	2.63e-02
2.50e+00	3.96e-02
2.90e+00	5.06e-02
3.22e+00	5.92e-02
3.61e+00	6.95e-02
3.91e+00	7.65e-02
4.16e+00	8.20e-02
4.50e+00	8.98e-02
4.76e+00	9.45e-02
5.10e+00	1.02e-01
5.42e+00	1.07e-01
5.61e+00	1.10e-01
5.98e+00	1.16e-01
6.21e+00	1.21e-01
6.53e+00	1.25e-01
6.83e+00	1.29e-01
7.16e+00	1.32e-01
7.33e+00	1.32e-01
];

% Remove Bath normalization.

LfDesnorm  = (E*I)/(P*R^2);
ResDesnorm = R;

TabWR(:,1) *= LfDesnorm;
TabWR(:,2) *= ResDesnorm;


TabUR(:,1) *= LfDesnorm;
TabUR(:,2) *= ResDesnorm;

TabVR(:,1) *= LfDesnorm;
TabVR(:,2) *= ResDesnorm;

TabUR

% Testing somethigns...

Vals = 0:0.05:1;

Vals = [Vals' Vals' Vals' Vals'];
Vals(:,2:4) = 0.0;

DataSize = size(Vals)(1);


for i = 1:DataSize
  % Current lf.
  lf = Vals(i,1);

  % Get W results.
  wid  = 1;
  found = false;
  for j = 2 : size(TabWR)(1)
    if (lf > TabWR(j-1,1) && lf < TabWR(j,1))
      wid = j;
      found = true;
      break;
    end
  end

  if (found) 
    delta = TabWR(wid,1) - TabWR(wid-1,1);
    w     = 0;
    w    += ((lf - TabWR(wid-1,1))/delta) * TabWR(wid,2);
    w    += ((TabWR(wid,1) - lf)/delta) * TabWR(wid-1,2);
  else
    w = 0;
  end

  % Get V results.
  vid  = 1;
  found = false;
  for j = 2 : size(TabVR)(1)
    if (lf > TabVR(j-1,1) && lf < TabVR(j,1))
      vid = j;
      found = true;
      break;
    end
  end

  if (found) 
    delta = TabVR(vid,1) - TabVR(vid-1,1);
    v     = 0;
    v    += ((lf - TabVR(vid-1,1))/delta) * TabVR(vid,2);
    v    += ((TabVR(vid,1) - lf)/delta) * TabVR(vid-1,2);
  else
    v = 0;
  end

  Vals(i,3) = v;

  % Get U results.
  uid  = 1;
  found = false;
  for j = 2 : size(TabUR)(1)
    if (lf > TabUR(j-1,1) && lf < TabUR(j,1))
      uid = j;
      found = true;
      break;
    end
  end

  if (found) 
    delta = TabUR(uid,1) - TabUR(uid-1,1);
    u     = 0;
    u    += ((lf - TabUR(uid-1,1))/delta) * TabUR(uid,2);
    u    += ((TabUR(uid,1) - lf)/delta) * TabUR(uid-1,2);
  else
    u = 0;
  end
 
  % Update table.
  Vals(i,2) = w;
  Vals(i,3) = v;
  Vals(i,4) = u;
end

Vals(:,1)
Vals(:,2)
Vals(:,3)
Vals(:,4)





