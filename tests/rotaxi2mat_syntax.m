function ok = test()

% Syntax checks
phi = 6*pi/7;

% Explicit rotation axis vector
v = [0.4 0.7 -1.2];
R = rotaxi2mat(v,phi);
R = rotaxi2mat(v.',phi);

% Shorthands for rotation axis
R = rotaxi2mat('x',phi);
R = rotaxi2mat('y',phi);
R = rotaxi2mat('z',phi);
R = rotaxi2mat('xy',phi);
R = rotaxi2mat('xz',phi);
R = rotaxi2mat('yz',phi);
R = rotaxi2mat('xyz',phi);

ok = true;
