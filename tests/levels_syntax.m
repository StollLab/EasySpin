function ok = test()

% Syntax tests

Sys.S = 1/2;
Sys.g = [2 3 4];

phi = 0.4564;
theta = 1.14343564;
chi = 0.16434345;
B = linspace(0,100,201);

% magnetic field vector, (phi,theta) orientation
E = levels(Sys,[phi theta],B);
E = levels(Sys,phi,theta,B);
[E,V] = levels(Sys,[phi theta],B);
[E,V] = levels(Sys,phi,theta,B);

% magnetic field vector, (phi,theta,chi) orientation
E = levels(Sys,[phi theta chi],B);
E = levels(Sys,phi,theta,B);
[E,V] = levels(Sys,[phi theta],B);
[E,V] = levels(Sys,phi,theta,B);

% magnetic field vector, string abbreviations for orientation
E = levels(Sys,'xy',B);
[E,V] = levels(Sys,'xy',B);

% magnetic field range [Bmin Bmax], string abbreviations for orientation
Brange = [30 100];
E = levels(Sys,'xy',Brange);
[E,V] = levels(Sys,'xy',Brange);

ok = true;
