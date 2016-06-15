% Hamiltonian symmetry using symm()
%==========================================================================
% A tutorial on Hamiltonian symmetry as determined
% by EasySpin's function symm().

% USAGE: Execute the code lines one after
% the other to see all intermediate results.

clear

% Take a single isotropic spin-1/2 system
Sys.S = 1/2;
Sys.g = 2.002;

% It's symmetry is O3, meaning that the
% spectrum of the associated spin Hamiltonian
% is independent of the magnetic field
% orientation
disp(['Point group:  ', symm(Sys)])

% For an axial spin system, the point group
% is Dinfh. The spectrum depends only on
% theta (angle between northpole and magnetic
% field direction.

Sys.g = [2 2 2.1];
disp(['Point group:  ', symm(Sys)])

% Note that in the previous example the rotation
% axis was along z. What if we set it along x?

Sys.g = [2.1 2 2];

% We get the same symmetry!
[PG,R] = symm(Sys);
disp(['Point group:  ',PG])

% The second output argument describes the frame
% in which this symmetry holds.

R

% The z axis of the symmtry frame is 

R(:,3)

% the original x axis.

% For an orthorhombic spin system we get
% D2h symmetry.

Sys.g = [1.9 2 2.1];
symm(Sys)

% This means that the spectrum
% depends on both polar angles of the magnetic
% field, but there is a symmetry relationship
% between all 8 octants. In powder simulations,
% only one octant

% Now we look at some cases with
% an additional nucleus. First we take an isotropic
% g and an axial hyperfine coupling.

Sys.g = 2;
Sys.Nucs = '1H';
Sys.A = [-1 -1 2]*10;

symm(Sys)

% The symmetry is determined by A alone. If we
% tilt the orientation of the A tensor, the
% symmetry doesn't change, but the symmetry frame
% does.

Sys.AFrame = [12 23 34]*pi/180;
[PG,R] = symm(Sys)

% The rotation axis is along

R(:,3)

% which corresponds to the angles [in radians]
[phi,theta] = vec2ang(R(:,3))

% When both g and A are axial and collinear,
% the symmetry is axial

Sys.AFrame = [0 0 0];
Sys.g = [2 2.1];

symm(Sys)

% If A is tilted along the third Euler angle
% (i.e. rotated around the z axis of the tensor) it remains axial.

Sys.AFrame = [0 0 56]*pi/180;
symm(Sys)

% But if it is tilted along the second Euler
% angle, i.e. away from the z axis, the
% symmetry gets lower.

Sys.AFrame = [0 29 0]*pi/180;
symm(Sys)

% C2h means that there is a 2-fold rotation
% axis and a mirror plane parallel to it. for
% a powder spectrum simulation, 2 octants have
% to be integrated.

% If g is axial and A orthorhombic, and only their z and x (y) axes
% do not coincide, the symmetry is the same.

Sys.A = [-0.5 -1.5 2]*10;
Sys.g = [2 2.1];
Sys.AFrame = [0 47 0]*pi/180;

symm(Sys)

% But if all three axes do no coincide, the symmetry
% is Ci, containing only the inversion centre. A full
% hemisphere (4 octants) has to be integrated for
% a powder spectrum simulation

Sys.AFrame = [0 7 13]*pi/180;
symm(Sys)


% If g and A are both orthorhombic, the symmetry is
% only orthorhombic if their frames coincide. Otherwise
% it is either C2h or Ci.

Sys.g = [1.9 2 2.1];
Sys.A = [-0.5 -1.5 2]*5 + 6;

Sys.AFrame = [0 0 0];
fprintf('Collinear frames: %s\n',symm(Sys));

Sys.AFrame = [0 0 12];
fprintf('Third Euler angle non-zero: %s\n',symm(Sys));

Sys.AFrame = [0 61 0];
fprintf('Second Euler angle non-zero: %s\n',symm(Sys));

Sys.AFrame = [0 -73 14];
fprintf('Last two Euler angle non-zero: %s\n',symm(Sys));

Sys.AFrame = [-12 -52 5];
fprintf('Three Euler angle non-zero: %s\n',symm(Sys));

% END *****************************
