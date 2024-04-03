% Compute a set orientations in a plane
%======================================================
% Computes N vectors in a plane normal to a given direction Z.

clear, clc, clf

% Define a plane by giving a direction normal to it
% a the number of vectors to compute
Z = [1; 1; 1];
N = 31;

% Construct two vectors in the plane perpendicular to Z
[phiZ,thetaZ] = vec2ang(Z);
y = ang2vec(phiZ,thetaZ+pi/2);
z = Z/norm(Z);
x = cross(y,z);

% Construct vectors in the plane as linear
% combinations of x and y.
chi = linspace(0,2*pi,N);  % angle in the plane
Vectors = x*cos(chi) + y*sin(chi);

% Plotting
for k = 1:N
  v = Vectors(:,k);
  line([0 v(1)],[0 v(2)],[0 v(3)]);
end
line([-1 1]*z(1),[-1 1]*z(2),[-1 1]*z(3),'Color','r','LineWidth',2);
axis equal
box on
