% Compute a set orientations in a plane
%======================================================
% Computes N vectors in a plane normal to a given direction Z.
% This can be done more easily using EasySpin's rotplane
% function.

% Define a plane by giving a direction normal to it
% a the number of vectors to compute
Z = [1; 1; 1];
N = 9;

% Construct two vectors in the plane perpendicular to Z
[phiZ,thetaZ] = vec2ang(Z);
y = ang2vec(phiZ,thetaZ+pi/2);
x = cross(y,Z/norm(Z));

% Construct vectors in the plane as linear
% combinations of x and y.
chi = linspace(0,2*pi,N);  % angle in the plane
Vectors = x*cos(chi) + y*sin(chi);
