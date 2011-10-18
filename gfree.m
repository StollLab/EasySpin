% gfree  Free electron g factor 
%
%   g = gfree
%   [g,sigma] = gfree
%
%   Returns the g value of the free
%   electron. sigma is the standard
%   uncertainty. Value taken from
%   D.Hanneke et al, Phys.Rev.Lett. 100,
%   120801(2008).

function [g,sigma] = gfree

% D. Hanneke, S. Fogwell, G. Gabrielse
% New Measurement of the Electron Magnetic Moment and the Fine Structure Constant
% Phys.Rev.Lett. 100, 120801 (2008)
%g =     1.00115965218073*2;
%sigma = 0.00000000000028*2;
g =     1.00115965218073*2;
sigma = 0.00000000000028*2;

% B.Odom, D.Hanneke, B.D'Urso, G.Gabrielse,
% New Measurement of the Electron Magnetic Moment
% Using a One-Electron Quantum Cyclotron.
% Phys.Rev.Lett. 97, 030801 (2006).
% g/2 = 1.001 159 652 180 85 (76) [0.76 ppt]

% 2006 CODATA recommended value:
% Value -2.002 319 304 3622
% Standard uncertainty 0.000 000 000 0015
% Relative standard uncertainty 7.4 x 10-13
% Concise form -2.002 319 304 3622(15)
%g =     2.0023193043622;
%sigma = 0.0000000000015;
