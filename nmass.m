% nmass   neutron mass 
%
%   mn = nmass
%   [mn,sigma] = nmass
%
%   Returns the mass of the neutron
%   in kg. sigma is the standard
%   uncertainty (2010 CODATA).

function [mn,sigma] = nmass

% 2010 CODATA value

%Value 	 1.674 927 351 x 10-27 kg
% Standard uncertainty 	 0.000 000 074 x 10-27 kg
%  Relative standard uncertainty 	  4.4 x 10-8
% Concise form 	 1.674 927 351(74) x 10-27 kg

mn =    1.674927351e-27;
sigma = 0.000000074e-27;
