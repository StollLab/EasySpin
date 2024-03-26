% EPR resonance fields of a Tb(IV) system, animated
%===============================================================================
% An animation showing resonance fields for a Tb(IV) ion as the spectrometer
% frequency is changed.

% This is from Fig. 5.7 (p. 337) from Abragam/Bleaney
% Electron Paramagnetic Resonance of Transition Metal Ions, after
% J. M. Baker, J. R. Chadwick, G. Garton, J. P. Hurrell
% E.p.r and endor of Tb4+ in thoria
% Proceedings of the Royal Society of London. Series A, Mathematical and Physical Sciences
% Vol. 286, No. 1406 (Jul. 20, 1965), pp. 352-365

clear, clc, clf

% Define spin and g value for S=7/2 Tb(IV)
Tb.S = 7/2;
Tb.g = 2.0136;

% Define the 4th- and 6th-order spin Hamiltonian parameters
b = -2527.53/60;
B40 = b; B44 = 5*b;
c = -24.84/1260;
B60 = c; B64 = -21*c;
Tb.B4 =     [B44 0 0 0 B40 0 0 0 0];
Tb.B6 = [0 0 B64 0 0 0 B60 0 0 0 0 0 0];

% Compute the energy level diagram using EasySpin's levels() function
B0max = 2500;
B0 = linspace(0,B0max,200); % mT
E = levels(Tb,'z',B0)/1e3; % MHz -> GHz

% Plot the energy level diagram
h = plot(B0/1e3,E);
set(h,'LineWidth',1.5);
xlabel('magnetic field (T)');
ylabel('energy (GHz)');

Exp.Range = [0 B0max];
Exp.SampleFrame = [0 0 0]; % crystal orientation in spectrometer

% Define a set of microwave frequencies (log scale)
mwFreq = 10.^linspace(1,log10(500),201); % GHz

% Loop over all frequencies
for iFreq = 1:numel(mwFreq)
  Exp.mwFreq = mwFreq(iFreq); % GHz
  
  % Compute the resonance fields using EasySpins's resfields_eig()
  resB = resfields_eig(Tb,Exp);
  
  % Delete old transition lines
  delete(findobj('Tag','resonance'));

  % Locate resonant transitions and plot transition lines
  % Loop over all resonance fields
  for iB = 1:numel(resB)
    % Determine all energies for a resonance field (linear interpolation)
    E_ = interp1(B0,E,resB(iB));
    
    % Calculate the offsets of all transitions from mwFreq
    deltaE = E_ - E_.'; % transition energies
    Eoffset = deltaE - Exp.mwFreq;
    
    % Find the two levels u and v which are closest to resonant
    [~,tr] = min(abs(Eoffset(:)));
    [u,v] = ind2sub(size(Eoffset),tr);
    
    % Plot lines, tag them (for deletion in next iteration)
    line([resB(iB) resB(iB)]/1e3, E_([u v]),'Color','r','Tag','resonance','LineWidth',1);
  end
  
  title(sprintf('%0.3g GHz',mwFreq(iFreq)));
  
  % Redraw
  drawnow limitrate
end
