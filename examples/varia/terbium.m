% Resonances of a Tb4+ system, animated
%==========================================================================
% A little animation showing resonance fields for a Tb4+ system
% as the spectrometer frequency is changed.

clear, clf

% Define the Tb Hamiltonian. It contains high-order spin terms.
b = -2527.53/60;
B40 = b; B44 = 5*b;
c = -24.84/1260;
B60 = c; B64 = -21*c;
Tb.S = 7/2;
Tb.g = 2.0136;
Tb.B4 =     [B44 0 0 0 B40 0 0 0 0];
Tb.B6 = [0 0 B64 0 0 0 B60 0 0 0 0 0 0];

% Compute the energy level diagram using EasySpin's levels() function
B0max = 2500;
B0 = linspace(0,B0max,200); % mT
E = levels(Tb,'z',B0)/1e3;

% Plot the energy level diagram
h = plot(B0/1e3,E,'k');
set(h,'LineWidth',1.5);
xlabel('magnetic field (T)');
ylabel('energy (GHz)');

Exp.Range = [0 B0max];
Exp.CrystalOrientation = [0 0 0]; % crystal orientation in spectrometer

% Define a set of microwave frequencies (log scale)
mwFreq = 10.^linspace(1,2.7,601); % GHz

% Loop over all frequencies
for iFreq = 1:numel(mwFreq)
  Exp.mwFreq = mwFreq(iFreq); % GHz
  
  % Compute the resonance fields using EasySpins's eigfields()
  resB = eigfields(Tb,Exp);
  
  % Delete old transition lines
  delete(findobj('Tag','resonance'));

  % Locate resonant transitions and plot transition lines
  % Loop over all resonance fields
  for k = 1:numel(resB)
    % Determine all energies for a resonance field (linear interpolation)
    E_ = interp1(B0,E,resB(k));
    
    % Calculate the offsets of all transitions from mwFreq
    deltaE = E_ - E_.'; % transition energies
    Eoffset = deltaE - Exp.mwFreq;
    
    % Find the two levels u and v which are closest to resonant
    [~,tr] = min(abs(Eoffset(:)));
    [u,v] = ind2sub(size(Eoffset),tr);
    
    % Plot lines, tag them (for deletion in next iteration)
    line([resB(k) resB(k)]/1e3, E_([u v]),'Color','r','Tag','resonance');
  end
  
  % Redraw
  drawnow
end
