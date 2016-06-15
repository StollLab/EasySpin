% Resonances of a Tb4+ system, animated
%==========================================================================
% A little animation showing resonance fields for a Tb4+ system
% as the spectrometer frequency is changed.

clear, clf

% Define the Th system. It contains high-order spin term!
b = -2527.53/60; c = -24.84/1260;
B40 = b; B44 = 5*b;
B60 = c; B64 = -21*c;
Tb.S = 7/2;
Tb.g = 2.0136;
Tb.B4 = [B44 0 0 0 B40];
Tb.B6 = [0 0 B64 0 0 0 B60];

% Computer the energy level diagram using EasySpin's level() function
N = 200;
maxB = 2500;
extB = linspace(0,maxB,N); % mT
E = levels(Tb,0,0,extB)/1e3;

% Plot the energy level diagram
h = plot(extB,E,'k');
set(h,'LineWidth',2);
xlabel('magnetic field [mT]');
ylabel('energy [GHz]');

Params.Range = [0 maxB];
Params.CrystalOrientation = [0 0 0]; % crystal orientation in spectrometer

% Avoid screen flickering
%set(gcf,'DoubleBuffer','on');

% Define a set of microwave frequencies
mwFreq = 10.^linspace(1,2.7,301);

% Loop over all frequencies
for iFreq = 1:length(mwFreq);
  Params.mwFreq = mwFreq(iFreq); % GHz
  
  % Compute the resonance fields using EasySpins's eigfields() function
  resB = eigfields(Tb,Params);
  
  % Convert eigenfields to indices
  idxB = 1 + (resB-extB(1))/(extB(2)-extB(1));
  
  % Delete old transition lines
  delete(findobj('Tag','resonance'));

  % Locate resonant transitions and plot transition lines
  NN = size(E,2);
  % Loop over all resonance fields
  for k= 1:length(idxB)
    % Determine all energies for a resonance field
    % (linear interpolation)
    d = rem(idxB(k),1);
    idx = fix(idxB(k));
    Ei = E(idx,:)*(1-d) + E(idx+1,:)*d;
    
    % Find the two levels u and v which are resonant
    EE = repmat(Ei,NN,1);
    EE = triu(EE-EE');
    EE = triu(EE-EE')-Params.mwFreq;
    [dum,tr] = min(abs(EE(:)));
    [u,v] = ind2sub([NN,NN],tr);
    
    % Plot lines, label them with tag (for deletion)
    line([resB(k) resB(k)], Ei([u v]),'Color','r','Tag','resonance');
  end
  
  % Redraw.
  drawnow
end

% end