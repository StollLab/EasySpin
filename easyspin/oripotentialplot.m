% oripotentialplot     Plot orientational potential
%
%    oripotentialplot(Potential)
%
%   Plots the orientational potential U(alpha,beta,gamma) specified in Potential.
%
%   Input:
%     Potential    array with L, M, K, and lambda values for the Wigner
%                  expansion representation of the potential:
%                  [L1 M1 K1 lambda; L2 M2 K2 lambda2; ...]
%
%   Example:
%     Potential = [2 2 1 1.4]; % L=2, M=2, K=1, lambda=1.4
%     oripotentialplot(Potential);

function oripotentialplot(Potential)

if nargin==0
  help(mfilename);
  return
end

if nargin>1
  error('Only one input argument allowed.');
end

% Check inputs
%-------------------------------------------------------------------------------
if ~isnumeric(Potential) || size(Potential,2)~=4
  error('Potential must be a Nx4 array');
end

% Internal settings
%-------------------------------------------------------------------------------
PlotStyle = 'sphere'; % 'rectangle', 'sphere'
plotPopulation = false;
plotAlphaBeta = true;

% Unpack potential and check potential parameters for compliance
%-------------------------------------------------------------------------------
% Remove zero terms
Potential(Potential(:,4)==0,:) = [];

Lp = Potential(:,1);
Mp = Potential(:,2);
Kp = Potential(:,3);
lam = Potential(:,4);

if ~isempty(lam)
  if any(Lp<0)
    error('L values of potential coefficients must be nonnegative.');
  end
  if any(abs(Kp)>Lp)
    error('L and K values of potential coefficients do not satisfy -L<=K<=L.');
  end
  if any(abs(Mp)>Lp)
    error('L and M values of potential coefficients do not satisfy -L<=M<=L.');
  end
  if any(Kp<0)
    error('Only nonnegative values of K are allowed.');
  end
  Mzero = Mp==0;
  if any(Kp(Mzero)<0)
    error('For potential terms with M=0, K must be nonnegative.');
  end
  if ~isreal(lam(Kp==0 & Mp==0))
    error('Potential coefficients for M=K=0 must be real-valued.');
  end
end

% Define potential energy and population function
%-------------------------------------------------------------------------------
kT = 2; % number is irrelevant
PotentialFunction = @(a,b,c) -kT*LMKsum(a,b,c,lam,Lp,Mp,Kp);
PopulationFunction = @(a,b,c) exp(-PotentialFunction(a,b,c)/kT);

normalizePopulation = true;
if normalizePopulation
  absTol = 1e-6;
  if Mzero
    Z = integral3(@(a,b,c)PopulationFunction(a,b,c).*sin(b),0,2*pi,0,pi,0,2*pi,'AbsTol',absTol);
  else
    Z = 2*pi*integral2(@(b,c)PopulationFunction(0,b,c).*sin(b),0,pi,0,2*pi,'AbsTol',absTol);
  end
  PopulationFunction = @(a,b,c) PopulationFunction(a,b,c)/Z;
end

% Set up UI components
%--------------------------------------------------------------------------
clf

h = uicontrol('Style','pushbutton','Position',[5 5 50 20],'String','ab/bc',...
  'Callback',@abbc_toggle,'ToolTip','Toggle between alpha/beta and beta/gamma display.');
h = uicontrol('Style','pushbutton','Position',[60 5 50 20],'String','pot/pop',...
  'Callback',@potpop_toggle,'ToolTip','Toggle between potential and equilibrium population.');

% Plotting
%--------------------------------------------------------------------------

a = rand*2*pi;
b = linspace(0,pi,51);
c = linspace(0,2*pi,101);

% Display list of potential coefficients
hAx = axes;
axes('Position',[0 0 1 1],'HitTest','off');
axis off
str = sprintf('Potential coefficients (L,M,K):\n');
for p = 1:numel(lam)
  str = sprintf('%s (%d,%d,%d):  %+0.4f  %+0.4f\n',str,Lp(p),Mp(p),Kp(p),real(lam(p)),imag(lam(p)));
end
text(0.02,0.98,str,'VerticalAlignment','top','FontSize',9);
axes(hAx);

updatedisplay;

function updatedisplay

switch PlotStyle
  case 'rectangle'
  
    if plotAlphaBeta
      [b_,a_] = meshgrid(b,a);
      if plotPopulation
        Pplot = PopulationFunction(a_,b_,c);
      else
        Pplot = PotentialFunction(a_,b_,c)/kT;
      end
      pcolor(a*180/pi,b*180/pi,Pplot.');
      ylabel('\beta / pi')
      xlabel('\alpha / pi')
    else
      [b_,c_] = meshgrid(b,c);
      if plotPopulations
        Pplot = PopulationFunction(0,b_,c_);
      else
        Pplot = PotentialFunction(0,b_,c_)/kT;
      end
      pcolor(c*180/pi,b*180/pi,Pplot.');
      ylabel('\gamma / pi')
      xlabel('\alpha / pi')
    end
    set(gca,'YDir','reverse');
    shading flat
    c = colorbar;
    c.Label.String = 'energy (kB*T)';
    title('Boltzmann population');
    set(gca,'XTick',0:45:360,'YTick',0:30:180);

  case 'sphere'
    % Plot on a sphere, (alpha,beta)

    GridSize = max(40,max(Lp)*2);
    symmgroup = 'C1';
    [Grid,tri] = sphgrid(symmgroup,GridSize);
    Vectors = Grid.vecs;

    if plotAlphaBeta
      [alpha,beta] = vec2ang(Vectors);
      if plotPopulation
        Pplot = PopulationFunction(alpha,beta,0);
        titlestr = 'Equilibrium population density (\alpha,\beta,0)';
      else
        Pplot = PotentialFunction(alpha,beta,0)/kT;
        titlestr = 'Potential energy (\alpha,\beta,0)';
      end
      vec = ang2vec(alpha,beta).';
    else
      [gamma,beta] = vec2ang(Vectors);
      if plotPopulation
        Pplot = PopulationFunction(0,beta,gamma);
        titlestr = 'Equilibrium population density (0,\beta,\gamma)';
      else
        Pplot = PotentialFunction(0,beta,gamma)/kT;
        titlestr = 'Potential energy (0,\beta,\gamma)';
      end
      vec = ang2vec(gamma,beta).';
    end

    h = trisurf(tri.idx,vec(:,1),vec(:,2),vec(:,3),Pplot);
    h.FaceAlpha = 0.85;

    % Add xyz axes and labels
    v = 1.2; d = 0.1;
    line([0 0],[0 0],[-v v],'Color','k');
    line([0 0],[-v v],[0 0],'Color','k');
    line([-v v],[0 0],[0 0],'Color','k');
    phi = linspace(0,2*pi);
    x = cos(phi)*0.99;
    y = sin(phi)*0.99;
    z = zeros(size(phi));
    line(x,y,z,'Color','k');
    line(y,z,x,'Color','k');
    line(z,x,y,'Color','k');
    text(v+d,0,0,'+x','VerticalAl','middle','Horizon','center');
    text(v+d,0,0,'+x','VerticalAl','middle','Horizon','center');
    text(0,v+d,0,'+y','VerticalAl','middle','Horizon','center');
    text(0,0,v+d,'+z','VerticalAl','middle','Horizon','center');
    text(-v-d,0,0,'-x','VerticalAl','middle','Horizon','center');
    text(0,-v-d,0,'-y','VerticalAl','middle','Horizon','center');
    text(0,0,-v-d,'-z','VerticalAl','middle','Horizon','center');

    axis equal
    shading interp
    axis off
    axis tight
    set(gca,'PlotBoxAspectRatioMode','manual');
    set(gca,'DataAspectRatioMode','manual');
    title(titlestr);
    if plotPopulation
      m = colormap('gray(256)');
      colormap(flipud(m));
      c = colorbar;
      c.Label.String = 'equilibrium population density (rad^{-3})';
      cl = get(gca,'CLim');
      cl(1) = 0;
      set(gca,'CLim',cl);
    else
      colormap('parula(256)');
      c = colorbar;
      c.Label.String = 'potential energy ({\itk}_B{\itT})';
      cl = get(gca,'CLim');
      set(gca,'CLim',cl);
    end
  
  otherwise
    
    error('Unknown value ''%s'' for PlotStyle.',PlotStyle);

end

end
%===============================================================================
% Function that calculates sum of LMK terms (real-valued)
%===============================================================================
function y = LMKsum(a,b,c,lam,Lp,Mp,Kp)
y = 0;
for p = 1:numel(lam)
  if lam(p)==0, continue; end
  if Mp(p)==0 && Kp(p)==0
    y = y + lam(p)*wignerd([Lp(p) 0 0],a,b,c);
  else
    ph = (-1)^(Kp(p)-Mp(p));
    y = y + lam(p) *wignerd([Lp(p)  Mp(p)  Kp(p)],a,b,c) + ...
         ph*lam(p)'*wignerd([Lp(p) -Mp(p) -Kp(p)],a,b,c);
  end
end

end

function abbc_toggle(hObject,event)
plotAlphaBeta = ~plotAlphaBeta;
updatedisplay();
end

function potpop_toggle(hObject,event)
plotPopulation = ~plotPopulation;
updatedisplay();
end

end
