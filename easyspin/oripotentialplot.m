% oripotentialplot     Plot orientational potential
%
%    oripotentialplot
%    oripotentialplot(Potential)
%
%   Plots the orientational potential U(alpha,beta,gamma) specified in Potential.
%   The potential is only plotted against alpha and beta; gamma is fixed at
%   zero.
%
%   Input:
%   Potential    array with L, M, K, and lambda values for the Wigner
%                expansion representation of the potential:
%                [L1 M1 K1 lambda; L2 M2 K2 lambda2; ...]
%
%   Example:
%     Potential = [2 2 1 1.4];
%     oripotentialplot(Potential);

function oripotentialplot(Potential)

% Internal settings
%-------------------------------------------------------------------------------
PlotStyle = 2; % 1: rectangular, 2: on a sphere
plotPopulation = false;
defaultPotential = [2 0 0 2];
%defaultPotential = [];

generatePotential = nargin==0;

if ~generatePotential
  if size(Potential,2)~=4
    error('Potential must be a Nx4 array');
  end
  
  % Remove zero terms
  Potential(Potential(:,4)==0,:) = [];
  
  Lp = Potential(:,1);
  Mp = Potential(:,2);
  Kp = Potential(:,3);
  lam = Potential(:,4);
end

% Generate a random potential if none is given
%-------------------------------------------------------------------------------
if generatePotential
  if ~isempty(defaultPotential)
    Lp = defaultPotential(1);
    Mp = defaultPotential(2);
    Kp = defaultPotential(3);
    lam = defaultPotential(4);
  else
    Lmax = 3;
    Leven = false;
    Mzero = false;
    realLambda = true;
    
    p = 0;
    if Leven, Lrange = 0:2:Lmax; else, Lrange = 0:Lmax; end
    for L_ = Lrange
      if Mzero, Mrange = 0; else, Mrange = -L_:L_; end
      for M_ = Mrange
        if M_<0, Kmin = 1; else, Kmin = 0; end
        for K_ = Kmin:L_
          p = p + 1;
          if M_==0 && K_==0
            lam(p) = 2*rand-1;
          else
            if realLambda
              lam(p) = 2*rand-1;
            else
              lam(p) = complex(2*rand-1,2*rand-1);
            end
          end
          Lp(p) = L_;
          Mp(p) = M_;
          Kp(p) = K_;
        end
      end
    end
  end
end

% Check potential parameters for compliance
%-------------------------------------------------------------------------------
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

kT = 2; % number is irrelevant
PotentialFunction = @(a,b,c) -kT*LMKsum(a,b,c,lam,Lp,Mp,Kp);
PopulationFunction =  @(a,b,c) exp(-PotentialFunction(a,b,c)/kT);

a = rand*2*pi;
b = linspace(0,pi,51);
c = linspace(0,2*pi,101);

normalizePopulation = false;
if normalizePopulation
  absTol = 1e-6;
  if Mzero
    Pint = integral3(@(a,b,c)PopulationFunction(a,b,c).*sin(b),0,2*pi,0,pi,0,2*pi,'AbsTol',absTol);
  else
    Pint = 2*pi*integral2(@(b,c)PopulationFunction(a,b,c).*sin(b),0,pi,0,2*pi,'AbsTol',absTol);
  end
else
  Pint = 1;
end


% Plotting
%--------------------------------------------------------------------------
clf

if PlotStyle==1
  [b_,a_] = meshgrid(b,a);
  if plotPopulations
    Pplot = 2*pi*PopulationFunction(a_,b_,c)/Pint;
  else
    Pplot = PotentialFunction(a_,b_,c)/kT;
  end
  
  pcolor(a*180/pi,b*180/pi,Pplot.');
  set(gca,'YDir','reverse');
  ylabel('beta / pi')
  xlabel('alpha / pi')
  shading flat
  colorbar
  title('Boltzmann population');
  set(gca,'XTick',0:45:360,'YTick',0:30:180);
  
elseif PlotStyle==2
  % Plot on a sphere, (alpha,beta)
  
  nKnots = max(30,max(Lp)*2);
  symmgroup = 'C1';
  Vectors = sphgrid(symmgroup,nKnots);
  [alpha,beta] = vec2ang(Vectors);
  if plotPopulation
    Pplot = PopulationFunction(alpha,beta,0);
    titlestr = 'Equilibrium population (alpha,beta,0)';
  else
    Pplot = PotentialFunction(alpha,beta,0)/kT;
    titlestr = 'Potential(alpha,beta,0) / kT';
  end
  vec = ang2vec(alpha,beta).';
  tri = sphtri(symmgroup,nKnots);
  
  % Display list of potential coefficients
  axes('Position',[0 0 1 1],'HitTest','off');
  axis off
  str = sprintf('Potential coefficients\n');
  for p = 1:numel(lam)
    str = sprintf('%s %d %d %d:  %f  %f\n',str,Lp(p),Mp(p),Kp(p),real(lam(p)),imag(lam(p)));
  end
  text(0.02,0.98,str,'VerticalAl','top','FontSize',9);
  
  axes
  h = trisurf(tri,vec(:,1),vec(:,2),vec(:,3),Pplot);
  h.FaceAlpha = 0.9;
  
  % Add xyz axes and labels
  v = 1.2; d = 0.1;
  line([0 0],[0 0],[-v v],'Color','k');
  line([0 0],[-v v],[0 0],'Color','k');
  line([-v v],[0 0],[0 0],'Color','k');
  text(v+d,0,0,'+x','VerticalAl','middle','Horizon','center');
  text(v+d,0,0,'+x','VerticalAl','middle','Horizon','center');
  text(0,v+d,0,'+y','VerticalAl','middle','Horizon','center');
  text(0,0,v+d,'+z','VerticalAl','middle','Horizon','center');
  text(-v-d,0,0,'-x','VerticalAl','middle','Horizon','center');
  text(0,-v-d,0,'-y','VerticalAl','middle','Horizon','center');
  text(0,0,-v-d,'-z','VerticalAl','middle','Horizon','center');
    
  view([115 22]);
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
    colorbar
    cl = get(gca,'CLim');
    cl(1) = 0;
    set(gca,'CLim',cl);
  else
    colormap('parula(256)');
    colorbar
    cl = get(gca,'CLim');
    set(gca,'CLim',cl);
  end
  
  
end

return

%===============================================================================
% General ordering potential function (real-valued)
%===============================================================================
function y = LMKsum(a,b,c,lam,Lp,Mp,Kp)
y = 0;
for p = 1:numel(lam)
  if Mp(p)==0 && Kp(p)==0
    y = y + lam(p)*wignerd([Lp(p) 0 0],a,b,c);
  else
    y = y + lam(p)*wignerd([Lp(p) Mp(p) Kp(p)],a,b,c) + ...
      (-1)^(Kp(p)-Mp(p))*lam(p)'*wignerd([Lp(p) -Mp(p) -Kp(p)],a,b,c);
  end
end

%{
% General ordering potential function (real-valued)
% (as implemented in startvec.m)
function u = LMKsum2(a,b,c)
u = 0;
for p = 1:numel(lambda)
  if lambda(p)==0, continue; end
  if Kp(p)==0 && Mp(p)==0
    u = u - wignerd([Lp(p) +Mp(p) +Kp(p)],b) * real(lambda(p));
  else
    u = u - 2*real(wignerd([Lp(p) +Mp(p) +Kp(p)],a,b,c) * lambda(p));
  end
end
%}
