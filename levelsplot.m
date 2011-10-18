% levelsplot    plot energy level diagram 
%
%  levelsplot(Sys,Ori,B);
%  levelsplot(Sys,Ori,B,mwFreq);
%
%    Sys        spin system structure
%    Ori        2-element vector [phi theta]
%               Euler angles for the magnetic field
%               alternatively, either 'x', 'y' or 'z'
%    B          field range [Bmin Bmax] in mT
%    mwFreq     spectrometer frequency in GHz
%
%  If mwFreq is given, resonances are drawn. Red
%  lines indicate allowed transitions, gray lines
%  forbidden ones. If the lines are terminated with
%  dots, the relative transition intensity is larger
%  than 1%.
%
%  Example:
%    Sys = struct('S',7/2,'g',[2 2 2],'D',[1 1 -2]*5e3);
%    levelsplot(Sys,[0;pi/3],[0 6e3],95);

function levelsplot(Sys,varargin)

switch (nargin)
case 5
  [Ori,B,mwFreq,N] = deal(varargin{:});
case 4
  [Ori,B,mwFreq] = deal(varargin{:});
case 3
  [Ori,B] = deal(varargin{:});
  mwFreq = inf;
case 2
  error('Third input (magnetic field range [Bmin Bmax]) missing.');
case 0
  help(mfilename);
  return;
otherwise
  error('Wrong number of input arguments!');
end

if (nargout<0), error('Not enough output arguments.'); end
if (nargout>0), error('Too many output arguments.'); end

if isstruct(Ori)
  error('Second input argument (Ori) wrong: can''t be a structure.');
end
if isstruct(B)
  error('Third input argument (B) wrong: can''t be a structure.');
end

if ischar(Ori)
  str ='xyz';
  idx = find(Ori==str);
  if numel(idx)~=1
    error('Orientation must be a 2-vector or either ''x'', ''y'' or ''z''');
  end
  philist = [0 pi/2 0];
  thetalist = [pi/2 pi/2 0];
  Ori = [philist(idx); thetalist(idx)];
end

if numel(Ori)==2
  phi = Ori(1);
  theta = Ori(2);
else
  error('Ori must be a two-element vector.');
end

N = 200;
%N = 1001;

switch numel(B)
  case 1
    B = [0 B];
    Bvec = linspace(B(1),B(2),N);
  case 2
    Bvec = linspace(B(1),B(2),N);
  otherwise
    Bvec = B;
end

E = levels(Sys,[phi;theta],Bvec);

plot(Bvec,E/1e3,'b');

AllowedColor = [1 0 0];
ForbiddenColor = [1 1 1]*0.8;
Threshold = 0.01;

if ~isfield(Sys,'S')
  nElectrons = 1;
  Sys.S = 1/2;
else
  nElectrons = numel(Sys.S);
end

if ~isfield(Sys,'g');
  Sys.g = ones(1,nElectrons)*2;
end

if isfinite(mwFreq)
  Opt = struct('Threshold',0,'Freq2Field',0);
  Par = struct('mwFreq',mwFreq,'Range',B([1 end]));
  Par.Orientations = [phi;theta;0];
  [eB,tp] = resfields(Sys,Par,Opt);
  if ~isempty(eB)
    if nElectrons==1
      % one electron spin: normalize amplitudes
      n = ang2vec(phi,theta);
      g = norm(diag(Sys.g)*n);
      c = g*bmagn/planck/1e9;
      tpMax = floor(Sys.S+0.5)*ceil(Sys.S+0.5)/4*c^2;
      %tp(tp>tpMax) = tpMax; % guard agains round-off problems
      tpMax = max(tp);
    else
      tpMax = max(tp);
    end
    if (tpMax>0), tp = tp/tpMax; end
    [tp,i] = sort(tp);
    eB = eB(i);
    n = ang2vec(phi,theta);
    idx = tril(ones(hsdim(Sys)));
    [F,G] = sham(Sys,n);

    for ieB = 1:length(eB)
      %if tp(ieB)<=Threshold; continue; end
      H = F+G*eB(ieB);
      E = sort(eig(H));
      EE = E(:,ones(1,length(E)));
      dE = triu(EE'-EE)-mwFreq*1e3;
      dE(logical(idx)) = NaN;
      [i,j] = find(abs(dE)==min(abs(dE(:))));
      h = line(eB(ieB)*[1 1],E([i(1),j(1)])/1e3);
      Color = tp(ieB)*AllowedColor + (1-tp(ieB))*ForbiddenColor;
      set(h,'Color',Color);
      if tp(ieB)>Threshold,
        set(h,'Marker','.');
      end
    end
  else
    xl = xlim; yl = ylim;
    h = text(xl(1),yl(1),' no resonances in range!');
    set(h,'Color','r','VerticalAl','bottom');
    % no resonance fields
  end
end

axis tight;
xlabel('magnetic field [mT]');
ylabel('energy [GHz]');

return
