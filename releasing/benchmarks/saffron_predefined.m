% saffron_predefined  Benchmark of saffron's predefined experiments
%
%   R = saffron_predefined
%   R = saffron_predefined(R0)
%
% Runs the predefined pulse EPR experiments (subspace engine) for a set of
% spin systems and variants, and returns run times and signals in R. With a
% baseline R0 (from an earlier run, e.g. on main), it prints the slowdown
% and the largest relative signal difference for each case.
%
% Run from the EasySpin root folder with easyspin/ on the path:
%   addpath('easyspin'); addpath('releasing/benchmarks');
%   R0 = saffron_predefined; save R0 R0       % on the baseline code
%   R = saffron_predefined(R0);               % on the new code

function R = saffron_predefined(R0)

if nargin<1, R0 = []; end

warning('off','all');

% New interface: Opt.Method, macros with Exp.tp + Opt.IdealPulses instead
% of Exp.ExciteWidth
newInterface = exist(fullfile(fileparts(which('saffron')),'private','sf_engine_subspace.m'),'file')==2;

% Spin systems
Sys = struct;
Sys(1).name = '14N(Q)+1H';
Sys(1).Sys = struct('Nucs','14N,1H','A',[1.2 0.3 0.4; 3 4 9],'Q',[0.8 0],'g',[2.0 2.01 2.02]);
Sys(2).name = '14N+1H';
Sys(2).Sys = struct('Nucs','14N,1H','A',[1.2 0.3 0.4; 3 4 9],'Q',[0 0],'g',[2.0 2.01 2.02]);
S3 = struct('Nucs','14N','A',[1.2 0.3 0.4],'Q',0,'g',[2.0 2.01 2.02]);
H = {[3 4 9],[1 2 5],[-2 -2 4],[0.5 0.8 3],[2 2 -4]};
for k = 1:3, S3 = nucspinadd(S3,'1H',H{k}); end
Sys(3).name = '14N+3x1H';
Sys(3).Sys = S3;
S5 = S3;
for k = 4:5, S5 = nucspinadd(S5,'1H',H{k}); end
Sys(4).name = '14N+5x1H';
Sys(4).Sys = S5;

tau = 0.136; dt = 0.016;
base.Field = 350;
base.tau = tau;
base.dt = dt;

Exps = {};
Exps{end+1} = {'2pESEEM',  setfield(rmfield(base,'tau'),'nPoints',512)}; %#ok<*SFLD>
Exps{end+1} = {'3pESEEM',  setfield(base,'nPoints',512)};
Exps{end+1} = {'4pESEEM',  setfield(base,'nPoints',512)};
Exps{end+1} = {'HYSCORE',  setfield(base,'nPoints',[128 128])};

cases = {};
for iS = 1:numel(Sys)
  for iE = 1:numel(Exps)
    name = Exps{iE}{1};
    if iS==4 && ~strcmp(name,'3pESEEM'), continue; end  % large system: 3p only
    E = Exps{iE}{2};
    E.Sequence = name;
    O = struct('GridSize',31);
    cases{end+1} = {sprintf('%-10s %-8s plain',Sys(iS).name,name),Sys(iS).Sys,E,O}; %#ok<*AGROW>
    if iS==2
      % orientation selection
      Eo = E; Eo.mwFreq = 9.8; Eo.ExciteWidth = 100;
      cases{end+1} = {sprintf('%-10s %-8s orisel',Sys(iS).name,name),Sys(iS).Sys,Eo,O};
      % relaxation
      Sr = Sys(iS).Sys; Sr.T1 = 5; Sr.T2 = 1;
      cases{end+1} = {sprintf('%-10s %-8s relax',Sys(iS).name,name),Sr,E,O};
      % large grid
      Og = O; Og.GridSize = 91;
      if ~strcmp(name,'HYSCORE')
        cases{end+1} = {sprintf('%-10s %-8s grid91',Sys(iS).name,name),Sys(iS).Sys,E,Og};
      end
    end
  end
end
% Mims ENDOR (new code uses a different RF model: timing only)
Se = struct('Nucs','1H','A',[2 3 9],'g',[2.0 2.01 2.02],'lwEndor',0.1);
Ee = struct('Field',350,'Sequence','MimsENDOR','tau',0.2,'Range',[10 20],'nPoints',1001);
cases{end+1} = {sprintf('%-10s %-8s plain','1H','MimsENDOR'),Se,Ee,struct('GridSize',31)};

nCases = numel(cases);
R.names = cell(1,nCases);
R.time = zeros(1,nCases);
R.y = cell(1,nCases);
for c = 1:nCases
  [name,S,E,O] = cases{c}{:};
  if newInterface
    [E,O] = translate(E,O);
  end
  saffron(S,E,O);  % warm up
  [~,y] = saffron(S,E,O);
  nReps = 5;
  tic
  for k = 1:nReps
    [~,y] = saffron(S,E,O);
  end
  R.time(c) = toc/nReps;
  R.names{c} = name;
  R.y{c} = y;
  line_ = sprintf('%-32s %8.3f s',name,R.time(c));
  if ~isempty(R0)
    k0 = find(strcmp(R0.names,name));
    if ~isempty(k0)
      y0 = R0.y{k0};
      if isequal(size(y0),size(y))
        d = max(abs(y(:)-y0(:)))/max(abs(y0(:)));
      else
        d = NaN;
      end
      line_ = sprintf('%s   base %8.3f s   ratio %5.2f   diff %8.1e',line_,R0.time(k0),R.time(c)/R0.time(k0),d);
    end
  end
  fprintf('%s\n',line_);
end

end

% Translate old-style inputs to the new interface
function [E,O] = translate(E,O)
if isfield(E,'ExciteWidth')
  E.tp = 1/E.ExciteWidth;  % ideal pulses with an excitation window 1/tp
  E = rmfield(E,'ExciteWidth');
  O.IdealPulses = true;
end
end
