function ok = test()

% Level population markers for spin-polarized and thermal triplet

Sys = struct('S',1,'D',[1500 300]);
Ori = 'xy';
B = [0 500];
mwFreq = 9.5;

hFig = figure('Visible','off');

% Spin-polarized triplet: compare with populations computed independently
pZF = [0.2 0.3 0.5];
Sys.initState = {pZF,'zerofield'};
levelsplot(Sys,Ori,B,mwFreq);
data = vertcat(findobj(hFig,'Tag','population').UserData);  % [level B pop]

[V0,E0] = eig(ham(Sys,[0 0 0]),'vector');
[~,idx] = sort(real(E0));
V0 = V0(:,idx);
rho = V0*diag(pZF)*V0';
[H0,muz] = ham(Sys,[1 1 0]/sqrt(2));
popRef = zeros(size(data,1),1);
for k = 1:size(data,1)
  [V,E] = eig(H0-muz*data(k,2),'vector');
  [~,idx] = sort(real(E));
  v = V(:,idx(data(k,1)));
  popRef(k) = real(v'*rho*v);
end
ok(1) = ~isempty(data) && areequal(data(:,3),popRef,1e-8,'abs');

% Thermal populations: not shown by default
Sys = rmfield(Sys,'initState');
Exp = struct('mwFreq',mwFreq,'Temperature',5);
levelsplot(Sys,Ori,B,Exp);
ok(2) = isempty(findobj(hFig,'Tag','population'));

% Thermal populations if requested: lower level more populated than upper level
levelsplot(Sys,Ori,B,Exp,struct('Populations',true));
data = vertcat(findobj(hFig,'Tag','population').UserData);
data = sortrows(data,[2 1]);
ok(3) = ~isempty(data) && all(data(1:2:end,3)>data(2:2:end,3));

% No populations in the high-temperature limit
levelsplot(Sys,Ori,B,mwFreq);
ok(4) = isempty(findobj(hFig,'Tag','population'));

close(hFig);
