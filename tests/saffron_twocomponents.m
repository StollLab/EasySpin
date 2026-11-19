function [ok,data] = test(opt,olddata) %#ok<FNDEF>

% First test: make sure two components work
Exp.Sequence = 'MimsENDOR';
Exp.Field = 325;  % mT
Exp.tau = 0.1;  % µs
Exp.rfFreq = larmorfrq('1H',Exp.Field) + [-1 1]*10;  % MHz

Sys1.Nucs = '1H';
Sys1.A = [1 3];  % MHz
Sys1.lwEndor = 0.1;

Sys2.Nucs = '1H';
Sys2.A = [7 9];
Sys2.lwEndor = 0.1;
Sys2.weight = 0.3;

Opt.Verbosity = 0;
[x,y] = saffron({Sys1,Sys2},Exp,Opt);
y = y/max(abs(y));

if opt.Display
  if ~isempty(olddata)
    subplot(3,1,[1 2]);
    plot(x,y,x,olddata.y);
    legend('new','old');
    subplot(3,1,3);
    plot(x,olddata.y-y);
    legend('old - new');
    xlabel('time (µs)');
  end
end

data.y = y;
data.x = x;

if ~isempty(olddata)
  ok(1) = areequal(olddata.y,y,1e-5,'abs');
else
  ok = [];
end

end
