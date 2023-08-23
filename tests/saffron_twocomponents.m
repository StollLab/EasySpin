function [ok,data] = test(opt,olddata)

% first test, make sure two components work
Exp.Sequence = 'MimsENDOR';
Exp.Field = 325;
Exp.tau = 0.1;
Exp.Range = larmorfrq('1H',Exp.Field) + [-1 1]*10;

Sys1.Nucs = '1H';
Sys1.A = [1 3];
Sys1.lwEndor = 0.1;
Sys2.Nucs = '1H';
Sys2.A = [7 9];
Sys2.lwEndor = 0.1;
Sys2.weight = 0.3;

Opt.Verbosity = 0;
[x,y] = saffron({Sys1,Sys2},Exp,Opt);
y = y/max(abs(y));

if (opt.Display)
  if ~isempty(olddata)
    subplot(3,1,[1 2]);
    plot(x,y,'r',x,olddata.y,'b');
    xlabel('time (µs)');
    legend('new','old');
    title(mfilename)
    subplot(3,1,3);
    plot(x,olddata.y-y);
    title('old - new');
  end
end

data.y = y;
data.x = x;

if ~isempty(olddata)
  ok(1) = areequal(olddata.y,y,1e-8,'abs');
  
  % second test - make sure we get an error if two components are combined
  % with user defined experiment
  Chirp90.Type = 'quartersin/linear';
  Chirp90.trise = 0.030;
  Chirp90.tp = 0.200;
  Chirp90.Flip = pi/2;
  Chirp90.Frequency = [-300 300]; % excitation band, GHz
  
  Exp2.Field = 1240; % run the experiment at Q band
  Exp2.mwFreq = 34.78;
  
  Exp2.Sequence = {Chirp90 0.5};
  
  Exp2.DetWindow = [-0.1 0.1];
  Exp2.DetPhase = 0;
  
  Opt.GridSize = 7;
  Opt.SimulationMode = 'thyme';
  
  try
    [~,~] = saffron({Sys1,Sys2},Exp2,Opt);
    ok(2) = false;
  catch
    ok(2) = true;
  end
else
  ok = [];
end



end