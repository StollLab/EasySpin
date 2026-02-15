function [ok,data] = test(opt,olddata)

% Symmetry test
  
Sy = struct('S',1/2,'Nucs','14N','g',[2.3 2.3 2],...
  'A',[44,49,54],'Q',-1*[-1,-1,2],'lwEndor',0);
Ex = struct('Range',[17 31],'Field',326.5,'nPoints',2^10,'mwFreq',9.7);
Si = struct('Threshold',1e-5,'GridSize',[20 5],'Intensity','on',...
  'Enhancement','off','Verbosity',opt.Verbosity,...
  'GridSymmetry','Ci','GridFrame',[0 0 0],'separate','transitions');

Sy.AFrame = pi/180*[0 -30 0];
Sy.QFrame = pi/180*[0 45 0];
[a,b,info] = salt(Sy,Ex,Si);

if opt.Display
  if ~isempty(olddata)
    subplot(3,1,[1 2]);
    plot(a,olddata.b,'k',a,b,'r');
    title('Low symmetry test (Ci, tilted A and Q)');
    subplot(3,1,3);
    plot(a,olddata.b-b);
    title('old - new');
  else
    plot(a,b);
    title('Low symmetry test (Ci, tilted A and Q)');
  end
  xlabel('frequency (MHz)');
end

data.b = b;

if isempty(olddata)
  ok = [];
else
  ok = areequal(olddata.b/max(olddata.b(:)),b/max(b(:)),1e-6,'abs');
end
