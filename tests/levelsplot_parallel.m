function ok = test()

% Parallel-mode transitions match resfields

Sys = struct('S',1,'D',[3000 500]);
Exp = struct('mwFreq',9.5,'mwMode','parallel');

hFig = figure('Visible','off');
levelsplot(Sys,'z',[0 600],Exp);
data = sortrows(vertcat(findobj(hFig,'Tag','transition').UserData),3);
close(hFig);

Exp.Range = [0 600];
Exp.SampleFrame = [0 0 0];
Opt = struct('Threshold',0,'Freq2Field',0);
[B,Int,~,Tr] = resfields(Sys,Exp,Opt);
Int = abs(Int)/max(abs(Int));
keep = Int>=1e-6;
ref = sortrows([Tr(keep,:) B(keep) Int(keep)],3);

ok = ~isempty(data) && size(data,1)==size(ref,1) && all(abs(data(:)-ref(:))<1e-6);
