function ok = test()

% Exp with only mwFreq gives the same transitions as numeric mwFreq

Sys = struct('S',1,'D',[3000 500]);
Exp = struct('mwFreq',9.5);

hFig = figure('Visible','off');
levelsplot(Sys,'xy',[0 600],9.5);
data1 = sortrows(vertcat(findobj(hFig,'Tag','transition').UserData));
levelsplot(Sys,'xy',[0 600],Exp);
data2 = sortrows(vertcat(findobj(hFig,'Tag','transition').UserData));
close(hFig);

ok = ~isempty(data1) && isequal(data1,data2);
