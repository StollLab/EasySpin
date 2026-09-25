function ok = test()

% Numeric mwFreq as fourth input draws transitions

Sys = struct('S',1,'D',[3000 500]);

hFig = figure('Visible','off');
levelsplot(Sys,'z',[0 600],9.5);
hTrans = findobj(hFig,'Tag','transition');
close(hFig);

ok = ~isempty(hTrans);
