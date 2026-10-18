function ok = test()

% Peaks at and beyond the range limits, and out-of-range counts

Range = [0 10];
nPoints = 11;  % points at 0,1,...,10
Pos = [-0.7 -0.3 0 10 10.3 10.7 NaN];

% 'round': half a point beyond each limit is still in range
[~,s,nOut] = makespec(Range,nPoints,Pos,[],'round');
ok(1) = s(1)==2 && s(end)==2 && sum(s)==4 && isequal(nOut,[1 1]);

% 'floor': nothing below the lower limit is included, up to one point above
% the upper limit is
[~,s,nOut] = makespec(Range,nPoints,Pos,[],'floor');
ok(2) = s(1)==1 && s(end)==3 && sum(s)==4 && isequal(nOut,[2 0]);

% 'split': nothing beyond either limit is included
[~,s,nOut] = makespec(Range,nPoints,Pos,[],'split');
ok(3) = s(1)==1 && s(end)==1 && sum(s)==2 && isequal(nOut,[2 2]);

% Warning only without the third output
warnState = warning('off','all');
lastwarn('');
[~,~] = makespec(Range,nPoints,20);
msg1 = lastwarn;
lastwarn('');
[~,~,~] = makespec(Range,nPoints,20);
msg2 = lastwarn;
warning(warnState);
ok(4) = ~isempty(msg1) && isempty(msg2);
