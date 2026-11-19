function ok = test() %#ok<FNDEF>

% Pathway parser tests

% Correct parsing for pathways containing 'c' (+1 or -1)
p.Flip = pi/2;
E = struct('Sequence',{{p 0.1 p 0.1 p 0.1}});
ok(1) = isequal(orders(E,'+c-'),[1 -1 -1; 1 1 -1]);
ok(2) = isequal(orders(E,'c0-'),[-1 0 -1; 1 0 -1]);

% An RF pulse keeps the coherence order: 'c' selects the one that is possible
rf = struct('Channel','rf','tp',0.02,'Flip',pi);
Er = struct('Sequence',{{p 0.1 rf 0.1 p 0.1}},'rfFreq',14);
ok(3) = isequal(orders(Er,'+c-'),[1 1 -1]);

% Error for nonexisting pathway and for invalid characters
ok(4) = throwsError(Er,'0c-','No coherence transfer pathway');
ok(5) = throwsError(E,'+x-','can only contain');

% Both engines: same signal as with the explicit pathways
Sys.g = [2.00 2.01 2.03];
Exp = E;
Exp.mwFreq = 9.5;  % GHz
Exp.Field = 338;  % mT
Exp.MolFrame = [0.2 0.8 0.1];
Exp.DetWindow = 0;
Exp.nPoints = 4;
Exp.Dim1 = {'d3' 0.02};
tol = 1e-10;
for method = {'subspace','fullspace'}
  Opt.Method = method{1};
  Exp.Pathways = '+c-';
  [~,yc] = saffron(Sys,Exp,Opt);
  Exp.Pathways = ['++-';'+--'];
  [~,yExplicit] = saffron(Sys,Exp,Opt);
  ok(end+1) = areequal(yc,yExplicit,tol,'rel') && any(yc(:)~=0);  %#ok<AGROW>
end

end

function co = orders(Exp,pathways)
Exp.Pathways = pathways;
seq = runprivate('sf_parsesequence',Exp,false);
co = unique(seq.co,'rows');  % alpha and beta populations both have order 0
end

function thrown = throwsError(Exp,pathways,text)
thrown = false;
try
  orders(Exp,pathways);
catch ME
  thrown = contains(ME.message,text);
end
end
