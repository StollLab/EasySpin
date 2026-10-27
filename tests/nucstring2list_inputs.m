function ok = test()

% Parsing of nuclear spin lists: char and string input, whitespace, empty
% entries, custom isotope mixtures, and invalid input

f = @(varargin) runprivate('nucstring2list',varargin{:});

% Regular lists, char and string input
ok(1) = isequal(f('63Cu,14N,14N,1H'),{'63Cu','14N','14N','1H'});
ok(2) = isequal(f("63Cu,14N"),{'63Cu','14N'});
ok(3) = isequal(f(' 1H , 14N ,'),{'1H','14N'});
ok(4) = isequal(f(',1H,,14N'),{'1H','14N'});
ok(5) = isequal(f('(63,65)Cu,N,35Cl','m'),{'(63,65)Cu','N','35Cl'});
ok(6) = isequal(f({'1H','14N'}),{'1H','14N'});

% Empty input gives an empty 1x0 list
emptyInputs = {'', "", ' ', ',', []};
for k = 1:numel(emptyInputs)
  L = f(emptyInputs{k});
  ok(end+1) = iscell(L) && isequal(size(L),[1 0]); %#ok<AGROW>
end

% Invalid input
badInputs = {'1H;14N', "1H;14N", ["1H","14N"], 14, ['1H';'2H']};
for k = 1:numel(badInputs)
  try
    f(badInputs{k});
    ok(end+1) = false; %#ok<AGROW>
  catch
    ok(end+1) = true; %#ok<AGROW>
  end
end

% Round trip with nuclist2string
Nucs = {'63Cu','14N','1H'};
ok(end+1) = isequal(f(runprivate('nuclist2string',Nucs)),Nucs);
