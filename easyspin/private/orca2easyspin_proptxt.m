% orca2easyspin_maintxt   Read EPR properties from main ORCA output file

function Sys = orca2easyspin_proptxt(propfilename,HyperfineCutoff)

fh = fopen(propfilename);
allLines = textscan(fh,'%s','whitespace','','delimiter','\n');
allLines = allLines{1};
fclose(fh);

if ~contains(allLines{2},"!PROPERTIES!")
  error('This is not a valid text-based ORCA property file.');
end

sections = find(cellfun(@(L)L(1)=='$',allLines));

S = [];
g = [];
gFrame = [];
E = [];  % total energy, hartree
A = [];
AFrame = [];
Q = [];
QFrame = [];

for s = 1:numel(sections)
  idx = sections(s);
  sectionTitle = allLines{idx};
  switch sectionTitle
    case '$ SCF_Energy'
      E = readvalue(allLines{idx+4});

    case '$ DFT_Energy'
      % not relevant

    case '$ Calculation_Info'
      multiplicity = readvalue(allLines{idx+4});
      S = (multiplicity-1)/2;

    case '$ SCF_Electric_Properties'
      % not relevant

    case '$ EPRNMR_GTensor'
      g_raw = readmatrix(allLines(idx+(8:10)));
      g_vecs = readmatrix(allLines(idx+(13:15)));
      g_vals = readvec(allLines(idx+18));
      g = g_vals;
      gFrame = eulang(g_vecs);

    case '$ EPRNMR_ATensor'
      % colon is missing after "Number of stored nuclei"
      valuestr = regexp(allLines{idx+4},'\d+$','match','once');
      nNuclei = str2double(valuestr);
      A = zeros(nNuclei,3);
      AFrame = zeros(nNuclei,3);
      i = idx+7;
      for n = 1:nNuclei
        A_raw = readmatrix(allLines(i+(6:8)));
        A_vecs = readmatrix(allLines(i+(11:13)));
        A_vals = readvec(allLines(i+16));
        A(n,:) = A_vals;
        if ~all(A_vecs==00,'all')
          AFrame(n,:) = eulang(A_vecs);
        else
          AFrame(n,:) = [0 0 0];
        end
        i = i+18;
      end

    case '$ EPRNMR_QTensor'
      % colon is missing after "Number of stored nuclei"
      valuestr = regexp(allLines{idx+4},'\d+$','match','once');
      nNuclei = str2double(valuestr);
      i = idx+8;
      Q = zeros(nNuclei,3);
      QFrame = zeros(nNuclei,3);
      for n = 1:nNuclei
        Q_raw = readmatrix(allLines(i+(6:8)));
        Q_vecs = readmatrix(allLines(i+(11:13)));
        Q_vals = readvec(allLines(i+16));
        Q(n,:) = Q_vals;
        QFrame(n,:) = eulang(Q_vecs);
        i = i+18;
      end

    otherwise
      fprintf('Skipping section %s\n',sectionTitle);

  end
end

Sys = struct;
if ~isempty(S), Sys.S = S; end
if ~isempty(g), Sys.g = g; end
if ~isempty(gFrame), Sys.gFrame = gFrame; end
if ~isempty(A), Sys.A = A; end
if ~isempty(AFrame), Sys.AFrame = AFrame; end
if ~isempty(Q), Sys.Q = Q; end
if ~isempty(QFrame), Sys.QFrame = QFrame; end
if ~isempty(E), Sys.Energy = E; end

end

function M = readmatrix(lines)
M(1,:) = sscanf(lines{1},'%g').';
M(2,:) = sscanf(lines{2},'%g').';
M(3,:) = sscanf(lines{3},'%g').';
M = M(:,2:end);
end

function val = readvalue(L)
i = find(L==':');
val = str2double(L(i+1:end));
end

function vec = readvec(lines)
vec(1,:) = sscanf(lines{1},'%g').';
vec = vec(2:end);
end
