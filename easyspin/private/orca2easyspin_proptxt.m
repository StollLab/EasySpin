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

for s = 1:numel(sections)
  idx = sections(s);
  sectionTitle = allLines{idx};
  switch sectionTitle
    case '$ SCF_Energy'
      L = allLines{idx+4};
      i = find(L==':');
      E = str2double(L(i+1:end));
    case '$ DFT_Energy'
      % nothing to read here
    case '$ Calculation_Info'
      L = allLines{idx+4};
      i = find(L==':');
      multiplicity = str2double(L(i+1:end));
      S = (multiplicity-1)/2;
    case '$ SCF_Electric_Propertoes'
      % nothing to read here
    case '$ EPRNMR_GTensor'
      g_raw = readmatrix(allLines(idx+(8:10)));
      [R,gpv] = eig(g_raw);
      g = diag(gpv).';
      gFrame = eulang(R);
    case '$ EPRNMR_ATensor'
      % not implemented yet
    case '$ EPRNMR_DTensor'
      % not implemented yet
  end
end

Sys = struct;
if ~isempty(S), Sys.S = S; end
if ~isempty(g), Sys.g = g; end
if ~isempty(gFrame), Sys.gFrame = gFrame; end
if ~isempty(E), Sys.Energy = E; end

end

function M = readmatrix(lines)
M(1,:) = sscanf(lines{1},'%g').';
M(2,:) = sscanf(lines{2},'%g').';
M(3,:) = sscanf(lines{3},'%g').';
M = M(:,2:end);
end
