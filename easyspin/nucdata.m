% nucdata  Nuclear spin data 
%
%   Spin = nucdata(Isotopes)
%   [Spin,gn] = nucdata(Isotopes)
%   [Spin,gn,qm] = nucdata(Isotopes)
%   [Spin,gn,qm,abund] = nucdata(Isotopes)
%
%   Returns nuclear spin, gyromagnetic
%   ratio, quadrupole moment and natural
%   abundance of one or several nuclei.
%
%   Isotopes is a string specifying the
%   nucleus, e.g. '1H', '13C', '63Cu', '191Ir'.
%   If Isotopes is a comma-separated list of nuclei
%   like '14N,14N,14N,1H,63Cu', vectors are returned.
%   Attention, case-sensitive!

% Reads data from external file:
%   private/nucdata.txt

function varargout = nucdata(Isotopes)

if (nargin==0), help(mfilename); return; end

if iscell(Isotopes)
  if numel(Isotopes)==1
    Isotopes = Isotopes{1};
  end
end
%--------------------------------------------------------------
global IsotopeList
persistent Elements IsotopeMatrix
% Read isotope data file only once per MATLAB session!
if isempty(IsotopeList)
  %disp('Loading nuclear isotope database...');
  % Determine data file name
  esPath = fileparts(which(mfilename));
  DataFile = [esPath filesep 'private' filesep 'isotopedata.txt'];
  if ~exist(DataFile,'file')
    error('Could not open nuclear data file %s',DataFile);
  end
  
  % Load data file
  fh = fopen(DataFile);
  C = textscan(fh,'%f %f %s %s %s %f %f %f %f','commentstyle','%');
  fclose(fh);

  [IsotopeList.Protons,IsotopeList.Nucleons,IsotopeList.Radioactive,...
   IsotopeList.Element,IsotopeList.Name,IsotopeList.Spins,...
   IsotopeList.gns,IsotopeList.Abundances,IsotopeList.qms] = C{:};
  for k = 1:numel(IsotopeList.Spins)
    IsotopeList.Symbols{k} = sprintf('%d%s',IsotopeList.Nucleons(k),IsotopeList.Element{k});
  end

  IsotopeMatrix = zeros(max(IsotopeList.Protons),max(IsotopeList.Nucleons));
  for k=1:numel(IsotopeList.Spins)
    idxP = IsotopeList.Protons(k);
    idxN = IsotopeList.Nucleons(k);
    if (idxN>0)
      IsotopeMatrix(idxP,idxN) = 1;
    end
  end
  
  elmidx = IsotopeList.Protons(k);
  for k = 1:numel(IsotopeList.Spins)
    Elements(elmidx).Symbol = IsotopeList.Element{k};
    Elements(elmidx).Name = IsotopeList.Name{k};
  end
  
end
%--------------------------------------------------------------

if isempty(Isotopes)
  varargout = {[],[],[],[],{}};
  varargout = varargout(1:nargout);
  return
end

if ~ischar(Isotopes) && ~iscell(Isotopes)
  error('Argument must be a string or a cell array!');
end

if ~iscell(Isotopes)
  Nucs = nucstring2list(Isotopes);
else
  Nucs = Isotopes;
end

for k = numel(Nucs):-1:1
  idx = find(strcmp(Nucs{k},IsotopeList.Symbols));
  if isempty(idx)
    % check if element without #nucleons is given
    RequestedElement = Nucs{k};
    RequestedElement(~isletter(RequestedElement)) = [];
    iidx = find(strcmp(RequestedElement,IsotopeList.Element));
    % compile list of available isotopes
    if ~isempty(iidx)
      Message = [];
      for iIsotope = 1:numel(iidx)
        Message = [Message sprintf(IsotopeList.Symbols{iidx(iIsotope)}) ', '];
      end
      Message = [sprintf('Problem in isotopes list, entry %d (''%s''): Please specify one of ',k,Nucs{k}) Message(1:end-2) '.'];
    else
      Message = sprintf('Problem in isotopes list, entry %d (''%s''): Unknown element ''%s''.',k,Nucs{k},RequestedElement);
    end
    error(Message);
  end
end

% Find nuclei and get spin and gn
for k = numel(Nucs):-1:1
  idx = find(strcmp(Nucs{k},IsotopeList.Symbols));
  Spin(k) = IsotopeList.Spins(idx);
  gn(k) = IsotopeList.gns(idx);
  qm(k) = IsotopeList.qms(idx);
  Abund(k) = IsotopeList.Abundances(idx)/100; % percent -> fraction between 0 and 1
end

if (nargout==0)
  fprintf('Isotope  Spin  gn   qm/b   abundance\n');
  for k = numel(Nucs):-1:1
    fprintf('%-5s    %g   %g  %g  %g\n',Nucs{k},Spin(k),gn(k),qm(k),Abund(k));
  end
end

varargout = {Spin,gn,qm,Abund,Nucs};
varargout = varargout(1:nargout);

return
