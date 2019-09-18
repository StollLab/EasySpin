%  md_readpsf  Read a PSF (protein structure file, in CHARMM or 
%                    X-PLOR format) to obtain nitroxide spin label 
%                    information.
%
%   psf = md_readpsf(FileName, ResName, LabelName, AtomNames);
%
%   Input:
%     FileName       character array
%                    Name of DCD file to be read.
%
%     SegName        character array
%                    Name of segment consisting of the protein and spin
%                    label.
%
%     ResName        character array
%                    Name of residue assigned to spin label side chain,
%                    e.g. "CYR1" is the default used by CHARMM-GUI.
%
%     LabelName      name of spin label, 'R1' or 'TOAC'
%
%     AtomNames      structure array
%                    Structure array containing the atom names used in the 
%                    PSF to refer to the following atoms in the nitroxide 
%                    molecule:
%
%                      R1:
%                                              ON (ONname)
%                                              |
%                                              NN (NNname)
%                                            /   \
%                                  (C1name) C1    C2 (C2name)
%                                           |     |
%                                 (C1Rname) C1R = C2R (C2Rname)
%                                           |
%                                 (C1Lname) C1L
%                                           |
%                                 (S1Lname) S1L
%                                          /
%                                (SGname) SG
%                                         |
%                                (CBname) CB
%                                         |
%                             (Nname) N - CA (CAname)
%
%                      TOAC:
%                                         ON (ONname)
%                                         |
%                                         NN (NNname)
%                                        /   \
%                             (CG1name) CG1  CG2 (CG2name)
%                                       |    |
%                             (CB1name) CB1  CB2 (CB2name)
%                                        \  /
%                             (Nname) N - CA (CAname)
%
%
%   Output:
%     psf            structure array
%                    NATOM: number of atoms to be included in simulation
%                    idx_Protein: indices of the protein and label atoms
%                    idx_ProteinCA: indices of protein's alpha carbon atoms
%                    idx_SpinLabel: indices of the spin label's atoms
%                    R1: 
%                      idx_ON: index corresponding to atom ONname
%                      idx_NN: " " NNname
%                      idx_C1: " " C1name
%                      idx_C2: " " C2name
%                      idx_C1R: " " C1Rname
%                      idx_C2R: " " C2Rname
%                      idx_C1L: " " C1Lname
%                      idx_S1L: " " S1Lname
%                      idx_SG: " " SGname
%                      idx_CB: " " CBname
%                      idx_CA: " " CAname
%                      idx_N: " " Nname
%                    TOAC: 
%                      idx_ON: index corresponding to atom ONname
%                      idx_NN: " " NNname
%                      idx_CG1: " " CG1name
%                      idx_CG2: " " CG2name
%                      idx_CB1: " " CB1name
%                      idx_CB2: " " CB2name
%                      idx_CA: " " CAname
%                      idx_N: " " Nname
%
% For more information on PSF file structure, see:
%  http://www.ks.uiuc.edu/Training/Tutorials/namd/namd-tutorial-unix-html/node23.html

function psf = md_readpsf(FileName, SegName, ResName, LabelName, AtomNames)

% Note:
% This function assumes the following order of sections in the PSF, which 
% is standard:
%   !NTITLE
%   !NATOM
%   !NBOND
%   ...

if ~ischar(FileName)
  error('FileName must be given as a character array.')
end

if ~strcmpi(FileName(end-3:end),'.psf')
  error('Please give the full filename including the ".psf" extension.')
end

if ~ischar(ResName)
  error('resName must be given as a character array.')
end


% check AtomNames

if any(~cellfun(@ischar,fieldnames(AtomNames)))
  error('All fields in AtomNames must be given as character arrays.')
end

switch LabelName
  case 'R1'
    if ~isfield(AtomNames,'ONname')||~isfield(AtomNames,'NNname')...
        ||~isfield(AtomNames,'C1name')||~isfield(AtomNames,'C2name')...
        ||~isfield(AtomNames,'C1Rname')||~isfield(AtomNames,'C2Rname')...
        ||~isfield(AtomNames,'C1Lname')||~isfield(AtomNames,'S1Lname')...
        ||~isfield(AtomNames,'SGname')||~isfield(AtomNames,'CBname')...
        ||~isfield(AtomNames,'CAname')||~isfield(AtomNames,'Nname')
      error('One or more required atom names is missing. Please check documentation.')
    end
  case 'TOAC'
    if ~isfield(AtomNames,'ONname')||~isfield(AtomNames,'NNname')...
        ||~isfield(AtomNames,'CGSname')||~isfield(AtomNames,'CGRname')...
        ||~isfield(AtomNames,'CBSname')||~isfield(AtomNames,'CBRname')...
        ||~isfield(AtomNames,'CAname')||~isfield(AtomNames,'Nname')
      error('One or more required atom names is missing. Please check documentation.')
    end
end


% open the file
FileID = fopen(FileName,'r');
if FileID<1, error('File "%s" could not be opened.', FileName); end
ensurefclose = onCleanup(@() fclose(FileID));

% proper PSF formatting dictates that the first line of the file should
% start with "PSF"
line = fgetl(FileID);
if ~ischar(line)||~strcmpi(line(1:3),'PSF')
  error('File "%s" does not have proper PSF format.', FileName)
end

% initialization
% -------------------------------------------------------------------------

% psf.isCHARMM = 0;  % CHARMM and X-PLOR formats only differ in how dihedral
%                    % force field terms are listed, which is irrelevant for
%                    % reading atom data

psf.NATOM = 0;  % this parameter should be checked against the trajectory 
                % file to ensure that the files are compatible
                
AtomFormat = '%10d %8s %8d %8s %8s %4s %14f%14f%8d %*f %*f';

reachedNTITLE = false;
reachedNATOM = false;

% Read the file
%--------------------------------------------------------------------------

while ~feof(FileID)
  line = strtrim(fgetl(FileID));
  
  if regexp(line, '.*\d.*!\w.*')
    
    line = strsplit(line);
    
    % each section should begin with a line in the following format:
    %    [nLines] ![section]...
    % where nLines gives the number of lines that the section occupies 
    % after this beginning line
    nLines = round(str2double(line{1}));
    section = line{2};
    
    % NTITLE section
    if ~isempty(strfind(section,'NTITLE')) %#ok
      reachedNTITLE = true;
    end
    
    % NATOM section
    if ~isempty(strfind(section,'NATOM'))  %#ok
      if ~reachedNTITLE
        error('Section ordering in "%s" is not standard. See documentation for proper formatting.', FileName)
      end
      reachedNATOM = true;
      psf.NATOM = nLines;
      
      % read entire NATOM section
      FileContents = textscan(FileID, AtomFormat, nLines);
      segmentNames = FileContents(2);
      segmentNames = segmentNames{1};
      residueNames = FileContents(4);
      residueNames = residueNames{1};
      atomNames = FileContents(5);
      atomNames = atomNames{1};
      mass = FileContents(8);
      psf.mass = mass{1};
      
      % pick first segment if not specified explicitly
      if isempty(SegName)
        SegName = segmentNames{1};
      end
      
      % filter for atoms belonging to the protein and spin label
      idx_ProteinLabel = strcmpi(segmentNames,SegName);
      if ~any(idx_ProteinLabel), error('SegName ''%s'' not found.', SegName), end
      segmentNames = segmentNames(idx_ProteinLabel);
      residueNames = residueNames(idx_ProteinLabel);
      atomNames = atomNames(idx_ProteinLabel);
      % the ProteinLabel indices will be read directly from the binary 
      % files, so convert to integer indices
      psf.idx_ProteinLabel = find(idx_ProteinLabel);
      
      % generate logical arrays to filter for alpha carbons and spin label atoms
      psf.idx_ProteinCA = strcmpi(segmentNames,SegName) & strcmpi(atomNames,'CA');
      idx_SpinLabel = strcmpi(segmentNames,SegName) & strncmpi(residueNames,ResName,4);
      psf.idx_SpinLabel = idx_SpinLabel;
      
      % locate spin label atoms
      findatomindex = @(name) strcmpi(atomNames(idx_SpinLabel),name);
      switch LabelName
        case 'R1'
          psf.idx_ON = findatomindex(AtomNames.ONname);
          psf.idx_NN = findatomindex(AtomNames.NNname);
          psf.idx_C1 = findatomindex(AtomNames.C1name);
          psf.idx_C2 = findatomindex(AtomNames.C2name);
          psf.idx_C1R = findatomindex(AtomNames.C1Rname);
          psf.idx_C2R = findatomindex(AtomNames.C2Rname);
          psf.idx_C1L = findatomindex(AtomNames.C1Lname);
          psf.idx_S1L = findatomindex(AtomNames.S1Lname);
          psf.idx_SG = findatomindex(AtomNames.SGname);
          psf.idx_CB = findatomindex(AtomNames.CBname);
          psf.idx_CA = findatomindex(AtomNames.CAname);
          psf.idx_N = findatomindex(AtomNames.Nname);
        case 'TOAC'
          psf.idx_ON = findatomindex(AtomNames.ONname);
          psf.idx_NN = findatomindex(AtomNames.NNname);
          psf.idx_CGS = findatomindex(AtomNames.CGSname);
          psf.idx_CGR = findatomindex(AtomNames.CGRname);
          psf.idx_CBS = findatomindex(AtomNames.CBSname);
          psf.idx_CBR = findatomindex(AtomNames.CBRname);
          psf.idx_CA = findatomindex(AtomNames.CAname);
          psf.idx_N = findatomindex(AtomNames.Nname);
      end

    end
    
    if ~isempty(strfind(section,'NBOND'))  %#ok
      if ~reachedNTITLE||~reachedNATOM
        error('Section ordering in "%s" is not standard. See documentation for proper formatting', FileName)
      end
      % we are done
      break
    end
  
  end
end

end
