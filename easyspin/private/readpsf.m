%  readpsf  Read a PSF (protein structure file, in CHARMM or X-PLOR format) 
%           to obtain nitroxide spin label information.
%
%   psf = readpsf(FileName, ResName, AtomNames);
%
%   Input:
%     FileName       character array
%                    Name of DCD file to be read.
%
%     ResName        character array
%                    Name of residue assigned to spin label side chain,
%                    e.g. "CYR1" is the default used by CHARMM-GUI.
%
%     AtomNames      structure array
%                    Structure array containing the atom names used in the 
%                    PSF to refer to the following atoms in the nitroxide 
%                    molecule:
%
%                                   O (OName)
%                                   |
%                                   N (NName)
%                                  / \
%                        (C1Name) C   C (C2Name)
%                                 |   |
%                                ... ...
%
%
%   Output:
%     psf            structure array
%                    NATOM: number of atoms to be included in simulation
%                    idx_SpinLabel: indices of the spin label's atoms
%                    idx_O: index corresponding to atom OName
%                    idx_N: index corresponding to atom NName
%                    idx_C1: index corresponding to atom C1Name
%                    idx_C2: index corresponding to atom C2Name
%
% For more information on PSF file structure, see:
%  http://www.ks.uiuc.edu/Training/Tutorials/namd/namd-tutorial-unix-html/node23.html

function psf = readpsf(FileName, ResName, AtomNames)

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

if ~endsWith(FileName,'.psf')
  error('Please give the full filename including the ".psf" extension.')
end

if ~ischar(ResName)
  error('resName must be given as a character array.')
end


% check AtomNames

if any(~cellfun(@ischar,fieldnames(AtomNames)))
  error('All fields in AtomNames must be given as character arrays.')
end

if ~isfield(AtomNames,'OName')||~isfield(AtomNames,'NName')...
    ||~isfield(AtomNames,'C1Name')||~isfield(AtomNames,'C2Name')
  error('One or more required atom names is missing. Please check documentation.')
end

% open the file
FileID = fopen(FileName, 'r');
if FileID<1, error('File "%s" could not be opened.', FileName); end
ensurefclose = onCleanup(@() fclose(FileID));

% proper PSF formatting dictates that the first line of the file should
% start with "PSF"
line = fgetl(FileID);
if ~ischar(line)||~startsWith(line,'PSF')
  error('File "%s" does not have proper PSF format. See documentation for details.', FileName)
end

% initialization
% -------------------------------------------------------------------------

% psf.isCHARMM = 0;  % CHARMM and X-PLOR formats only differ in how dihedral
%                    % force field terms are listed, which is irrelevant for
%                    % reading atom data

psf.NATOM = 0;  % this parameter should be checked against the trajectory 
                % file to ensure that the files are compatible
                
AtomFormat = '%10d %8s %8d %8s %8s %4s %14f%14f%8d %*f %*f';

reachedNTITLE = 0;
reachedNATOM = 0;

% read the file
% -------------------------------------------------------------------------

while ~feof(FileID)
  line = strtrim(fgetl(FileID));
  
  if regexp(line, '.*\d.*!\w.*')
    
    line = strsplit(line);
    
    % each section should begin with a line in the following format:
    %    [nLines] ![section]...
    % where nLines gives the number of lines that the section occupies 
    % after this beginning line
    nLines = str2num(line{1});
    section = line{2};
    
    if contains(section,'NTITLE')
      % skip this section
      reachedNTITLE = 1;
    end
    
    if contains(section,'NATOM')
      if ~reachedNTITLE
        error('Section ordering in "%s" is not standard. See documentation for proper formatting.', FileName)
      end
      reachedNATOM = 1;
      psf.NATOM = nLines;
      
      % read NATOM section
      data = textscan(FileID, AtomFormat, nLines);
      residue_names = data(4);
      residue_names = residue_names{1};
      atom_names = data(5);
      atom_names = atom_names{1};
      
      % obtain logical arrays of spin label's and nitroxide coordinate 
      % system's atoms
      idx_SpinLabel = strncmpi(residue_names,ResName,4);
      idx_O = strcmpi(atom_names(idx_SpinLabel),AtomNames.OName);
      idx_N = strcmpi(atom_names(idx_SpinLabel),AtomNames.NName);
      idx_C1 = strcmpi(atom_names(idx_SpinLabel),AtomNames.C1Name);
      idx_C2 = strcmpi(atom_names(idx_SpinLabel),AtomNames.C2Name);
      
      % convert to integer indices
      psf.idx_SpinLabel = find(idx_SpinLabel);
      psf.idx_O = nonzeros(psf.idx_SpinLabel.*idx_O);
      psf.idx_N = nonzeros(psf.idx_SpinLabel.*idx_N);
      psf.idx_C1 = nonzeros(psf.idx_SpinLabel.*idx_C1);
      psf.idx_C2 = nonzeros(psf.idx_SpinLabel.*idx_C2);
    end
    
    if contains(section,'NBOND')
      if ~reachedNTITLE||~reachedNATOM
        error('Section ordering in "%s" is not standard. See documentation for proper formatting', FileName)
      end
      % we are done
      break
    end
  
  end
end

end