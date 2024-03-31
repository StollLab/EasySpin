# EasySpin release packaging and publishing protocol

## Things that should have been done before tagging a commit
- fix as many bugs as possible (see issue list on GitHub)
- update release notes in `docsrc/releases.html`
- update incompatibility notes in `docsrc/releases.html`
- update installation instructions in `docsrc/installation.html`
- if they have changed, compile mex files on all platforms and 
  put `*.mex*` into the private folder
- add new functions to `docsrc/funcsalphabet.html`
- add new functions to `docsrc/funcscategory.html`
- add new functions to `easyspin/functionSignatures.json`

## Testing
- run complete test suite (`tests/estest`) on all platforms available
- run all examples

## Compilation (happens automatically on tag)
- for bugfix release: switch to stable branch; 
  for new release: switch to default branch
- pull and merge all changes from remote repository

Documentation (happens automatically on tag)
- recompile LaTeX source for formulas in documentation by running
  `docbuilder.pl`

Release tagging (happens automatically on tag)
- set release version tag (`ReleaseID` in `esbuild.m`)
- set expiry date (`ExpiryData` in 
  `releasing/esbuild.m`) - this value can be adapted in `config.pl`
- make sure all files that contain `$ReleaseID$`, `$ReleaseDate$`, `$ExpiryDate$` are
  explicitly listed in `esbuild.m` (in source tree)
  (currently `easyspin_info.m`, `eschecker.m`)

## Build (happens automatically on tag)
Execute `esbuild` script in MATLAB, which does the following:
- set release tag in documentation
- set release tag in `easyspin_info.m`
- set release date in `eschecker.m` and `easyspin_info.m`
- generate all `p` files using oldest supported Matlab on any platform
- extract file help from all m files
- set release info file
- delete all `*.m` files in private subdirectory
- delete all `*.dll` files
- check all `*.c` files for absence of `//`
- generate example index (`mkexamples.pl`)
- zip all release files
- zip entire development tree
- archive release and development zip files

## Deployment: Web site update (happens automatically on tag)
- copy new release zip file into website folder
- delete old documentation and examples subfolders
- copy documentation and examples subfolders from new release into website folder
- update `index.html`: version number
- update `download.html`: version number and name of zip file
- update `versions.txt`: version number and name of zip file (needed for update check)

## Deployment: Announce
- major and minor release: announce to forum
- bug release: no announcement
