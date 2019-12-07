# EasySpin release packaging and publishing protocol

## Things that should have been done before tagging a commit
- fix all bugs (see issue list on bitbucket)
- update release notes in `docsrc/releases.html`
- update incompatibility notes in `docsrc/releases.html`
- update installation instructions in `docsrc/installation.html`
- If they have changed, compile mex files on all platforms and 
  put `*.mex*` into the private folder
- add new functions to `docsrc/funcsalphabet.html`
- add new functions to `docsrc/funcscategory.html`

## Testing
- run complete test suite (tests/estest) on all platforms available
- run all examples

## Compilation (happens automatically on tag)
- for bugfix release: switch to stable branch
  for new release: switch to default branch
- pull and merge all changes from remote repository

Documentation (happens automatically on tag)
- recompile LaTeX source for formulas in documentation by running
  `docbuilder.pl`

Release tagging (happens automatically on tag)
- set release version tag (`ReleaseID` in `esbuild.m`)
- set expiry and horizon dates (`ExpiryData` and `HorizonDate` in 
  `releasing/esbuild.m`) - these values can be adapted in `config.pl`
- make sure all files that contain `$ReleaseID$`, `$ReleaseDate$`, `$ExpiryDate$` are
  explicitly listed in `esbuild.m` (in source tree)
  (currently `easyspininfo.m`, `eschecker.m`)

## Build (happens automatically on tag)
Execute esbuild in MATLAB, done by script `esbuild.m`:
- set release tag in documentation
- set release tag in `easyspininfo.m`
- set release date in `eschecker.m` and `easyspininfo.m`
- generate all `p` files using oldest supported Matlab on any platform
- extract file help from all m files
- set release info file
- delete all `*.m` files in private subdirectory
- delete all `*.dll` files
- check all `*.c` files for absence of `//`
- remove Mercurial folders and files
- generate example index
- zip all release files
- zip entire development tree
- archive zip and development file

## Deployment: Web site update (happens automatically on tag)
- copy new release zip file into website folder
- delete old documentation and examples subfolders
- copy documentation and examples subfolders from new release into website folder
- update `version.html`: version number (needed for update check)
- update `index.html`: version number
- update `download.html`: version number and name of zip file
- update `versions.html`: version number and name of zip file

## Deployment: Announce
- notify people that reported crucial bugs
- major release: announce to email list
- major and minor release: announce to forum
