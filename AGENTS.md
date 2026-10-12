# EasySpin development guide for coding agents

EasySpin is a MATLAB toolbox for simulating and analyzing EPR spectra and magnetometry data. This file covers what agents need. Read [tests/README.md](tests/README.md) only when writing new tests or regression tests, and [releasing/README.md](releasing/README.md) only for build or release work. [CONTRIBUTING.md](CONTRIBUTING.md) is for human contributors and is covered here.

## Layout

- `easyspin/`: public toolbox functions, one per file. `functionSignatures.json` holds tab-completion signatures for the public functions.
- `easyspin/private/`: internal helpers that only `easyspin/*.m` can call. The release build deletes the `.m` files in this folder, keeping only p-code.
- `tests/`: one test per file, named `<function>_<description>.m`. Regression reference data is in `tests/data/`, and test helpers such as `areequal` are in `tests/private/`.
- `docsrc/`: hand-written HTML documentation (see [Documentation](#documentation)).
- `examples/`: example scripts grouped by topic.
- `releasing/`: build and doc scripts (`docbuilder.pl`, `buildsystem/esbuild.m`). The build is triggered on tag by `.github/workflows/continuous_deployment.yml`.
- `documentation/` (built docs) and `scratch/` are gitignored.

## Physics conventions

- Units: couplings (`A`, `D`, `J`, `ee`, `Q`, `nn`, ...) in MHz, `Exp.mwFreq` in GHz, magnetic fields in mT, Euler angles in radians, temperatures in K.
- Microwave polarization is fixed in the lab frame.
- Positive B points along z(Lab), negative B points along -z(Lab).
- Spin systems can also have angular momentum, via `Sys.L`.

## Code

### Style

- **Target MATLAB R2021b.** Don't use newer language features. `chkmlver` enforces R2021b as the minimum, so version checks for older releases are dead code.
- Do not use functions from toolboxes outside core MATLAB.
- Use 2-space indentation and camelCase names. Readability matters more than performance or cleverness.
- Cite the literature (DOI and equation number) when you implement equations from a paper.
- Use analytical derivatives in toolbox code. Finite differences are only for tests, as independent checks.
- Share code between related functions (e.g. the matrix solvers `resfields`/`resfreqs_matrix` and the perturbation solvers `*_perturb`), but keep private helpers few, with simple interfaces (e.g. pass `Sys` and arrays over all orientations). Keep few-line steps inline instead of adding one-line helpers.

### Structure and shared helpers

- Every public function starts with a help block and `if nargin==0, help(mfilename); return; end`, unless it is a GUI function (see [Documentation](#documentation)).
- Inputs come as structures: `Sys` (spin system), `Exp` (experiment), `Opt` (options). Common helpers:
  - `validatespinsys` normalizes and checks `Sys`.
  - `validate_exp` and the `p_*.m` helpers (`p_sweeprange`, `p_temperature`, `p_sampletype`, ...) parse shared `Exp` fields consistently across simulation functions. Reuse them instead of re-parsing fields.
  - `adddefaults(User,Default)` merges defaults into a structure, and `parseoption` parses enumerated string options.
  - `logmsg(level,fmt,...)` writes log output, controlled by `Opt.Verbosity`.
- Simulation functions (`pepper`, `garlic`, `chili`, `salt`, `saffron`, ...) loop over components and isotopologues through `compisoloop`, which calls the function recursively with `Sys.singleiso`.
- MEX sources (`*.c`) in `easyspin/private/` come with precompiled binaries for Windows, Linux and macOS (`.mexw64`, `.mexa64`, `.mexmaci64`). If you change C code, recompile with `easyspin compile` and avoid `//` comments.

## Tests

### Writing tests

- Tests return a logical scalar or array (one element per subtest). They must not plot unless `opt.Display` is true.
- Use `areequal(a,b,tol,'rel')` for numeric comparisons.
- Reach private helpers through `runprivate('name',args...)`.
- Keep tests short and fast, with spin systems as small as possible.

### Running tests

Run tests in MATLAB with `easyspin/` on the path. The runner is `estest` ([easyspin/estest.m](easyspin/estest.m)).

```matlab
estest pepper_c2h       % one test
estest pepper* garlic*  % all tests matching any of the patterns
estest pepper* -t       % with timings
estest pepper_c2h -d    % display results, for investigating a failure
estest blochsteady_simple -r   % regenerate regression reference data
```

From the shell (MATLAB is on PATH), MATLAB startup dominates the run time, so put all patterns in a single call:

```sh
matlab -batch "addpath('easyspin'); estest resfields* pepper* eigfields*"
```

- **Never run the full suite (plain `estest`).** It takes very long. Run only the patterns for the functions you changed, plus the callers that are affected. Note that `estest * -t` is parsed as `estest()*(-t)` and runs the full suite.
- Tests for `cardamom`, `saffron` and `spidyan` are very slow. Avoid running them unless there are edits to those functions.
- Only regenerate reference data (`-r`) when you are sure the new results are correct.

## Documentation

- **Header help:** every public function starts with a help block: `% name  one-line description`, then usage, inputs and outputs.
- **Doc pages:** `docsrc/` has one hand-written HTML page per public function (`pepper.html` etc.). `docsrc/funcsalphabet.html` and `docsrc/funcscategory.html` list all functions.
- **Audience:** the docs are for users without advanced math. Use plain language and practical rules of thumb (e.g. "keep the FWHM below about 15°"), and leave derivations to code comments and tests.
- **Tab completion:** `easyspin/functionSignatures.json` holds the signatures of the public functions. The test `function_signatures_completeness` checks that it is complete.
- **Release notes:** add an entry to `docsrc/releases.html` only for significant user-visible behavior changes and incompatibilities. Skip small feature extensions, latent-bug fixes and stricter validation of unusual input, and ask when unsure. Keep entries concise.

## Checklist: adding or changing a public function

- Add or update tests in `tests/`.
- Update the header help, the doc page and the user guides in `docsrc/`.
- When changing the interface, update `easyspin/functionSignatures.json`.
- For a new function, also add it to `docsrc/funcsalphabet.html`, `docsrc/funcscategory.html` and `easyspin/functionSignatures.json`.
- Add a release note only if the change meets the bar above.

## Commits

- Do not automatically commit anything.
- Use short, specific, lowercase commit messages, prefixed with the affected function where that applies (`resfreqs_perturb: fix bug in g strain calculation`). Reference issues as `(closes #N)`. Keep separate lines of work in separate commits, unless changes are minimal.
