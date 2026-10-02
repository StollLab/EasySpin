# CLAUDE.md

EasySpin is a MATLAB toolbox for simulating and analyzing EPR spectra and magnetometry data. See [README.md](README.md), [CONTRIBUTING.md](CONTRIBUTING.md), [tests/README.md](tests/README.md) and [releasing/README.md](releasing/README.md) for the full guides.

## Layout

- `easyspin/`: public toolbox functions, one per file. `functionSignatures.json` holds tab-completion signatures for the public functions, and the test `function_signatures_completeness` checks that it is complete.
- `easyspin/private/`: internal helpers that only `easyspin/*.m` can call. The release build deletes the `.m` files in this folder, keeping only p-code.
- `tests/`: one test per file, named `<function>_<description>.m`. Regression reference data is in `tests/data/`, and test helpers such as `areequal` are in `tests/private/`. Tests reach private helpers through `runprivate('name',args...)`.
- `docsrc/`: hand-written HTML documentation, one page per function (`pepper.html` etc.), plus `funcsalphabet.html`, `funcscategory.html` and `releases.html`.
- `examples/`: example scripts grouped by topic.
- `releasing/`: build and doc scripts (`docbuilder.pl`, `buildsystem/esbuild.m`). The build is triggered on tag by `.github/workflows/continuous_deployment.yml`.
- `documentation/` (built docs) and `scratch/` are gitignored.

## Running tests

Run tests in MATLAB with `easyspin/` on the path. The runner is `estest` ([easyspin/estest.m](easyspin/estest.m)).

```matlab
estest pepper_c2h       % one test
estest pepper*          % all tests matching a pattern
estest pepper* -t       % with timings
estest blochsteady_simple -r   % regenerate regression reference data
```

From the shell (MATLAB is on PATH):

```sh
matlab -batch "addpath('easyspin'); estest pepper*"
```

- **Never run the full suite (plain `estest`).** It takes very long. Run only the patterns for the functions you changed, plus the callers that are affected.
- Only regenerate reference data (`-r`) when you are sure the new results are correct.
- Tests return a logical scalar or array (one element per subtest). They must not plot unless `opt.Display` is true. Use `areequal(a,b,tol,'rel')` for numeric comparisons.

- For private functions, use `runprivate`.
- Keep tests short. Tests should run fast.
- If a test involves a spin system, keep it as small as possible.
- Tests for `cardamom`, `saffron` and `spidyan` are very slow. Avoid running them unless there are edits to those functions.


## Code conventions

- **Target MATLAB R2021b.** Don't use newer language features. `chkmlver` enforces R2021b as the minimum, so version checks for older releases are dead code.
- Do not use functions from toolboxes outside core MATLAB.
- Use 2-space indentation and camelCase names. Readability matters more than performance or cleverness.
- Cite the literature (DOI and equation number) when you implement equations from a paper.
- Every public function starts with a help block (`% name  one-line description`, then usage, inputs and outputs) and `if nargin==0, help(mfilename); return; end`, unless it is a GUI function.
- Inputs come as structures: `Sys` (spin system), `Exp` (experiment), `Opt` (options). Common helpers:
  - `validatespinsys` normalizes and checks `Sys`.
  - `validate_exp` and the `p_*.m` helpers (`p_sweeprange`, `p_temperature`, `p_sampletype`, ...) parse shared `Exp` fields consistently across simulation functions. Reuse them instead of re-parsing fields.
  - `adddefaults(User,Default)` merges defaults into a structure, and `parseoption` parses enumerated string options.
  - `logmsg(level,fmt,...)` writes log output, controlled by `Opt.Verbosity`.
- Simulation functions (`pepper`, `garlic`, `chili`, `salt`, `saffron`, ...) loop over components and isotopologues through `compisoloop`, which calls the function recursively with `Sys.singleiso`.
- MEX sources (`*.c`) in `easyspin/private/` come with precompiled binaries for Windows, Linux and macOS (`.mexw64`, `.mexa64`, `.mexmaci64`). If you change C code, recompile with `easyspin_compile` and avoid `//` comments.

## When adding or changing a public function

- Add or update tests in `tests/`.
- Update the doc page in `docsrc/` and the header help.
- For a new function, also add it to `docsrc/funcsalphabet.html`, `docsrc/funcscategory.html` and `easyspin/functionSignatures.json`.
- When changing a function interface, also update `easyspin/functionSignatures.json`.
- Note user-visible changes and incompatibilities in `docsrc/releases.html`. Keep entries concise, only mention significant behavior changes that are relevant to the user.

## Physics

- Microwave polarization is fixed in the lab frame.
- Positive B points along z(Lab), negative B points along -z(Lab).
- Spin systems can also have angular momentum, via `Sys.L`.

## Commits

- Do not automatically commit anything.
- Use short, specific, lowercase commit messages, prefixed with the affected function where that applies (`resfreqs_perturb: fix bug in g strain calculation`). Reference issues as `(closes #N)`. Keep separate lines of work in separate commits.
