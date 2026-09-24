# Testing framework in EasySpin

The folder `tests/` contains the test suite for EasySpin. Each test is a separate function file. The test runner is `estest`.

A file counts as a test if its name contains an underscore (`_`). Other files, such as `estestexamples.m`, are ignored. Test files are named `<function>_<description>.m`, where `<function>` is the EasySpin function being tested, for example `pepper_c2h.m` or `sop_spinonehalf.m`.

## Running tests

Run `estest` from any folder, with EasySpin on the MATLAB path. `estest.m` is located in `easyspin/`, but is not included in releases.

Arguments starting with `-` are options. All other arguments are test name patterns. A pattern matches a test name exactly, unless it contains the wildcard `*`. You can give several patterns, and the options can go in any position.

```matlab
estest                  % run all tests (takes a long time!)
estest -t               % run all tests and report timings
estest pepper_c2h       % run only the test pepper_c2h
estest pepper*          % run all tests whose names start with pepper
estest *crystal*        % run all tests whose names contain crystal
estest pepper* garlic*  % run all pepper and garlic tests
estest pepper* -t       % run all pepper tests and report timings
estest pepper* -d       % run all pepper tests and plot their results
```

The function syntax works too, for example `estest('pepper*','-t')`. You can also pass a cell array of test names: `estest({'pepper_c2h','sop_spinonehalf'})`.

| Option | Effect |
| ------ | ------ |
| `-d` | Display. Tests plot or print their results. `estest` pauses after each test; press a key to continue. |
| `-t` | Report the time taken by each test, plus the total time and the 10 slowest tests. With `-d`, the time spent waiting for a keypress is not counted, but the time spent plotting is. |
| `-r` | Regenerate and store reference data for regression tests (see below). |
| `-c` | Report code coverage of the EasySpin functions. |
| `-l` | With `-c`, also list the lines that the tests didn't cover. |

Options can be combined, for example `-tc`.

### Output

`estest` prints one line per test. Each line shows the test type (`direct` or `regression`) and the outcome:

- `pass`: all subtests returned `true`.
- `failed`: at least one subtest returned `false`. The line lists the numbers of the failed subtests.
- `crashed`: the test threw an error. `estest` catches it, prints the error report, and moves on to the next test.
- `not tested`: the test returned an empty result, for example because a regression test has no reference data yet.

At the end, `estest` lists the failed and crashed tests again and prints a summary. Test names are links that open the file in the MATLAB editor. The link "rerun failed and crashed tests" runs only those tests again.

To get the results programmatically, request an output: `out = estest('pepper*')`. `out.outcomes` holds one code per test: 0 pass, 1 failed, 2 crashed, 3 not tested. `out.Results` holds the details.

## Writing a test

A test function uses one of these three signatures. The name of the function inside the file doesn't matter; `estest` calls it by its file name.

```matlab
function ok = test()                    % direct test
function ok = test(opt)                 % direct test that responds to options
function [ok,data] = test(opt,refdata)  % regression test with stored reference data
```

`ok` is either a single logical or a logical array with one element per subtest. The test passes if all elements are `true`. Returning an array rather than combining the subtests into one value lets `estest` report which subtests failed.

`opt` is a structure with these fields:

- `opt.Display`: `true` if the test should plot or print its results (`-d`). Otherwise, the test must not plot anything.
- `opt.Regenerate`: `true` if reference data is being regenerated (`-r`).
- `opt.Verbosity`: 1 with `-d` and 0 otherwise. You can pass it on to EasySpin functions to get extra log output.

### Regression tests

Some tests compare against a reference result that was calculated earlier, for example a stored simulation. For these tests:

1. The test calculates its result and returns it as `data`.
2. `estest` stores `data` in `tests/data/<testname>.mat`.
3. On later runs, `estest` loads that file and passes its contents to the test as `refdata`. The test compares its new result against it.

`estest` stores reference data only when called with `-r`. When you write a new regression test, run it once with `-r` to create the reference data:

```matlab
estest pepper_axiallw -r
```

With `-r`, `refdata` is empty. In that case, the test should return `ok = []`, so that it is reported as `not tested`. Only regenerate reference data if you are sure the current results are correct. Otherwise you overwrite a good reference with a wrong one.

See `blochsteady_simple.m` for a complete example.

### Comparing numbers

For numerical comparisons, use `areequal` (in `tests/private/`):

```matlab
ok = areequal(value,refvalue);             % true if exactly equal
ok = areequal(value,refvalue,1e-9,'rel');  % true if all differences are within 1e-9 times max(abs(refvalue))
ok = areequal(value,refvalue,1e-2,'abs');  % true if all differences are within 1e-2
```

Both inputs must be numeric arrays of the same size. If they aren't, `areequal` throws an error, and the test is reported as crashed. Prefer the relative mode (`'rel'`), because it shows how tight the test is without knowing the scale of the values.

## Tips for writing good tests

- Keep tests short and simple.
- Write many tests, so that all execution paths through a function are covered.
- Keep tests fast. If a test involves heavy calculations, reduce the number of points and similar settings. Use `estest -t` to find slow tests.
- If a test contains several subtests, return their outcomes as a logical array rather than combining them into one value. That way `estest` can report which subtests failed.
- Test functions must not plot anything unless `opt.Display` is `true`.
