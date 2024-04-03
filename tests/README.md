# Testing framework in EasySpin

The folder `tests/` contains the test suite for EasySpin. Each file in this folder is a test or group of tests. The function that runs tests is `estest`.

## Running tests

To run tests, use the function `estest`. Its first input is the name of the test or group of tests. The second input are options. Here are a few example

```matlab
estest           % run all tests (could take a long time!)
estest so        % run all tests starting with so
estest pepper    % run all tests starting with pepper
estest pepper d  % run all tests starting with pepper, plotting test results if implemented
estest pepper t  % run all tests starting with pepper, reporting timings
```

There are several possible outcomes of each test function call. It could crash, which `estest` will catch and then proceed to the next test. If it doesn't crash, it returns an array of booleans containing `true` or `false` depending on whether the associated subtests passed or failed. For a test function to pass, all subtests must return `true`.

There are two different types of tests, those that are self-sufficient, and those that need external reference data. The latter are stored in `tests/data`. Both types of tests are run by `estest`.

`estest` prints a summary of test outcomes to the command window. Filenames of failed tests are linked such that they can be opened in the MATLAB editor with one click.

## Writing a test

Use one of the following three test function definitions
Unit test function syntax:

```matlab
  function ok = test()          % simple direct test
  function ok = test(opt)       % direct test that responds to options
  function [ok,data] = test(opt,refdata) % regression test with external reference data
```

Options are passed to test functions in `opt`. The most relevant is `opt.Display`. Test functions can plot test results if this variable is set to `true` and should not plot anything if it is set to `false`.

If a test needs external reference data, for example, a reference simulation, the test function should generate that data and return it in the second output argument `data`. `estest` will store the data in the `tests/data` folder and supply it the next time the test functions is called in the second input argument, `refdata`.

If you call `estest` with the `r` option, then `estest` will not record the results of the test, but only store the regenerated reference data obtained from the test function. Example

```matlab
estest pepper_edgecase r
```

There is an important utility function, `areequal()`, that tests should use for numerical comparisons. Here are some examples:

```matlab
ok = areequal(value,refvalue,1e-9,'rel');  % true if value and refvalue are within 1e-9 relative of each other
ok = areequal(value,refvalue,1e-2,'abs');  % true if value and refvalue are within 1e-9 absolute of each other
```

Preferably, use the relative mode (`'rel'`) so that it is easy to see how tight the test is without knowing the scale of the signals.

## Tips for writing good tests

- Write short tests. Keep them as simple as possible.
- Write many tests, to cover all possible execution pathways through a function.
- Write fast tests. If a test involves intense calculations, reduce the numnber of points etc. to keep the time cost at a minimum.
- If you run multiple subtests, store test outcomes in an Boolean array, and return that. `estest` can make use of this granular information.
