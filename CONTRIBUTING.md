# Contributing

Please always look at this file *before* contributing to EasySpin.
In here you the coding standards and contributing guidelines as well as all the steps necessary to set up your own development environment. 

## Acceptable contributions

The following contributions are accepted:

- fix bugs for existing functions
- enhance the API or implementation of an existing function
- adds a function that is only slightly modified from a StackOverflow answer
- is tested

In the case of adding a new function or feature:

- NOT break any previously existing features and functions
- all tests must pass
- new tests must be added
- documentation must be added

## Coding style

## Set up instructions

You're going to need a few things set up correctly in order to be able to contribute to EasySpin.

First of all, you are going to need a GitHub account (if you don't already have one).
Then, you will need to know the basics of working with `git`.
If you are new to using version control software, it might be easier to use one of the graphical interfaces like [SourceTree](https://www.sourcetreeapp.com/) or [GitHub desktop](https://desktop.github.com/).
This guide assumes you have a basic understanding of `git` and know how to make a commit and create local branches.

#### The following instructions show you how to create a pull request using *SourceTree*:

1. Fork the EasySpin repo
    - go to https://github.com/StollLab/EasySpin
    - find the `Fork` button in the top right corner and click on it
    - this creates a copy of EasySpin for your GitHub (e.g. `yourusername/EasySpin`) 
2. Clone your fork to your computer
    - make sure to use the adress from *your* EasySpin repository
3. Add the original EasySpin repository to your sources so that you can get updates
    - find and click the `Settings` button in the top right corner of the tab in SourceTree
    - add a new remote:
        - Remote name: EasySpin
        - URL/path: https://github.com/StollLab/EasySpin
    - you should now have two remote repository paths
4. Link the `EasySpin/master` branch to your local `master` branch
    - right click your master branch
    - in the popup menu select `Track remote branch` and then `EasySpin/master`
    - you now can pull updates from the EasySpin repository, time to get working
5. Create a new local branch
    - give the branch a good descriptive name
6. Start coding away and creating commits to that branch
7. Once you have made a commit, you can push that commit to your forked repository
8. Create a pull request
    - go to the github website of your forked repository
    - if you do this within an hour of pushing a commit, there will be a button for that on the front page. If you can't find it, you have to click on the `Pull requests` tab and create a `New pull request`
    - Start the name of the pull 





<!-- 3. Create a branch

4. Run `npm install`
5. Run `npm t && npm run build`. If everything works, then you're ready to make changes.
6. Run `npm run test:watch`. See that it's watching your file system for changes.
7. Make your changes and try to make the tests pass. If you can't or need help then commit what you have with `--no-verify` and make a PR
8. If you get things working, add your changed files with `git add` and run `npm run commit` to get an interactive prompt for creating a commit message that follows [our standards](https://github.com/stevemao/conventional-changelog-angular/blob/master/convention.md). You'll notice that there are git hooks in place which will run testing, linting, etc. (unless you commit with `--no-verify`).
9. Push your changes to your fork with `git push`
10. Create a pull request.
11. Iterate on the solution.
12. Get merged! 🎉 🎊

## Commit messages

We follow a convention for our commit messages, to learn about why and how, see [this free egghead.io series](http://kcd.im/write-oss) -->
