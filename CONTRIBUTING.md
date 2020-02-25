# Contributing

Please always look at this file *before* contributing to EasySpin.

In this document you can find the [contributing](#Acceptable-contributions) and [coding](#coding-style)  guidelines, how to come up with good [commit messages and branch names](#Commit-messages-and-branch-names) as well as all the steps necessary to set up your own [development environment](#Forking-EasySpin), [staying up-to-date](#Keeping-your-EasySpin-fork-up-to-date) and [creating pull requests](#Creating-a-pull-request). 

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

## Commit messages and branch names

## Forking EasySpin

You're going to need a few things set up correctly in order to be able to contribute to EasySpin.

First of all, you are going to need a GitHub account (if you don't already have one).
Then, you will need to know the basics of working with `git`.
If you are new to using version control software, it might be easier to use one of the graphical interfaces like [SourceTree](https://www.sourcetreeapp.com/) or [GitHub desktop](https://desktop.github.com/).
This guide assumes you have a basic understanding of `git` and know how to make a commit and create local branches.

1. Fork the EasySpin repo
    - go to https://github.com/StollLab/EasySpin
    - find the `Fork` button in the top right corner and click on it
    <p align="center">
    <img width="700" src="docsrc/img/contributing-forking.jpg">
    </p>
    - this creates a copy of EasySpin for your GitHub account (e.g. `yourusername/EasySpin`) 
2. Clone the new repository with your git tool of choice
    <p align="center">
    <img width="300" src="docsrc/img/contributing-cloning.jpg">
    </p>


#### The following instructions show you how to create a pull request using *SourceTree*:

1. Fork the EasySpin repo
2. Clone your fork to your computer
    - make sure to use the adress from *your* EasySpin repository
    <p align="center">
    <img width="400" src="docsrc/img/contributing-checkout.jpg">
    </p>
3. Make sure your repository stays up to date by adding the original EasySpin (the upstream) repository
    - find and click the `Settings` button in the top right corner of the tab in SourceTree
    - add a new remote:
        - Remote name: upstream
        - URL/path: https://github.com/StollLab/EasySpin
        <p align="center">
        <img width="400" src="docsrc/img/contributing-add_easyspin_remote.jpg">
        </p>
    - you should now have two remote repository paths
4. Link the `EasySpin/master` branch to your local `master` branch
    - fetch from all repositories via the `Fetch` button
    <p align="center">
    <img width="600" src="docsrc/img/contributing-fetch.jpg">
    </p>
    - right click your master branch
    - in the popup menu select `Track remote branch` and then `upstream/master`
    <p align="center">
    <img width="500" src="docsrc/img/contributing-tracking.png">
    </p>
    - you now can pull updates from the EasySpin repository - time to get working on your code!
5. Create a new local branch
    - give it a good descriptive name (see the [examples](#Commit-messages-and-branch-names))
6. Start coding away and creating commits to that branch
7. Once you have made a commit, you can push that commit to your forked repository

#### If you prefer the CLI of `git`, this is how you fork EasySpin:
1. Fork the EasySpin repo
2. Clone your fork to your computer
    - first navigate to the folder where you want to create a copy of EasySpin
    - make sure to use the adress from *your* EasySpin repository
        ```git
        git clone git@github.com:yourusername/EasySpin.git
        ```
3. Stay up-to-date by addding the original EasySpin repository (the upstream repository) as a remote
    ```git
    git remote add upstream git@github.com:StollLab/EasySpin.git
    git fetch upstream
    ```
4. Have your local master branch track the EasySpin/master branch
    ```git
    git branch --set-upstream-to=upstream/master master
    ```
5. Create a new local branch
    - give it a good descriptive name (see the [examples](#Commit-messages-and-branch-names))
        ```git
        git checkout -b yourgood/branchname
        ```
6. Start coding away and creating commits to that branch
7. Once you have made a commit, you can push to your forked repository
    ```git
    git push origin yourgood/branchname
    ```
    or if you want to make your life easier in the future
    ```git
    git push --set-upstram origin yourgood/branchname
    ```
    and from now on, as long as you are on this branch, all you have to do to push changes to your forked repository is 
    ```git
    git push
    ```

## Keeping your EasySpin fork up-to-date
If you just forked your EasySpin, your fork will be up-to-date with the source repository (which we called upstream before). 
But if you have been working on your project for a while, there might have been a lot of changes to the original repo since you forked it, changes that put your repo out of sync.
Fortunately, it is very easy to incorporate those into your fork.

#### Using *SourceTree*:
1. checkout the master branch
    <p align="center">
    <img width="300" src="docsrc/img/contributing-syncing1.png">
    </p>
2. pull updates from github.com/stolllab/EasySpin.git
    <p align="center">
    <img width="270" src="docsrc/img/contributing-syncing2.png">
    </p>
    - now your local `master` branch is up to date
3. Let's make sure that your forked repository on GitHub is aware of those changes as well and push to it
    - make sure you push to `origin`(your repository) and not to `upstream` (the EasySpin repo)
    <p align="center">
    <img width="350" src="docsrc/img/contributing-syncing3.png">
    </p>
4. Now checkout the branch you were working on and merge the `master` branch into it
    <p align="center">
    <img width="270" src="docsrc/img/contributing-syncing4.png">
    </p>
    - if everything went well, there won't be any merge-conflicts and you now have the newest version of EasySpin to continue working with!


#### Using the CLI of `git`:
1. switch to your local `master` branch (which we set to track `upstream/master`) and fetch
    ```git
    git checkout master
    git fetch upstream
    ```
2. merge changes from the upstream repository into your `master` branch and push the new commits to your own repository
    ```git
    git merge upstream/master
    git push origin master
    ```
3. now switch back to the branch that you have been working on and merge the `master` branch into it to have the newest code available (alternatively you can use the `rebase` command, if you know how)
    ```git
    git checkout yourgood/branchname
    git merge master
    ```


## Creating a pull request
1. Create a pull request
    - go to the github website of your forked repository
    - if you do this within an hour of pushing a commit, there will be a button for that on the front page. If you can't find it, you have to click on the `Pull requests` tab and create a `New pull request`
    - Start the name of the pull request with `WIP:` to let everyone know that it is work in progress and that your code is not yet ready to be merged
    - Make sure the pull request is comparing across forks, it should look something like this:
    <p align="center">
    <img width="600" src="docsrc/img/contributing-new_PR.jpg">
    </p>
    - use the following settings:
        - Source repository: StollLab/EasySpin
        - base: master
        - head repository: yourusername/EasySpin
        - compare: the name of the branch you created
    - provide a description of what your PR does and why you think it should be added to EasySpin
    - optionally, you can also assign reviewers, assignees and labels
    - if everything looks okay, create the pull request!
2. The EasySpin community can now see what you are doing and track your progress
3. Once you are at a point where you think your contribution is ready to be merged into EasySpin, you can remove the `WIP` from the title and request a repository maintainer to review and merge your PR
    - PRs will only be merged if they pass all the tests
    - it is always a good idea to pull updates from EasySpin to your local master branch and to include these in your branch via rebase or a merge
