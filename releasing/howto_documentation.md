# Making a local copy of the EasySpin documentation

On release, the EasySpin documentation is automatically built from the the html files in the `/docsrc/` folder.

If you change the documentation you might want to verify that your
modifications were implemented correctly. For this you will have to build the documentation locally (don't worry, the built documentation will notbe added to the repository, as it is excluded through the `.gitignore` file).

**Adding math in the documentation:** If you want to include equations in the documentation
you have to do so through the use of LaTeX markup, which will be
automatically converted to a png file.

## Requirements
In order to build the documentation you need to have perl and latex installed.

If you are on Windows:
* Install perl (e.g. get it from http://strawberryperl.com/)
* Verify perl is on the search path. To do so, open a terminal and run `perl -v`. If you get an error message, you probably need to add the installation path manually to the Windows path defintion ([see here for how](https://docs.alfresco.com/4.2/tasks/fot-addpath.html)).
* Install LaTeX (e.g. from https://miktex.org/)
* Make sure it was correctly added to the search path by running `pdftex -version`. If you get an error message, add the directory that contains the LaTex files to the search path.

If you are on a Unix system:
* Install perl with `sudo apt-get install perl`
* Install LaTeX with `sudo apt-get install texlive-latex-extra` and `sudo apt-get install latex2html`

## Building the documentation
Once you have made sure, that you have working versions of perl and LaTeX, open a terminal *in* the `releasing` folder.
Then run:
```
perl docbuilder.pl
```
This will create the compiled version of the EasySpin documentation in the `documentation` folder. 






