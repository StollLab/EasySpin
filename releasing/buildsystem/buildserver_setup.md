# Setting up the build server for CI/CD on a Unix machine

## 1. Install and set up `perl`:

```bash
sudo apt-get install perl
```

Install packages:
```bash
apt-get install libgmp-dev
sudo cpan Net::SSH::Perl
```

If building failes with `"Could not make: unknown error"` the `make` command might be missing. Install it with 
```bash
sudo apt-get install build-essential
```
If it still does not install, force the install even if dependencies are missing or not passing all tests: 
```bash
sudo cpan install -f Net::SSH::Perl
```
 
 
## 2. Set up SSH authentication

1. Make sure `ssh` is installed
	```bash
	sudo apt-get install ssh 
	```
2. Copy key files `hostmonster_rsa` and `github_rsa` (need to create one for user) to: `~/.ssh/`.

3. Set the ssh agent up to run automatically when somebody logs in. 
Add the following to `~/.bash_profile`, if `.bash_profile` does not exist, just create a new one:
 	```bash
	SSH_ENV="$HOME/.ssh/environment"

	function start_agent {
		echo "Initialising new SSH agent..."
		/usr/bin/ssh-agent | sed 's/^echo/#echo/' > "${SSH_ENV}"
		echo succeeded
		chmod 600 "${SSH_ENV}"
		. "${SSH_ENV}" > /dev/null
		/usr/bin/ssh-add;
	}

	# Source SSH settings, if applicable

	if [ -f "${SSH_ENV}" ]; then
		. "${SSH_ENV}" > /dev/null
		ps -ef | grep ${SSH_AGENT_PID} | grep ssh-agent$ > /dev/null || {
		start_agent;
	}
	else
		start_agent;
	fi
	```
			
   Now every time you login you should get a message that the SSH agent is being initialised (additionally the Github action builder also loads the `bash_profile` to make sure everything worked)
4. To verify that the keyfiles work, load them with 
	```bash
	ssh-add ~/.ssh/hostmonster_rsa
	ssh-add ~/.ssh/github_rsa
	```
	If you get an error message, it might be necessary to change read write permissions:
	```bash
	chmod 600 ~/.ssh/hostmonster_rsa
	chmod 600 ~/.ssh/github_rsa
	```
		
	
## 3. Install TexLive
```bash
sudo apt-get install texlive-latex-extra
sudo apt-get install latex2html
```

## 4. Have MATLAB installed and verify that you can run it
```bash
matlab -help
```

## 5. Install and setup up git
```bash
sudo apt-get install git
```
	

## 6. Setup the buildsys
1. create a directory of your choosing. e.g. 
	```bash
	mkdir ~/esbuildsys
	```
2. copy the files `config.pl`, `build.pl`, `publish.pl`, `esbuild.m` and `docbuilder.pl` from the EasySpin repository into the folder that you just created.

3. Adjust the settings in `config.pl` if needed (see `config.pl` for details).

# How to build and publish EasySpin versions
## Packaging - `build.pl`
- Build all EasySpin versions that are not in the `easyspin-builds` directory and upload the most recent versions for each release channel with
	```bash
	perl build.pl
	```

- Build a specific tag
	```bash
	perl build.pl 5.3.1
	```
	This only builds the specified EasySpin version, but does not upload anything.
		
## Uploading to easyspin.org - `publish.pl`
By default, the automatic build system (which gets triggered through Github Actions), makes sure that the most recent EasySpin versions are on the webserver. If you need to upload a different version, do so with `publish.pl`. This will not only upload the packaged EasySpin, but also make sure that all the links on the website are updated accordingly.

- The `publish.pl` script *must* be called with an argument that specifies the version that you want to upload.

- The corresponding easyspin-TAG.zip file must be in the build directory (specified with config.pl)
- This example uploads the 5.2.22 version: 
	```bash
	perl publish.pl 5.2.22
	```
If you want to change the release channel, this can be specified in `config.pl`, e.g.:
```perl
my $stableversion = 5;
```
As of December 6 2019, all EasySpin major versions 5 are interpreted as meant for the release channel stable.

In this example, `publish.pl` will update all entries on the website that are tagged with stable with the version that is being uploaded.

### Development versions
Development versions _need_ to end with `alpha`/`beta`/`dev` and a number, e.g: `easyspin-6.0.0-dev.3`, `easyspin-5.2.24-beta.19`.

They are all interpreted as belonging to the development channel and everything that is tagged with development is updated on the website, regardless of the major version.

To upload 6.0.0-dev.18, you would type :
```bash
perl publish.pl 6.0.0-dev.18
```

### Uploading a manually created build
It is possible to upload a build that does not conform to the semantic versioning, e.g. `easyspin-evolve.zip`	- this could be an experimental build, coming from a specific branch. 
In such a case (when no semantic versioning is available) the releasechannel for such a file always is `experimental`.

The corresponding packaged EasySpin (`easyspin-evolve.zip`) _must_ be present in the build directory!
Upload with:
	```bash
	perl publish.pl evolve
	```