# Settings

# directories
our $parentdir = '.'; # the directory where the buildsys is supposed to run
our $serverdir = '~/public_html/easyspin/test/'; # easyspin directory on the server

our $builddir = $parentdir.'/easyspin-builds/'; # directory that builds are created in, will be created if not present

# temporary directories:
our $uploaddir =  $parentdir.'/upload/'; # directory that will be used to upload files to easyspin.org, temporary, and will be deleted aftera successful upload
our $repodir = $parentdir.'/easyspin-temp/'; # where easyspin is cloned into

# settings of the easyspin webserver
our $hostname = "easyspin.org"; # URL of easyspin server
our $username = "easyspin"; # username for SSH

our @htmlfiles = ("index.html","download.html","version.html"); # files on the webserver that will be updated

# assign major versions of easyspin to releasechannels
# releasechannels are only used for sorting through the html files and replacing links and version numbers
our $stableversion = 5; # major version that corresponds to the stable releasechannel
our $stabletag = 'stable'; # keyword of this channel that is used for tagging in the html files

our $defaultversion = 6;
our $defaulttag = 'default';

our $devtag = 'development'; # keyword for all development builds, example for such a build easyspin-6.0.0-dev.3.zip
our $exptag = 'experimental'; # keyword for all builds that do not conform to semantic versioning, example for such a build; easyspin-evole.zip

our $channelfordocumentation = $stabletag; # decides which release channel is used for the online documentation

our @cutoffversion = (5, 1, 0); # versions lower than this (major,minor,patch) version are ignored and will not be automatically built

# name and location of the esbuild script
our $esbuild = "./esbuild.m";





