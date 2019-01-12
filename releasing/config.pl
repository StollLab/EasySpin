# Settings

# directories
our $parentdir = '.'; # topdirectory 
our $builddir = $parentdir.'/easyspin-builds/'; # directory that builds are created in
our $uploaddir =  $parentdir.'/upload/'; # directory that will be used to upload files to easyspin.org
our $serverdir = '~/public_html/easyspin/test/'; # parentfolder on the server, used to copy files from $uploaddir into
our $repodir = $parentdir.'/easyspin-temp/';

# tagged main versions, needs to be updated for each version
our $stableversion = 5;
our $stabletag = 'stable';

our $defaultversion = 6;
our $defaulttag = 'default';

our $devtag = 'development';
our $exptag = 'experimental';

our $channelfordocumentation = $stabletag;

our @cutoffversion = (5, 1, 0); # dont build versions lower than this (main,minor) version

# esbuild.m location
our $esbuild = "esbuild.m";

our $username = "easyspin";
our $hostname = "easyspin.org";

our @htmlfiles = ("index.html","download.html","version.html");

