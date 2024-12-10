# Settings

# directories
our $SourceDir = '.'; # the directory where the buildsys is supposed to run
our $ServerDir = '~/public_html/easyspin/'; # easyspin directory on the server

our $BuildsDir = $SourceDir.'/easyspin-builds/'; # directory that builds are created in, will be created if not present

# temporary directories:
our $UploadDir =  $SourceDir.'/upload/'; # directory that will be used to upload files to easyspin.org, temporary, and will be deleted aftera successful upload
our $TempRepoDir = $SourceDir.'/easyspin-temp/'; # where easyspin is cloned into

# settings of the easyspin webserver
our $hostname = "easyspin.org"; # URL of easyspin server
our $username = "easyspin"; # username for SSH
 
our $KeyWebserver = '~/.ssh/easyspin-org_rsa'; # private key to log into easyspin.org
our $KeyGitHub = '~/.ssh/github-easyspin_rsa'; # private key to log into GitHub

our @HTMLfiles = ("index.html","download.html","versions.txt"); # files on the webserver that will be updated

# assign major versions of easyspin to releasechannels
# releasechannels are only used for sorting through the html files and replacing links and version numbers
our $StableMajorVersion = 6; # major version that corresponds to the stable releasechannel
our $KeyForStableVersion = 'stable'; # keyword of this channel that is used for tagging in the html files
our $MonthsToExpireStable = 24; # months the copy of easyspin will be licenced for

our $DefaultMajorVersion = 6;
our $KeyForDefaultVersion = 'default';
our $MonthsToExpireDefault = 24;

our $KeyForDeveloperVersion = 'development'; # keyword for all development builds, example for such a build easyspin-6.0.0-dev.3.zip
our $KeyForExperimentalVersion = 'experimental'; # keyword for all builds that do not conform to semantic versioning, example for such a build; easyspin-evole.zip
our $MonthsToExpireDeveloper = 6;

our $ChannelForDocumentation = $KeyForStableVersion; # decides which release channel is used for the online documentation

our @VersionCutoff = (5, 2, 0); # versions lower than this (major,minor,patch) version are ignored and will not be automatically built
