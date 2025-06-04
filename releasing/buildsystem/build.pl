use strict;
use warnings;

# other dependencies
use Fcntl ':flock'; # for locking on system level
use Net::SSH::Perl; # to use ssh
use Cwd; # to get working directory
my $WorkingDir = getcwd; # get current working directory

# variables imported from config.pl
our ($SourceDir, $BuildsDir, $TempRepoDir); # directories
our ($StableMajorVersion,  $KeyForStableVersion, $MonthsToExpireStable); # settings for the stable versions
our ($DefaultMajorVersion, $KeyForDefaultVersion, $MonthsToExpireDefault); # settings for the default versions
our ($KeyForDeveloperVersion, $KeyForExperimentalVersion, $MonthsToExpireDeveloper); # settings for the experimental versions
our (@VersionCutoff, $KeyGitHub); # some other settings

print 'Reading config file\n';

require './config.pl';  # load the configuration file

# Settings
# ---------------------------------------------------------------------------------
# Options for lile locking - to ensure only one instance is running
our $NumberOfAttempts = 3; # number of Attempts to obtain a lock
our $WaitTime = 90; # time to wait between Attempts in seconds

# Creating a lock File 
# ---------------------------------------------------------------------------------
print 'Creating lock file \n';
my $LockFilename = "build.lock";
open (my $LockFile,'>'.$SourceDir.'/'.$LockFilename) or die $!;

my $Attempts = 0;
my $LockObtained = 0;

while ($Attempts < $NumberOfAttempts) {
    $LockObtained = flock $LockFile, LOCK_EX|LOCK_NB;
    if ($LockObtained) {
        last;
    }
    ++$Attempts;
    print("Another instance of build.pl appears to be running, trying again in $WaitTime seconds.\n");
    sleep($WaitTime);
}

if ($LockObtained) {
    print "Created lock File. \n";
}
else {
    print "Can not obtain lock, exiting. Maybe a build got stuck? \n";
    exit;
}

# Set up build environment
# ---------------------------------------------------------------------------------
system(qq(ssh-add $KeyGitHub)); # private key to log into GitHub, needs to be adapted to specific user

# delete and reinitialize temporary directory of EasySpin if a previous build crashed
if (-e "$TempRepoDir") {
    system("rm -rf $TempRepoDir");
}
system("mkdir $TempRepoDir");
system(qq(git clone git\@github.com:StollLab/EasySpin.git $TempRepoDir));

# Create the directory where builds are stored if not already available
unless (-e "$BuildsDir") {
    system("mkdir $BuildsDir");
}

# Delete temporary build directory where p files are encoded in if a previous build crashed
system("rm -rf /tmp/easyspin*");


# Process tag
# -----------------------------------------------------------------
my @TagsToBuild = ();

my @NewestVersion = (0, 0, 0);  # (stable, default, dev)

my $callTags = qq(git --git-dir=$TempRepoDir/.git tag);  # read tagfile
my @TagFile = `$callTags`;

# compute the numeric value of the cutoff version
my $NumericCutoff = 100000*$VersionCutoff[0]+1000*$VersionCutoff[1]+$VersionCutoff[2];

unless ($ARGV[0]) {
    # If called without an argument, all missing versions are built

    my @AvailableTags = ();

    # Process the hg AvailableTags File and grab the NumericVersion numbers, including the dev NumericVersions, but not 'tip'
    foreach (@TagFile) {
        my @vBuildID = ($_ =~ m/(v?)(\d+)\.(\d+)\.(\d+)(.*?)\s/);

        if (@vBuildID and $vBuildID[1]){
            my $NumericVersion = 100000*$vBuildID[1]+1000*$vBuildID[2]+$vBuildID[3];

            if ($NumericVersion >= $NumericCutoff) {
                my $ID = "$vBuildID[0]$vBuildID[1].$vBuildID[2].$vBuildID[3]$vBuildID[4]";
                push @AvailableTags,$ID;
            }

        }
    }

    # Get highest NumericVersion number for the three branches

    # Scan through all the AvailableTags, and compare them to the NewestVersions
    foreach (@AvailableTags) {
        my @BuildID = ($_ =~ m/v?(\d+).(\d+).(\d+)(.*)/); # match major, minor, patch and everything that follows
        my $NumericVersion = 100000*$BuildID[0]+1000*$BuildID[1]+$BuildID[2];
        
        if ($BuildID[3]){ # check if the currently processed tag is a developer NumericVersion
            my @DevVersion = ($BuildID[3] =~ m/-?([a-z]+).*?(\d{1,3})/);
            # read whether tag is alpha beta or dev
            if ($DevVersion[0] eq 'alpha') {
                $NumericVersion = $NumericVersion + 0.2 
            }
            elsif ($DevVersion[0] eq 'beta') {
                $NumericVersion = $NumericVersion + 0.3
            }
            elsif ($DevVersion[0] eq 'dev') {
                $NumericVersion = $NumericVersion + 0.1
            }

            # in case the tag also contains a numeric value, eg. dev.3, alpha1
            if ($DevVersion[1]) {
                $NumericVersion = $NumericVersion + 0.0001*$DevVersion[1];
            }

            # Update NewestVersion if necessary
            $NewestVersion[2] = $NumericVersion if $NumericVersion > $NewestVersion[2];
        }
        # if major version corresponds to the default branch:
        elsif ($BuildID[0] eq $DefaultMajorVersion) {
            $NewestVersion[1] = $NumericVersion if $NumericVersion > $NewestVersion[1];
        }
        # if major version corresponds to the stable branch:
        elsif ($BuildID[0] eq $StableMajorVersion) {
            $NewestVersion[0] = $NumericVersion if $NumericVersion > $NewestVersion[0];
        }
    }

    print "The most recent NumericVersions are @NewestVersion \n";

    # Get the currently available builds in the build directory
    opendir(my $FilesInBuildDir,$BuildsDir);
    my @AvailableBuilds = readdir($FilesInBuildDir);

    my $zipFiles = '';
    for my $File (@AvailableBuilds) {
        my @BuildID = ($File =~ m/(v?)(\d+).(\d+).(\d+)(.*).zip/);
        if (@BuildID){
            my $ID = "$BuildID[0]$BuildID[1].$BuildID[2].$BuildID[3]$BuildID[4],";
            $zipFiles = join( "", $zipFiles, $ID);
        }
    }

    print "The following builds are already in the build directory: $zipFiles \n";

    # Identify which version need to be built
    foreach (@AvailableTags) {
        unless ($zipFiles =~ m/$_,/) {
            push @TagsToBuild, $_;
        }
    }
}
else {
    # If build.pl is called with a commandline argument
    my $cmdLineArgument = $ARGV[0];
    my $TagExists = 0;

    # Try to match the commandline argument against the semantic versioning
    my @SemanticBuildID = ($cmdLineArgument =~ m/v?(\d+)\.(\d+)\.(\d+)(.*?)/);

    # If argument corresponds to semantic versioning, make sure version is newer than cutoff version
    if (@SemanticBuildID){
        my $NumericVersion = 100000*$SemanticBuildID[0]+1000*$SemanticBuildID[1]+$SemanticBuildID[2];

        if ($NumericVersion < $NumericCutoff) {
            die "Only Easyspin versions starting from $VersionCutoff[0].$VersionCutoff[1].$VersionCutoff[2] can be built using this script \n";
        }

    }

    # Check if provided tag actually exists in tag file
    foreach (@TagFile) {
        if ($_ =~ m/\b$cmdLineArgument\b/) {
            $TagExists = 1  ;
        }
    }

    # Error if tag is not existent
    unless ($TagExists) {
        die "The tag '$cmdLineArgument' was not found \n";
    }
    push @TagsToBuild, $cmdLineArgument;
}

if (@TagsToBuild == 0) {
    print("No new versions to build. \n");
}
else {
    print("The following versions will now be built: @TagsToBuild \n");
}


# Loop over all tags that should be built
# ---------------------------------------------------------------------------------
foreach (@TagsToBuild) {
    my $thisTag = $_;

    print "Building git tag $thisTag \n";

    # Update git repository
    # ---------------------------------------------------------------------------------
    # Clean cloned repo
    system("git --git-dir=$TempRepoDir/.git --work-tree=$TempRepoDir clean -f");

    # Update to tag that is being built
    system("git --git-dir=$TempRepoDir/.git --work-tree=$TempRepoDir -c advise.detachedHead=false checkout $thisTag ");

    # Generate HTML file that contains list of examples
    # ---------------------------------------------------------------------------------
    chdir($TempRepoDir.'releasing/');
    print "Creating list with examples. \n";
    system('perl mkexamples.pl');
    chdir($WorkingDir);

    # Build documentation (compile math with LaTeX etc.)
    # ---------------------------------------------------------------------------------
    print "Building documentation.\n";
    system("perl docbuilder.pl");

    # Write esbuild_config.m 
    # ---------------------------------------------------------------------------------
    my $vMatchPattern = '(v)(\d+.\d+.*)';
    my $ReleaseID;
 
    my @ShortTag = ($thisTag =~ m/$vMatchPattern/);
    
    print "ShortTag[0]: $ShortTag[0]\n";
    print "ShortTag[1]: $ShortTag[1]\n";
 
    if ($ShortTag[1]) {
        $ReleaseID = $ShortTag[1];
    }
    else {
        $ReleaseID = $thisTag;
    }
    
    my $ReleaseChannel;
    my $MonthsToExpire;

    my $MatchVersionPattern = '(v?)(\d+).(\d+).(\d+)-?([a-z]+)?[-.]?(\d+)?';
    my @thisBuildID = ($thisTag =~ m/$MatchVersionPattern/);
    if ($thisBuildID[1]) {
        if ($thisBuildID[4]) {
            $ReleaseChannel = $KeyForDeveloperVersion;
            $MonthsToExpire = $MonthsToExpireDeveloper;
        }
        elsif ($thisBuildID[1] eq $StableMajorVersion) {
            $ReleaseChannel = $KeyForStableVersion;
            $MonthsToExpire = $MonthsToExpireStable;
        }
        elsif ($thisBuildID[1] eq $DefaultMajorVersion) {
            $ReleaseChannel = $KeyForDefaultVersion;
            $MonthsToExpire = $MonthsToExpireDefault;
        }
    }
    else {  # if tag does not follow semantic versioning, e.g. easyspin-evolve.zip
        $ReleaseChannel = $KeyForExperimentalVersion;
        $MonthsToExpire = $MonthsToExpireDeveloper;
    }

    my $SourceDir = "$TempRepoDir";

    print("Writing esbuild config file\n");
    print "  ReleaseID:       $ReleaseID\n";
    print "  Release channel: $ReleaseChannel\n";
    print "  Months:          $MonthsToExpire\n";
    print "  Source dir:      $TempRepoDir\n";
    print "  Zip dest. dir:   $BuildsDir\n";
    
    my $esbuildconfigfile = "esbuild_config.m";
    open(my $Output,'>'.$esbuildconfigfile) or die("Cannot open $esbuildconfigfile!");
    print $Output "releaseID = '$ReleaseID';\n";
    print $Output "releaseChannel = '$ReleaseChannel';\n";
    print $Output "monthsToExpiry = $MonthsToExpire;\n";
    print $Output "sourceDir = '$TempRepoDir';\n";
    print $Output "zipDestDir = '$BuildsDir';\n";
    close($Output) or die("Cannot close $_!");
    
    
    # Call Matlab to run esbuild.m. This generates a packaged zip file
    # ---------------------------------------------------------------------------------
    my $MatlabOptions = '-nosplash -nodesktop -nodisplay';
    my $MatlabTarget = qq(-r "run('esbuild.m');exit;");
    print("Running MATLAB build\n");
    system('matlab '.$MatlabOptions." ".$MatlabTarget);
    system("rm $esbuildconfigfile");
    print("MATLAB build completed\n");


    # Publish zip file
    # ---------------------------------------------------------------------------------
    # Decide wether to upload version, the first conditional checks if the build version follows semantic versioning or if not (e.g. easyspin-evolve)
    my $uploadBuild = 0;
    if ($thisBuildID[1]) {
        # Translate semantic version, compare to newest version and decide wether it needs to be uploaded
        my $NumericVersion = 100000*$thisBuildID[1]+1000*$thisBuildID[2]+$thisBuildID[3];

        if ($thisBuildID[4]){ # check if is an developer NumericVersion
            if ($thisBuildID[4] eq 'alpha') {
                $NumericVersion = $NumericVersion + 0.2 
            }
            elsif ($thisBuildID[4] eq 'beta') {
                $NumericVersion = $NumericVersion + 0.3
            }
            elsif ($thisBuildID[4] eq 'dev') {
                $NumericVersion = $NumericVersion + 0.1
            }
            if ($thisBuildID[5]) {
                $NumericVersion = $NumericVersion + 0.0001*$thisBuildID[5];
            }
        }
        foreach my $VersionToCompare (@NewestVersion) {
            if ($NumericVersion == $VersionToCompare) {
                $uploadBuild = 1;
            }
        }
    }
    
    # Upload the current build by calling publish.pl
    if ($uploadBuild) {
        print "Calling publish script for upload.\n";
        system(qq(perl publish.pl $thisTag));
    }
    else {
        print "No upload.\n";
    }
}

# Clean up temporary EasySpin directories
# ---------------------------------------------------------------------------------
if (-e "$TempRepoDir") {
    system("rm -rf $TempRepoDir");
}

# Remove lock file and exit
# ---------------------------------------------------------------------------------
print "removing lock file \n";
close $LockFile;
system('rm '.$SourceDir.'/'.$LockFilename);

# ---------------------------------------------------------------------------------
print "Finished.\n";
