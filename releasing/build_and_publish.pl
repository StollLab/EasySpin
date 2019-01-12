use strict;
use warnings;

# other dependencies
use Fcntl ':flock'; # for locking on system level
use Net::SSH::Perl; # to use ssh

# ---------------------------------------------------------------------------------
# Settings

# directories
my $parentdir = '.'; # topdirectory 
my $builddir = $parentdir.'/easyspin-builds/'; # directory that builds are created in
my $uploaddir =  $parentdir.'/upload/'; # directory that will be used to upload files to easyspin.org
my $serverdir = '~/public_html/easyspin/test/'; # parentfolder on the server, used to copy files from $uploaddir into
my $repodir = $parentdir.'/easyspin/';

# tagged main versions, needs to be updated for each version
my $stableversion = 5;
my $defaultversion = 6;
my @cutoffversion = (5, 1, 0); # dont build versions lower than this (main,minor) version

# esbuild.m location
my $esbuild = "esbuild.m";

my $username = "easyspin";
my $hostname = "easyspin.org";

my $WebServerLogin = $username."@".$hostname;

my @htmlfiles = ("index.html","download.html","version.html");

# options for file locking - to ensure only one instance is running:
my $numberallowedattempts = 3; # number of attempts to obtain a lock
my $waittime = 90; # time to wait between attempts in seconds

# ---------------------------------------------------------------------------------
# creating a lock file 
my $lockfilename = "build.lock";
open (my $lockfile,'>'.$parentdir.'/'.$lockfilename) or die $!;

my $attempts = 0;
my $lockobtained = 0;

while ($attempts < $numberallowedattempts) {
    $lockobtained = flock $lockfile, LOCK_EX|LOCK_NB;
    if ($lockobtained) {
        last;
    }
    ++$attempts;
    print("Another instance of build.pl appears to be running, trying again in $waittime seconds.\n");
    sleep($waittime);
}

if ($lockobtained) {
    print "Created lock file. \n";
}
else {
    print "Can not obtain lock, exiting. Maybe a build got stuck? \n";
    exit;
}

# ---------------------------------------------------------------------------------
# Ensure ssh client is running and add keys to keychain - for this to work, the following must be added to ~/.bash_profile (create a new file if it doesnt exist):

# SSH_ENV="$HOME/.ssh/environment"

# function start_agent {
#     echo "Initialising new SSH agent..."
#     /usr/bin/ssh-agent | sed 's/^echo/#echo/' > "${SSH_ENV}"
#     echo succeeded
#     chmod 600 "${SSH_ENV}"
#     . "${SSH_ENV}" > /dev/null
#     /usr/bin/ssh-add;
# }

# # Source SSH settings, if applicable

# if [ -f "${SSH_ENV}" ]; then
#     . "${SSH_ENV}" > /dev/null
#     #ps ${SSH_AGENT_PID} doesn't work under cywgin
#     ps -ef | grep ${SSH_AGENT_PID} | grep ssh-agent$ > /dev/null || {
#         start_agent;
#     }
# else
#     start_agent;
# fi
# 
# The above code will make sure that the ssh agent is running when somebody logs in

system("ssh-add ~/.ssh/hostmonster_rsa"); # private key to log into hostmonster.com
system("ssh-add ~/.ssh/no_phrase_rsa"); # private key to log into bitbucket, needs to be adapted specific user

# ---------------------------------------------------------------------------------
# Get Tags and compare to locally available zip files, then get the list versions to build and loop over
# ---------------------------------------------------------------------------------
# updating and pulling from bitbucket.com  
system('hg pull -u -R '.$repodir);

my @tagstobuild = ();

my @newestversion = (0, 0, 0);  # (stable, default, dev)

my $callTags = qq(hg tags -R $repodir);
my @tagfile = `$callTags`;

my $NumericCutoff = 100000*$cutoffversion[0]+1000*$cutoffversion[1]+$cutoffversion[2];

unless ($ARGV[0]) {

    my @tags = ();

    # process the hg tags file and grab the version numbers, including the dev versions, but not 'tip'
    foreach (@tagfile) {
        my @buildID = ($_ =~ m/(\d+)\.(\d+)\.(\d+)(.*?)\s/);

        if (@buildID and $buildID[0]){
            my $version = 100000*$buildID[0]+1000*$buildID[1]+$buildID[2];

            if ($version >= $NumericCutoff) {
            my $ID = "$buildID[0].$buildID[1].$buildID[2]$buildID[3]";
            push @tags,$ID;
            }

        }
    }

    # get highest version number for the three branches

    # scan through all the tags, and compare them to the newestversions
    foreach (@tags) {
        my @buildID = ($_ =~ m/(\d+).(\d+).(\d+)(.*)/); # match major, minor, patch and everything that follows
        my $version = 100000*$buildID[0]+1000*$buildID[1]+$buildID[2];
        
        if ($buildID[3]){ # check if is an developer version
            my @devVersion = ($buildID[3] =~ m/-?([a-z]+).*(\d+)/);
            if ($devVersion[0] eq 'alpha') {
                $version = $version + 0.2 
            }
            elsif ($devVersion[0] eq 'beta') {
                $version = $version + 0.3
            }
            elsif ($devVersion[0] eq 'dev') {
                $version = $version + 0.1
            }

            if ($devVersion[1]) {
                $version = $version + 0.01*$devVersion[1];
            }

            $newestversion[2] = $version if $version > $newestversion[2];
        }
        elsif ($buildID[0] eq $defaultversion) {
            $newestversion[1] = $version if $version > $newestversion[1];
        }
        elsif ($buildID[0] eq $stableversion) {
            $newestversion[0] = $version if $version > $newestversion[0];
        }
    }

    print "The most recent versions are @newestversion \n";

    #  Get the currently available builds in the build directory
    opendir(my $buildfiles,$builddir);
    my @availablebuilds = readdir($buildfiles);

    my $zipfiles = '';
    for my $file (@availablebuilds) {
        my @buildID = ($file =~ m/(\d+).(\d+).(\d+)(.*).zip/);
        if (@buildID){
            my $ID = "$buildID[0].$buildID[1].$buildID[2]$buildID[3],";
            $zipfiles = join( "", $zipfiles, $ID);
        }
    }

    print "The following builds are already in the build directory: $zipfiles \n";

    # Identify which versions need to be built
    foreach (@tags) {
        unless ($zipfiles =~ m/$_,/) {
            push @tagstobuild, $_;
        }
    }
}
else {
    my $cmdlinearg = $ARGV[0];
    my $tagexists = 0;

    my @SemanticBuildID = ($cmdlinearg =~ m/(\d+)\.(\d+)\.(\d+)(.*?)/);
    print @SemanticBuildID;

    if (@SemanticBuildID){
        my $version = 100000*$SemanticBuildID[0]+1000*$SemanticBuildID[1]+$SemanticBuildID[2];

        if ($version < $NumericCutoff) {
            die "Only Easyspin versions starting from $cutoffversion[0].$cutoffversion[1].$cutoffversion[2] can be built using this script \n";
        }

    }

    foreach (@tagfile) {
        if ($_ =~ m/\b$cmdlinearg\b/) {
            $tagexists = 1  ;
        }
    }

    unless ($tagexists) {
        die "the tag '$cmdlinearg' was not found \n";
    }
    push @tagstobuild, $cmdlinearg;
}

print("The following versions will be built: @tagstobuild \n");

# ---------------------------------------------------------------------------------
# Processes all the tagstobuild
# Uploads happen only if the newest version of a release channel is being built
my $MatchPattern = '(\d+).(\d+).(\d+)-?([a-z]+)?[-.]?(\d+)?';

foreach (@tagstobuild) {
    my $thisBuild = $_;
    my @thisBuildID = ($thisBuild =~ m/$MatchPattern/);

    print "Building $thisBuild \n";

    system('hg purge -R '.$repodir);

    system("hg update $thisBuild -R $repodir -C");

    # ---------------------------------------------------------------------------------
    # Build .svg versions of formulae in the html documentation and replace the tags in html files
    print "Converting LaTeX markup in documentation to svg and modifying html files. \n";
    system('perl mkformulas.pl');

    # ---------------------------------------------------------------------------------
    # Update ReleaseID (and betaVersion) in esbuild.m 
    my $findReleaseID = "ReleaseID(.*?); \% major.minor.patch"; # pattern to find ReleaseID 
    my $thisReleaseID = "ReleaseID = \'".$thisBuild."\'; \% major.minor.patch"; # Update ReleaseID

    my $esbuildNew = $esbuild.'new';

    print("Updating ReleaseID in esbuild.m\n");
    open(my $input,'<'.$esbuild) or die("Cannot open $esbuild!");
    open(my $output,'>'.$esbuildNew) or die("Cannot open $esbuildNew!");
    while (<$input>) {
        $_ =~ s/$findReleaseID/$thisReleaseID/g;
        # do we need this? --------------------------------------
        # if ($dev){
        #     $_ =~ s/betaVersion(.*?);/betaVersion = false;/g;
        #     }
        # else {
        #     $_ =~ s/betaVersion(.*?);/betaVersion = false;/g;
        # }
        # -------------------------------------------------------
        print $output $_;    
    }

    close($input) or die("Cannot close $_!");
    close($output) or die("Cannot close $_!");

    system('mv '.$esbuildNew.' '.$esbuild);

    # ---------------------------------------------------------------------------------
    # Needs to be reworked - use pgrep call to Matlab does not work on my virtual machine linux
    # ---------------------------------------------------------------------------------
    # Call Matlab to run esbuild.m
    my $MatlabOptions = '-nosplash -nodesktop -nodisplay';
    my $MatlabTarget = qq(-r "run('$esbuild');exit;");

    # Open build dirrectory and get current number of files
    opendir(my $tempdir,$builddir) or die("Can't open $builddir: $!\n");

    my $EntriesInBuilddir = () = readdir($tempdir);
    my $CurrentEntriesInBuilddir = $EntriesInBuilddir;

    closedir($tempdir);

    print("Triggering Matlab build \n");
    system('matlab.exe '.$MatlabOptions." ".$MatlabTarget);

    print("Waiting for Matlab to finish building EasySpin $thisBuild \n");

    while ($CurrentEntriesInBuilddir == $EntriesInBuilddir) {
        opendir(my $tempdir,$builddir);
        $CurrentEntriesInBuilddir = () = readdir($tempdir);
        closedir($tempdir);
        sleep 5;
    }
    sleep 10;

    # ---------------------------------------------------------------------------------
    # Translate semantic versioning and compare decived wether it needs to be uploaded

    my $version = 100000*$thisBuildID[0]+1000*$thisBuildID[1]+$thisBuildID[2];

    if ($thisBuildID[3]){ # check if is an developer version
        if ($thisBuildID[3] eq 'alpha') {
            $version = $version + 0.2 
        }
        elsif ($thisBuildID[3] eq 'beta') {
            $version = $version + 0.3
        }
        elsif ($thisBuildID[3] eq 'dev') {
            $version = $version + 0.1
        }
        if ($thisBuildID[4]) {
            $version = $version + 0.01*$thisBuildID[4];
        }
    }

    my $upload = 0;
    foreach my $VersionToCompare (@newestversion) {
        if ($version == $VersionToCompare) {
            $upload = 1;
        }
    }

    if ($upload) {

        # copy easyspin.zip file to upload directory
        my $zipfilename = 'easyspin-'.$thisBuild.'.zip';
        print("Copying $zipfilename to upload directory. \n");
        system('cp '.$builddir.$zipfilename.' '.$uploaddir.$zipfilename);


        my $releasechannel = '';

        print "Updating website and uploading for $thisBuild \n";

        if ($thisBuildID[3]) {
            $releasechannel = 'development';
        }
        elsif ($thisBuildID[0] eq $stableversion) {
            $releasechannel = 'stable';
        }
        elsif ($thisBuildID[0] eq $defaultversion) {
            $releasechannel = 'default';
        }

        
        # ---------------------------------------------------------------------------------
        # get html files from easyspin.org and update them with the new version tags and zipfile names

        # regexp to find versions and links to zip files in the html files
        my $oldzipfile = '<!--'.$releasechannel.'-->(.*?)<!--zip-->';
        my $newzipfile = '<!--'.$releasechannel.'--><a href="easyspin-'.$thisBuild.'.zip"><!--zip-->';

        my $oldversion = '<!--'.$releasechannel.'-->(.*?)<!--version-->';
        my $newversion = '<!--'.$releasechannel.'-->'.$thisBuild.'<!--version-->';

        # loop over all the html files

        foreach (@htmlfiles) {
            my $currentfile = $_;
            print("Getting $currentfile \n");
            # download the current zipfile from easyspin.org
            system('scp '.$WebServerLogin.':'.$serverdir.$currentfile.' '.$uploaddir.$currentfile.'.bak');
            
            print("Updating index.html \n");
            open(my $inputhtml,'<'.$uploaddir.$currentfile.'.bak') or die("Cannot open $currentfile.bak!");
            open(my $outputhtml,'>'.$uploaddir.$currentfile) or die("Cannot open $currentfile!");
            while (<$inputhtml>) {
                $_ =~ s/$oldzipfile/$newzipfile/g;
                $_ =~ s/$oldversion/$newversion/g;
                print $outputhtml $_;    
            }

            close($inputhtml) or die("Cannot close $inputhtml!");
            close($outputhtml) or die("Cannot close $outputhtml!");  
        }

        print("Deleting backup versions of html files \n");
        system('rm '.$uploaddir.'*.bak');

        # upload entire upload directory to easyspin org and then clean it
        print("Uploading all new files to easyspin.org \n");
        system('scp '.$uploaddir.'* '.$WebServerLogin.':'.$serverdir);

        # clear upload directory
        print("Clean upload directory \n");
        system("rm ".$uploaddir."*");

        # ---------------------------------------------------------------------------------
        # SSH into the server, unzip the build and extract documentation
        # only happens for the 'stable' release channel!
        if ($releasechannel eq 'stable') {
            print("Logging into easyspin.org via SSH \n");
            my $ssh = Net::SSH::Perl->new($hostname);
            $ssh -> login("$username");

            my $changedirectory = "cd ".$serverdir." \n";
            my $removefolders = "rm -r ./documentation ./examples \n";
            my $unzipdocumentation = qq(unzip -qq $zipfilename 'easyspin-$thisBuild/documentation/*' -d ./tmp/ \n);
            my $unzipexamples = qq(unzip -qq $zipfilename 'easyspin-$thisBuild/examples/*' -d ./tmp/ \n);
            my $movefiles = qq(mv ./tmp/easyspin-$thisBuild/* ./ \n);
            my $rmtemp = qq(rm -r ./tmp \n);

            print("Unzipping new stable version and updating documentation and examples \n");
            my $cmd = $changedirectory.$removefolders.$unzipexamples.$unzipdocumentation.$movefiles.$rmtemp;

            my ($stdout,$stderr,$exit) = $ssh->cmd("$cmd");
        }
    }
}

# ---------------------------------------------------------------------------------
# Discard changes and purge untracked files
print "Cleaning up repository \n";
system("hg update -r default -C -R $repodir");
system("hg purge -R $repodir");

# ---------------------------------------------------------------------------------
# Clean up lockfile and exit
\print "removing lockfile \n";
close $lockfile;
system('rm '.$parentdir.'/'.$lockfilename);

print "All finished.\n";