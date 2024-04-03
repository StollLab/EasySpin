use strict;
use warnings;

# other dependencies
use Fcntl ':flock'; # for locking on system level
use Net::SSH::Perl; # to use SSH protocol

my $Build;
if ($ARGV[0]){
    $Build = $ARGV[0];
}
else {
    die("publish.pl must be called with an argument that specifies the version to upload. \n");
}

print "Publishing\n";
print "------------------------------------------------------\n";
print "Build: $Build.\n";

# variables imported from config.pl
our ($SourceDir, $BuildsDir, $UploadDir, $ServerDir, $StableMajorVersion, $DefaultMajorVersion, $KeyForStableVersion, $KeyForDefaultVersion, $KeyForDeveloperVersion, $KeyForExperimentalVersion, $ChannelForDocumentation, $username, $hostname, @HTMLfiles, $KeyWebserver);

print "Loading config file.\n";
require './config.pl'; # load configuration file

# settings ------------------------------------------------------------------
# options for file locking - to ensure only one instance is running:
my $NumberOfAttempts = 3; # number of Attempts to obtain a lock
my $WaitTime = 90; # time to wait between Attempts in seconds

my $WebServerLogin = $username."@".$hostname;

# Creating a lock file 
# ---------------------------------------------------------------------------------
print "Creating lock file.\n";
my $LockFilename = "upload.lock";
open (my $LockFile,'>'.$SourceDir.'/'.$LockFilename) or die $!;

my $Attempts = 0;
my $LockObtained = 0;

while ($Attempts < $NumberOfAttempts) {
    $LockObtained = flock $LockFile, LOCK_EX|LOCK_NB;
    if ($LockObtained) {
        last;
    }
    ++$Attempts;
    print("Another instance of publish.pl appears to be running, trying again in $WaitTime seconds.\n");
    sleep($WaitTime);
}

if ($LockObtained) {
    print "Created lock file. \n";
}
else {
    print "Cannot obtain lock, exiting. \n";
    exit;
}

# Add key to hostmonster to keychain
# ---------------------------------------------------------------------------------
system("ssh-add $KeyWebserver"); # private key to log into hostmonster.com


# Set up environment
# ---------------------------------------------------------------------------------
if (-e "$UploadDir") {
    system("rm -r $UploadDir");
}

system("mkdir $UploadDir");

# Determine the release channel
# ---------------------------------------------------------------------------------
print "Determining release channel from build ID.\n";
my $ReleaseChannel;
my $MatchPattern = '(\d+).(\d+).(\d+)-?([a-z]+)?[-.]?(\d+)?';

my @BuildID = ($Build =~ m/$MatchPattern/);

if ($BuildID[0]) {
    if ($BuildID[3]) {
        $ReleaseChannel = $KeyForDeveloperVersion;
    }
    elsif ($BuildID[0] eq $StableMajorVersion) {
        $ReleaseChannel = $KeyForStableVersion;
    }
    elsif ($BuildID[0] eq $DefaultMajorVersion) {
        $ReleaseChannel = $KeyForDefaultVersion;
    }
}
else {
    # if tag does not follow semantic versioning, e.g. easyspin-evolve.zip
    $ReleaseChannel = $KeyForExperimentalVersion;
}
print "  Release channel: $ReleaseChannel.\n";

my $ZipBuild;

if ($BuildID[0]) {
    my $vMatchPattern = 'v(.*)';
    my @ShortTag = ($Build =~ m/$vMatchPattern/);
    if (@ShortTag[0]) {
        $ZipBuild = $Build;
        $Build = @ShortTag[0];
    }
    else {
        $ZipBuild = $Build;
    }
}
else {
    $ZipBuild = $Build;
}

# Look for zip file with the provided tag
# ---------------------------------------------------------------------------------
my $zipFileName = 'easyspin-'.$ZipBuild.'.zip';
my $NewzipFileName = 'easyspin-'.$Build.'.zip';

if (-e "$BuildsDir$zipFileName") {
    print("Copying $zipFileName to upload directory. \n");
    system('cp '.$BuildsDir.$zipFileName.' '.$UploadDir.$NewzipFileName);
}
else {
    die("$zipFileName does not exist in $BuildsDir \n");
}

# Regexp to find versions and links to zip files in the html files
# ---------------------------------------------------------------------------------
my $findLinkTozipFile = '<!--'.$ReleaseChannel.'zip-->';
my $matchzipFile = "easyspin-(.*?).zip";
my $replacezipFile = "easyspin-$Build.zip";my $replaceLinkTozipFile = '<!--'.$ReleaseChannel.'--><a href="easyspin-'.$Build.'.zip"><!--zip-->';

my $matchOldVersion = '<!--'.$ReleaseChannel.'-->(.*?)<!--version-->';
my $replaceOldVersion = '<!--'.$ReleaseChannel.'-->'.$Build.'<!--version-->';

my $matchInVersionsFile = "$ReleaseChannel:.*";
my $replaceInVersionsFile = "$ReleaseChannel:$Build";

# Get html files from easyspin.org and update them with the new version tags and zipfile names
# ---------------------------------------------------------------------------------
print "Update HTML files.\n";
foreach (@HTMLfiles) {
    my $CurrentFile = $_;
    print("Getting $CurrentFile \n");

    # download the current zipfile from easyspin.org
    system('scp '.$WebServerLogin.':'.$ServerDir.$CurrentFile.' '.$UploadDir.$CurrentFile.'.bak');
    
    # scan through the html file and replace strings
    print("Updating $CurrentFile \n");
    open(my $InputHTML,'<'.$UploadDir.$CurrentFile.'.bak') or die("Cannot open $CurrentFile.bak!");
    open(my $OutputHTML,'>'.$UploadDir.$CurrentFile) or die("Cannot open $CurrentFile!");
    while (<$InputHTML>) {

        if ($_ =~ m/$findLinkTozipFile/) {
            $_ =~ s/$matchzipFile/$replacezipFile/g;
        }
        $_ =~ s/$matchOldVersion/$replaceOldVersion/g;
        $_ =~ s/$matchInVersionsFile/$replaceInVersionsFile/g;
        print $OutputHTML $_;    
    }

    close($InputHTML) or die("Cannot close $InputHTML!");
    close($OutputHTML) or die("Cannot close $OutputHTML!");  
}

print("Deleting backup versions of html files \n");
system('rm '.$UploadDir.'*.bak');

# Upload entire upload directory to easyspin org and then clean it
# ---------------------------------------------------------------------------------
print("Uploading all new files to easyspin.org \n");

system('scp '.$UploadDir.'* '.$WebServerLogin.':'.$ServerDir);

# clear upload directory
print("Clear upload directory \n");
system("rm -R $UploadDir");

# SSH into the server, unzip the build and extract documentation
# ---------------------------------------------------------------------------------
# only happens for the release channel that is specified with $ChannelForDocumentation in the config file

if ($ReleaseChannel eq $ChannelForDocumentation) {
    print("Logging into easyspin.org via SSH \n");
    my $SSHSession = Net::SSH::Perl->new($hostname);
    $SSHSession -> login("$username");

    my $changeDir = "cd ".$ServerDir." \n";
    my $rmFolders = "rm -rf ./documentation ./examples \n";
    my $unzipDoc = qq(unzip -qq $NewzipFileName 'easyspin-$Build/documentation/*' -d ./tmp/ \n);
    my $unzipExamples = qq(unzip -qq $NewzipFileName 'easyspin-$Build/examples/*' -d ./tmp/ \n);
    my $moveFiles = qq(cp -r ./tmp/easyspin-$Build/* ./ \n);
    my $rmTempDir = qq(rm -r ./tmp \n);

    print("Unzipping new stable version and updating documentation and examples \n");
    my $IssueCmd = $changeDir.$rmFolders.$unzipExamples.$unzipDoc.$moveFiles.$rmTempDir;

    # send command
    my ($STDOut,$STDErr,$Exit) = $SSHSession->cmd($IssueCmd);

}

# Clean up lock file and exit
# ---------------------------------------------------------------------------------
print "removing lock file \n";
close $LockFile;
system('rm '.$SourceDir.'/'.$LockFilename);

# ---------------------------------------------------------------------------------
print "All finished.\n";
