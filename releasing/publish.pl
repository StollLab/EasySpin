use strict;
use warnings;

# other dependencies
use Fcntl ':flock'; # for locking on system level
use Net::SSH::Perl; # to use ssh

my $Build;

if ($ARGV[0]){
    $Build = $ARGV[0];
}
else {
    die("publish.pl must be called with an argument that specifies the version to upload. \n");
}


# values imported from config.pl
our ($parentdir, $builddir, $uploaddir, $serverdir, $repodir, $stableversion, $defaultversion, $stabletag, $defaulttag, $devtag, $exptag, @cutoffversion, $channelfordocumentation, $username, $hostname, @htmlfiles);
require './config.pl';

# settings ------------------------------------------------------------------
# options for file locking - to ensure only one instance is running:
my $numberallowedattempts = 3; # number of attempts to obtain a lock
my $waittime = 90; # time to wait between attempts in seconds

my $WebServerLogin = $username."@".$hostname;

# ---------------------------------------------------------------------------------
# creating a lock file 
my $lockfilename = "upload.lock";
open (my $lockfile,'>'.$parentdir.'/'.$lockfilename) or die $!;

my $attempts = 0;
my $lockobtained = 0;

while ($attempts < $numberallowedattempts) {
    $lockobtained = flock $lockfile, LOCK_EX|LOCK_NB;
    if ($lockobtained) {
        last;
    }
    ++$attempts;
    print("Another instance of publish.pl appears to be running, trying again in $waittime seconds.\n");
    sleep($waittime);
}

if ($lockobtained) {
    print "Created lock file. \n";
}
else {
    print "Can not obtain lock, exiting. \n";
    exit;
}

# ---------------------------------------------------------------------------------
# add keys to keychain
system("ssh-add ~/.ssh/hostmonster_rsa"); # private key to log into hostmonster.com

if (-e "$uploaddir") {
    system("rm -R $uploaddir");
}

system("mkdir $uploaddir");

my $zipfilename = 'easyspin-'.$Build.'.zip';

if (-e "$builddir$zipfilename") {
    print("Copying $zipfilename to upload directory. \n");
    system('cp '.$builddir.$zipfilename.' '.$uploaddir.$zipfilename);
}
else {
    die("$zipfilename does not exist in $builddir \n");
}

my $releasechannel;
my $MatchPattern = '(\d+).(\d+).(\d+)-?([a-z]+)?[-.]?(\d+)?';

my @BuildID = ($Build =~ m/$MatchPattern/);

if ($BuildID[0]) {
    if ($BuildID[3]) {
        $releasechannel = $devtag;
    }
    elsif ($BuildID[0] eq $stableversion) {
        $releasechannel = $stabletag;
    }
    elsif ($BuildID[0] eq $defaultversion) {
        $releasechannel = $defaulttag;
    }
}
else {
    $releasechannel = $exptag;
}

# regexp to find versions and links to zip files in the html files
my $oldzipfile = '<!--'.$releasechannel.'-->(.*?)<!--zip-->';
my $newzipfile = '<!--'.$releasechannel.'--><a href="easyspin-'.$Build.'.zip"><!--zip-->';

my $oldversion = '<!--'.$releasechannel.'-->(.*?)<!--version-->';
my $newversion = '<!--'.$releasechannel.'-->'.$Build.'<!--version-->';


# ---------------------------------------------------------------------------------
# get html files from easyspin.org and update them with the new version tags and zipfile names

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
print("Clear upload directory \n");
system("rm -R $uploaddir");

# ---------------------------------------------------------------------------------
# SSH into the server, unzip the build and extract documentation
# only happens for the 'stable' release channel!
if ($releasechannel eq $channelfordocumentation) {
    print("Logging into easyspin.org via SSH \n");
    my $ssh = Net::SSH::Perl->new($hostname);
    $ssh -> login("$username");

    my $changedirectory = "cd ".$serverdir." \n";
    my $removefolders = "rm -r ./documentation ./examples \n";
    my $unzipdocumentation = qq(unzip -qq $zipfilename 'easyspin-$Build/documentation/*' -d ./tmp/ \n);
    my $unzipexamples = qq(unzip -qq $zipfilename 'easyspin-$Build/examples/*' -d ./tmp/ \n);
    my $movefiles = qq(mv ./tmp/easyspin-$Build/* ./ \n);
    my $rmtemp = qq(rm -r ./tmp \n);

    print("Unzipping new stable version and updating documentation and examples \n");
    my $cmd = $changedirectory.$removefolders.$unzipexamples.$unzipdocumentation.$movefiles.$rmtemp;

    my ($stdout,$stderr,$exit) = $ssh->cmd("$cmd");
}

# ---------------------------------------------------------------------------------
# Clean up lockfile and exit
print "removing lockfile \n";
close $lockfile;
system('rm '.$parentdir.'/'.$lockfilename);

print "All finished.\n";