#!/usr/bin/perl

#----------------------------------------------------------
# Generates a file scripts/funclist.txt containing
# all EasySpin function names, i.e., m files
# from the toolbox directory.
#----------------------------------------------------------

# Get list of EasySpin functions
chdir('../easyspin');
@files = glob('*.m');
map { s/\.m//} @files;

# Generate file with function list
open(D,'>../scripts/funclist.txt');
$string = join("\n",@files);
print D $string;

# Generate list for doc function input form
#foreach $st (@files) {
#  print 'funList[++i]="'.$st.'"',";\n";
#}

$length = scalar(@files);
print "$length EasySpin functions\n";
