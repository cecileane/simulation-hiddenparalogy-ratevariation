#!/usr/bin/perl

# Perl script to run iqtree + astral for each locus.
# this is *not* part of the TICR pipeline
# adapted from Solis-Lemus, Yang and Ane (2016) scripts:
# https://github.com/crsl4/InconsistencySpeciesTreeGeneFlow/blob/master/scripts/estGeneTrees/raxml.pl

# model: HKY + Gamma (for rate variation across sites)

# usage:
#
# iqtree.pl --seqdir=xxx/yyy --iqtreedir=xxx --iqtreedir=xxx
#
# other options:
# --numCores = number of cores to use by iqtree. default 6
# --No boostrap 
# --numboot = number of bootstrap reps. default 100
# --convert2phylip (default) or --noconvert2phylip to convert nexus input files to phylip format
# --doastral (default) or --nodoastral to do or not do ASTRAL at the end.
#
# the script will create a log file in iqtreedir/iqtree.pl.log

# warning: assuming all paths (seqdir and iqtreedir) are relative paths, and linux/mac

use Getopt::Long;
use File::Path qw( make_path );
use strict;
use warnings;
use Carp;

use Cwd;
#use lib '/u/c/l/claudia/lib/perl/lib/site_perl/'; # we need this because I had to install locally the Statistics module
#use Statistics::R;


# ================= parameters ======================
my $currentdir = `pwd`;
chomp $currentdir; # remove new line at the end of the current path 
my $boot = 1; # boot=0 is not implemented, in fact!
my $numboot = 100;
my $seed_iqtree;
my $numCores = 1; # The script got multi-processed through simulation_iqtree.jl 
my $seqdir;   # directory where sequences are
my $phylipdir;
my $iqtreedir;
# my $astraldir;
my $convertphylip = 1;
my $doastral = 1;  
# my $astral = $currentdir . '/executables/astral.5.7.8.jar'; # Old astral version 
# my $astral = $currentdir . '/executables/astral-pro3'; # Astral-pro output trees in substitution unit 
# my $astral = $currentdir . '/executables/wastral'; 
# Above: weights input gene trees based on branch length and output branch length in CU and PP
my $iqtree = $currentdir . '/executables/iqtree2';

# -------------- read arguments from command-line -----------------------
GetOptions( 'numboot=i' => \$numboot,
        'seed_iqtree=i' => \$seed_iqtree,
	    'boot!' => \$boot,
	    'numCores=i' => \$numCores,
	    'seqdir=s' => \$seqdir,
	    'iqtreedir=s' => \$iqtreedir,
	    # 'astraldir=s' => \$astraldir,
	    'convert2phylip!' => \$convertphylip,
	    # 'doastral!' => \$doastral,
        'numCores=i' => \$numCores
    );

die "seqdir not defined or not a directory" if (!(defined $seqdir) or !(-d $seqdir));
# sequence directory should have been uncompressed. If not, add this:
# system("tar -zxvf $seqdir");
if ($convertphylip) {
    $phylipdir = $seqdir;
    $phylipdir =~ s/[^\/]+$//; # remove last directory
    $phylipdir .= "phylip";    # replace by 'phylip'
    print "directory for phylip files: $phylipdir\n";
    make_path $phylipdir unless(-d $phylipdir);
}

die "a directory for IQtree output should be specified with the --iqtreedir option" if (not defined $iqtreedir);
die ("iqtreedir should be only 1 level up\n") if ($iqtreedir =~ /\//);
make_path $iqtreedir unless(-d $iqtreedir);
my $logfile = "$iqtreedir/iqtree.pl.log";

# die "a directory for ASTRAL output should be specified with the --astraldir option" if ($doastral and (not defined $astraldir));
# $astraldir = "astral" if !defined($astraldir);
# make_path $astraldir if !(-d $astraldir);

system("date > $logfile");
system("hostname >> $logfile");

#-----------------------------------------------#
#  get list of input (nexus) files              #
#-----------------------------------------------#

#my @gene = glob("$seqdir/*");
chdir($seqdir) or die "can't go to sequence directory";
my $genefiles = `ls`;
my @genes = split(/\n/, $genefiles);
my @generoots;
print(@genes);
foreach my $gene (@genes){
    my $generoot = $gene;
    $generoot =~ s/\.\w{3}//;
    push @generoots, $generoot;
}
my $nloci = scalar(@genes);
chdir($currentdir) or die "can't go back to original directory";

#-----------------------------------------------#
#  convert nexus to phylip files                #
#  interleaved format: do *not* repeat taxon names
#-----------------------------------------------#

if ($convertphylip) {
  for my $ig (0 .. $#genes){
    my $infn = "$seqdir/${genes[$ig]}";
    print($infn);
    my $oufn = "$phylipdir/${generoots[$ig]}.phy";
    my $read = 0;
    my $removeNames = 0; my $nReadNames = 0;
    my $ntax = 0;
    my $nchar = 0;
    open my $FHi, $infn or die "can't open NEXUS gene sequence file";
    open my $FHo, ">", $oufn or die "can't open PHYLIP gene sequence file";
    while (<$FHi>){
	  if ($read){
	    if (/^\s*;/){ last; } # end of alignment
      if (/^\s*$/){ next; } # don't write blank lines
      if ($removeNames){
          if (/^[^\s]+\s+(.*)/) { print $FHo "$1\n"; }
      } else {                    print $FHo $_;
          $nReadNames++;
          if ($nReadNames == $ntax){ $removeNames = 1; }
      }
	  }
	  if ($read==0){
	    if (/ntax\s*=\s*(\d+)/i){
	    	$ntax = $1;
	    	print $FHo " $ntax ";
	    }
	    if (/nchar\s*=\s*(\d+)/i){
		    $nchar = $1;
		    if ($ntax==0){ print "problem in file $infn: found nchar before ntax\n"}
		    print $FHo "$nchar\n";
	    }
	  }
	  if (/^\s*matrix/i){
	    $read=1;
	    if ($ntax==0 or $nchar==0){
		    print "problem in file $infn: was unable to find ntax ($ntax) or nchar ($nchar)\n";
	    }
	  }
    }
    close $FHi;
    close $FHo;
  }
  $seqdir = $phylipdir;
  for my $ig (0..$#genes){
      $genes[$ig] = $generoots[$ig]. ".phy";
  }
}

#-----------------------------------------------#
#  Run Iqtree  -- No boostrap to increase speed #
#-----------------------------------------------#
# Below codes run iqtree using each individual file. 
# Notice: 
# This code doesn't work for now since we changed the file names for seq-gen outputs 
# open FHlog, ">> $logfile" or die "Cannot open log file $logfile: $!\n";
# chdir($iqtreedir) or die ("can't go to iqtree directory $iqtreedir\n");
# for my $ig (0 .. $#genes){
#     my $infn = "$currentdir/$seqdir/${genes[$ig]}"; 
#     my $oufn = "${generoots[$ig]}";     
#     my $str = "$iqtree -s $infn -m HKY+G -T $numCores -pre $oufn";
#     print FHlog "starting IQ-TREE for gene $ig...\n";
#     print FHlog "$str\n";
#     system($str);
# }
# chdir($currentdir)
# system("date >> $logfile");
# close FHlog;

# Here, use the directory as the input file to speed the process: 
# This code is also robust to different file names. 
open FHlog, ">> $logfile" or die "Cannot open log file $logfile: $!\n";
chdir($iqtreedir) or die ("can't go to iqtree directory $iqtreedir\n");
my $iqtree_input = "$seqdir";
# The below iqtree command is hard-coded 
my $iqtreecmd = "$iqtree -S $iqtree_input -m HKY+G -T $numCores --seed $seed_iqtree -pre gene -B 1000"; 

system($iqtreecmd);
chdir($currentdir);
system("date >> $logfile");
close FHlog;

#-----------------------------------------------#
#  create mapping file                          #
#-----------------------------------------------#
# create a mapping file between input files and output trees

my $treefile   = "$iqtreedir/gene.treefile";   # single file with all trees
my $mappingOUT = "$iqtreedir/mapping.csv";

open(my $mfh, ">", $mappingOUT) or die "Cannot open $mappingOUT: $!\n";
print $mfh "phy_file,tree\n";   # CSV header

# collect all phy files in sorted order
my @phyfiles = sort glob("$seqdir/*.phy");

# open the treefile
open(my $tfh, "<", $treefile) or die "Cannot open $treefile: $!\n";

my $i = 0;
while (my $tree = <$tfh>) {
    chomp $tree;
    if ($i > $#phyfiles) {
        warn "More trees in $treefile than phy files in $seqdir!\n";
        last;
    }
    my $phyfile = $phyfiles[$i];
    print $mfh "$phyfile : \"$tree\"\n";   # quote tree to protect commas
    $i++;
}

close $tfh;
close $mfh;

print "Mapping file written to $mappingOUT\n";

#-----------------------------------------------#
#  restructure output files                     #
#-----------------------------------------------#
# delete intermediate files ending in ckp.gz (checkpoint files) created by IQ-tree: 
my @files_checkpoint = glob("$iqtreedir/*.ckp.gz");
if (@files_checkpoint ) {
    print "Deleting .ckp.gz files in $iqtreedir...\n";
    unlink @files_checkpoint or warn "Failed to delete .ckp.gz files: $!";
    print "Successfully removed the .ckp.gz files.\n";
}

# delete intermediate files ending in .mldist created by IQ-tree: 
my @files_dist = glob("$iqtreedir/*.mldist");
if (@files_dist) {
    print "Deleting .mldist files in $iqtreedir...\n";
    unlink @files_dist or warn "Failed to delete .mldist files: $!";
    print "Successfully removed the .mldist files.\n";
}

# delete intermediate files ending in .bionj created by IQ-tree: 
my @files_bionj = glob("$iqtreedir/*.bionj");
if (@files_bionj) {
    print "Deleting .bionj files in $iqtreedir...\n";
    unlink @files_bionj or warn "Failed to delete .bionj files: $!";
    print "Successfully removed the .bionj files.\n";
}

# delete intermediate files ending in .parstree created by IQ-tree: 
my @files_parstree = glob("$iqtreedir/*.parstree");
if (@files_parstree) {
    print "Deleting .parstree files in $iqtreedir...\n";
    unlink @files_parstree or warn "Failed to delete .parstree files: $!";
    print "Successfully removed the .mldist files.\n";
}

# delete intermediate files ending in .nex created by IQ-tree: 
my @files_nex = glob("$iqtreedir/*.nex");
if (@files_nex) {
    print "Deleting .nex files in $iqtreedir...\n";
    unlink @files_nex or warn "Failed to delete .nex files: $!";
    print "Successfully removed the .nex files.\n";
}

# delete intermediate files ending in .nex created by IQ-tree: 
my @files_pl_log = glob("$iqtreedir/*.pl.log");
if (@files_pl_log) {
    print "Deleting .pl.log files in $iqtreedir...\n";
    unlink @files_pl_log or warn "Failed to delete .pl.log files: $!";
    print "Successfully removed the .pl.log files.\n";
}

# delete intermediate files ending in .nex created by IQ-tree: 
my @files_log = glob("$iqtreedir/*.log");
if (@files_log) {
    print "Deleting .log files in $iqtreedir...\n";
    unlink @files_log or warn "Failed to delete log files: $!";
    print "Successfully removed the log files.\n";
}

# delete intermediate files ending in .nex created by IQ-tree: 
my @files_iqtree = glob("$iqtreedir/*.iqtree");
if (@files_iqtree) {
    print "Deleting .iqtree files in $iqtreedir...\n";
    unlink @files_iqtree or warn "Failed to delete .iqtree files: $!";
    print "Successfully removed the .iqtree files.\n";
}

# Now the output files include .iqtree, .treefile, and .log files 
# create file listing all best trees: one line per gene
my $iqtreeOUT = "$iqtreedir/besttrees.tre";
`cat $iqtreedir/gene*.treefile > $iqtreeOUT`;

exit(0); 

# # ----------------------------------------------#
# #   astral analysis                             #
# # ----------------------------------------------#
# $astraldir = "astral" if !defined($astraldir);
# #my $bsfile =  "$astraldir/BSlistfiles";
# my $astralLOG =  "$astraldir/astral.screenlog";
# my $astralOUT =  "$astraldir/astral.tre";

# # `ls -d $bootpath/* > $bsfile`;

# # Below is for astral.5.7.8.jar
# # my $astralcmd = "$astral -i $iqtreeOUT -b $bsfile -r $numboot -o $astralOUT > $astralLOG 2>&1";

# # Below is for weighted-astral (wastral) and astral-pro (astral-pro3)
# my $astralcmd = "$astral -i $iqtreeOUT -u 1 -o $astralOUT > $astralLOG 2>&1"; 
 
# open FHlog, ">> $logfile";
# if ($doastral){
#     print FHlog "running astral:\n";
# } else {
#     print FHlog "astral could be run with:\n";
# }
# print FHlog "$astralcmd\n";
# close FHlog;
# if ($doastral){
#     system($astralcmd);
# }