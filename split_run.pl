#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my $script = shift;
my $split = shift;
die "Usage: perl $0 commands.sh split#\n" unless $split;

my $line = `wc -l $script`;
my $line =~ /(\d+).*/;
my $line = $1;
my $chunk_size = sprintf "%d", $line/$split;

use POSIX;
open IN, "$script" or die;
while(<IN>){
  my $current_split = $script.".split".POSIX::ceil($./$chunk_size);
  open OU, ">>$current_split" or die;
  print OU "$_";
  close OU;
}

my @scripts = split /\n/, `ls *split*`;
foreach(@scripts){
  print "nohup sh $_ > $_.o 2> $_.e &\n";
}









