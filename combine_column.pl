#!/usr/bin/perl -w
use strict;

##########################################################################
# Peter Ung @ UMich 
#
# Add selected column(s) from second file to the end of each line of File 1
# File 1 and 2 must have same number of line.
#
# v1.0	07.08.31 
# v2.0	15.03.10 -- added -c comma option

##########################################################################

die "\n    Usage: x.pl <file 1> <file 2> <col no. of F2> ....
           File 1 and 2 must have same number of row.
           Optional: -t  -- use tab instead of space to separate columns.
	             -c  -- use comma instead of space to separate columns.\n\n"
  if @ARGV < 3;

##########################################################################
my $tab   = 0;		# use tab instead of space to separate column
my $comma = 0;		# use comma instead of space to separate

for (my $i = 0; $i <= $#ARGV; $i++) {
  if ($ARGV[$i] =~ /-t/) {
    $tab = 1;
    delete $ARGV[$i];
  }
  elsif($ARGV[$i] =~ /-c/) {
    $comma = 1;
    delete $ARGV[$i];
  }
}

my $name1 = shift @ARGV;
my $name2 = shift @ARGV;

foreach (@ARGV) { die "\n    Error: No column number.\n\n" if /\D+/; }

my @file_1 = &openFile($name1);
my @file_2 = &openFile($name2);

die "$name1 and $name2 do not have same number of line\n" 
  if ($#file_1 != $#file_2);

my @out = ();
for (my $i = 0; $i <= $#file_1; $i++) {
  my $row1 = $file_1[$i];
  my @row1 = @$row1;
  my $row2 = $file_2[$i];
  my @row2 = @$row2;

  for (my $j = 0; $j <= $#ARGV; $j++) {
    push @row1, $row2[$ARGV[$j]-1] if $ARGV[$j] =~ /\d+/;
  }
  push @out, [@row1];
}

foreach my $file (@out) {
  my $count = 0;
  my @file = @$file;
  foreach my $item (@file) {
    print $item;
    print ',' if $comma and $count < $#file;
    print "\t" if $tab and $count < $#file;
    print ' ' if !$tab and !$comma and $count < $#file;
    $count++;
  }
  print "\n";
}

#######################################################
sub openFile {
  my @out= ();
  open FILE, "$_[0]" 
    or die "\n    Error: $_[0] not found; cannot open file.\n\n";
  while (<FILE>) {
    chomp();
    my @row = split;
    push @out, [@row];
  }
  close FILE;
  return @out;
}

