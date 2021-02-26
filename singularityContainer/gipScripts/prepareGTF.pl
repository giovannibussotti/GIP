#!/usr/bin/env perl
#############################################################################
# giptools                                                                  #
#                                                                           #
# Authors: Giovanni Bussotti                                                #
# Copyright (c) 2021  Institut Pasteur                                      #
#                                                                           #
#                                                                           #
# This program is free software: you can redistribute it and/or modify      #
# it under the terms of the GNU General Public License as                   #
# published by the Free Software Foundation, either version 3 of the        #
# License, or (at your option) any later version.                           #
#                                                                           #
# This program is distributed in the hope that it will be useful,           #
# but WITHOUT ANY WARRANTY; without even the implied warranty of            #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
# GNU General Public License for more details.                              #
#                                                                           #
# You should have received a copy of the GNU General Public License         #
# along with this program.  If not, see <https://www.gnu.org/licenses/>.    #
#                                                                           #
#############################################################################

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use List::Util qw(max min);
use Data::Dumper;

my $help = 0;
my ($gtf,$out);
my $result = GetOptions ("gtf=s"           => \$gtf,
	                     "out=s"           => \$out,
                         "help"            => \$help);
$help = 1 unless ($result);
if ((! $help) && (! defined $gtf)) {
  warn "Error! The input is missing. Please specify the -gtf parameter\n";
  $help = 1;
}
if ((! $help) && (! defined $out)) {
  warn "Error! The input is missing. Please specify the -out parameter\n";
  $help = 1;
}
pod2usage(-exitval => 0, -verbose => 2, -noperldoc => 1) if ($help);

#open differently if compressed
my $gOut = `file $gtf`;
if ($gOut =~/compressed/){
   open( GTF, "gunzip -c $gtf |" ) or die $!;
   } else {
   open( GTF, "<$gtf" ) or die $!;
}
open(OUTGE , ">${out}.ge") or die $!;
open(OUTEX , ">${out}.ex") or die $!;
open(OUTRE , ">${out}.rest") or die $!;


#read
my (%exDerivedFromGe , %geDerivedFromEx);
while(<GTF>){
	if($_=~/^(\S+)\s+\S+\s+gene\s+(\S+)\s+(\S+)\s+.\s+(\S+)\s+.\s+gene_id \"([^\"]+)\"/){
		$exDerivedFromGe{$5}{'chr'}    = $1;
		$exDerivedFromGe{$5}{'start'}  = $2;
		$exDerivedFromGe{$5}{'end'}    = $3;
		$exDerivedFromGe{$5}{'strand'} = $4;
		print OUTGE $_;
	} elsif ($_=~/^(\S+)\s+\S+\s+exon\s+(\S+)\s+(\S+)\s+.\s+(\S+)\s+.\s+gene_id \"([^\"]+)\"/){
		my $chr    = $1;
		my $tmpStart = $2;
		my $tmpEnd   = $3;
		my $strand = $4;
		my $ge     = $5;
		#initialize
		if(! exists $geDerivedFromEx{$ge}){
			$geDerivedFromEx{$5}{'start'}  = "inf";
			$geDerivedFromEx{$5}{'end'}    = 0;
		}
		$geDerivedFromEx{$5}{'chr'}    = $chr;
		$geDerivedFromEx{$5}{'start'}  = min($geDerivedFromEx{$5}{'start'} , $tmpStart);
		$geDerivedFromEx{$5}{'end'}    = max($geDerivedFromEx{$5}{'end'} , $tmpEnd);
		$geDerivedFromEx{$5}{'strand'} = $strand;
		print OUTEX $_;
	} elsif ($_=~/^(\S+)\s+\S+\s+(CDS)/){ #|(start_codon)|(stop_codon)|(3UTR)|(5UTR)
		print OUTRE $_;
	}
}
close GTF;
close OUTGE;
close OUTEX;
close OUTRE;


#complement
open(OUTEX , ">>${out}.ex") or die $!;
foreach my $ge (keys %exDerivedFromGe){
	if (! exists $geDerivedFromEx{$ge}){
		print OUTEX "$exDerivedFromGe{$ge}{chr}\tderivedFromGe\texon\t$exDerivedFromGe{$ge}{start}\t$exDerivedFromGe{$ge}{end}\t.\t$exDerivedFromGe{$ge}{strand}\t.\tgene_id \"$ge\"; transcript_id \"$ge\";\n"
	}
}
close OUTEX;

open(OUTGE , ">>${out}.ge") or die $!;
foreach my $ge (keys %geDerivedFromEx){ 
	if (! exists $exDerivedFromGe{$ge}){
		print OUTGE "$geDerivedFromEx{$ge}{chr}\tderivedFromEx\tgene\t$geDerivedFromEx{$ge}{start}\t$geDerivedFromEx{$ge}{end}\t.\t$geDerivedFromEx{$ge}{strand}\t.\tgene_id \"$ge\"; transcript_id \"$ge\";\n"
	}
}
close OUTGE;

#append
system("cat ${out}.ex ${out}.ge >> ${out}.rest");
system("sort -k1,1 -k4,4n -T _tmp/_sort ${out}.rest > ${out}.ann");
system("rm -rf ${out}.ex ${out}.rest");


=head1 NAME

prepareGTF.pl - extract gene and exons coordinates from GTF

=head1 SYNOPSIS

prepareGTF.pl [options]

Options:
-gtf  fileName
-out  fileName

=head1 OPTIONS

=over 8

=item B<-gtf> I<fileName>

input annotation file in GTF format

=item B<-out> I<fileName>

output name

=back

=head1 DESCRIPTION

Given a GTF file the script separates gene and exons annotations

Then if some genes do not have any annotated exon the script will derive exon coordinates from the gene coordinates

Similarly if some exons do not have any annotated gene the script will derive gene coordinates from the exon coordinates

The output is the gene only annotation file (.ge) and an additional annotation file including genes, exons and if available CDS (.ann file)  

=cut
