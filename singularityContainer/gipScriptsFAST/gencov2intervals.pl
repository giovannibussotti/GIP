#!/usr/bin/env perl
#Author: Giovanni Bussotti, giovanni.bussotti@pasteur.fr
use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);

my $help = 0;
my ($gcov,$step);
my $result = GetOptions ("gcov=s" => \$gcov,
			 "step=i" => \$step,
                         "help"   => \$help);
$help = 1 unless ($result);
if ((! $help) && (! defined $gcov)) {
  warn "Error! The gcov is missing. Please specify the -gcov parameter\n";
  $help = 1;
}
if ((! $help) && (! defined $step)) {
  warn "Error! The step is missing. Please specify the -step parameter\n";
  $help = 1;
}
pod2usage(-exitval => 0, -verbose => 2, -noperldoc => 1) if ($help);


sub readLine {
	my($row)=@_;
	my ($chr,$pos,$score);
	chomp $row;
	if($row=~/^(\S+)\s+(\S+)\s+(\S+)/){
                $chr   = $1;
                $pos   = $2;
                $score = $3;
        }
	return ($chr,$pos,$score);
}
sub mean {
    	my(@values)=@_;
	my $mean= sum(@values)/scalar(@values);
	return $mean;
}
sub median {
    my @a = sort {$a <=> $b} @_;
    my $length = scalar @a;
    return undef unless $length;
    ($length % 2)
        ? $a[$length/2]
        : ($a[$length/2] + $a[$length/2-1]) / 2.0;
}


my $i = -1;
open(F,"<$gcov") or die "Error. Cannot open $gcov \n$!\n";
my ($chromosome , @positions, @scores);
print "chromosome\tstart\tend\tmean\tmedian\n";
while(<F>){
	#entering
	if ($i == -1) {
		my($chr,$pos,$score)=readLine($_);
		push(@positions,$pos);
		push(@scores,$score);
		$chromosome = $chr;
		$i=1;
	}	

	#looping
	while (($i < $step)&&(my $row = <F>)){
		$i++;
		#my $row = <F>;
		my ($chr,$pos,$score) = readLine($row);
		#chr change in the loop
		if ($chr ne $chromosome){
			print "$chromosome\t" . min(@positions) . "\t" . max(@positions) . "\t" . mean(@scores) . "\t" . median(@scores) . "\n";
        		@positions = ();
        		@scores = ();
			push(@positions,$pos);
                	push(@scores,$score);
			$chromosome = $chr;
        		$i = 1 ;
			next;
		}
		
		push(@positions,$pos);
		push(@scores,$score);
		$chromosome = $chr;	
	}
	
	#printing
	print "$chromosome\t" . min(@positions) . "\t" . max(@positions) . "\t" . mean(@scores) . "\t" . median(@scores) . "\n";
	@positions = ();
	@scores = ();
	$i = -1 ;
}
close F;





=head1 NAME

gencov2intervals.pl - convert the output of bedtools genomecov (run with -d option) to bins and compute the mean and the median for each of them

=head1 SYNOPSIS

gencov2intervals.pl [options]

Options:
-gcov fileName
-step integer

=head1 OPTIONS

=over 8

=item B<-gcov> I<fileName>

output file of bedtools genomecov run with the option -d

=item B<-step> I<integer>

bin size

=back

=head1 DESCRIPTION

To load a bedtools genome coverage track to circos it is needed to compress it into bins and do not report each and every position (otherwise it gets stuck). Even when bedtools is run with option -bga, the file is too big to be read by circos. 

This script reads the output of bedtools genomecov (run with -d option) reporting a score for each nucleotide. Then it loops over each chromosome and reports an average score for each bin. The bins size is defined by the option -step. If the chromosome ends in the middle of a bin, the last bin length will be shorter (truncated to the position of the last nucleotide of the chromosome)   

=cut

