#!/usr/bin/perl -w
use strict;
use warnings;
push (@INC,"/bin/");
use Pod::Usage;
use Getopt::Long;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
require "customPerlLib.pl";
my $help                   = 0;
my $gtf; 
# parsing parameters
my $result = GetOptions ("gtf=s"  => \$gtf,
                         "help"   => \$help);
$help = 1 unless ($result);

if ((! $help) && (! defined $gtf)) {
  warn "Error! -gtf parameter is missing. Specify an exon gtf file";
  $help = 1;
}
pod2usage(-exitval => 0, -verbose => 2, -noperldoc => 1) if ($help);

my %infoExons         = readingGTFexons("$gtf");
my %infoTranscripts   = gtfExons2Transcripts(%infoExons);

#initializing
my (%geneInfo , %allGeneIds);
foreach my $transcript_id (keys %infoExons){
  my $gene_id = $infoTranscripts{$transcript_id}{'gene_id'};
  $allGeneIds{$gene_id}        = 1;
  $geneInfo{$gene_id}{'end'}   = 0;
  $geneInfo{$gene_id}{'start'} = 10000000000000000000000000000000000000000000000000000000000000000000000;
}

#Taking the smallest start and the biggest end
foreach my $transcript_id (keys %infoExons){
	my $gene_id = $infoTranscripts{$transcript_id}{'gene_id'};	
	my $txstart = $infoTranscripts{$transcript_id}{'start'};
	my $txend   = $infoTranscripts{$transcript_id}{'end'};
	
	$geneInfo{$gene_id}{'chr'}	= $infoTranscripts{$transcript_id}{'chr'};
	$geneInfo{$gene_id}{'source'}	= $infoTranscripts{$transcript_id}{'source'};
	$geneInfo{$gene_id}{'feature'}  = "gene";
	$geneInfo{$gene_id}{'start'}	= min ($geneInfo{$gene_id}{'start'} , $txstart);
	$geneInfo{$gene_id}{'end'}	= max ($geneInfo{$gene_id}{'end'}   , $txend);	
	$geneInfo{$gene_id}{'score'}	= $infoTranscripts{$transcript_id}{'score'};
	$geneInfo{$gene_id}{'strand'}	= $infoTranscripts{$transcript_id}{'strand'};
	$geneInfo{$gene_id}{'frame'}	= $infoTranscripts{$transcript_id}{'frame'};
	$geneInfo{$gene_id}{'gene_id'}	= $infoTranscripts{$transcript_id}{'gene_id'};
	$geneInfo{$gene_id}{'transcript_id'}	= $infoTranscripts{$transcript_id}{'gene_id'};	
}
#printing
printOut();



##########################
sub printOut {
  foreach my $gene_id (keys %allGeneIds){
      	my $chr      = $geneInfo{$gene_id}{'chr'}; 
      	my $source   = $geneInfo{$gene_id}{'source'};
      	my $feature  = "gene";
      	my $start    = $geneInfo{$gene_id}{'start'};
      	my $end      = $geneInfo{$gene_id}{'end'};
      	my $score    = $geneInfo{$gene_id}{'score'};
     	my $strand   = $geneInfo{$gene_id}{'strand'};
      	my $frame    = $geneInfo{$gene_id}{'frame'};
	print "$chr\t$source\t$feature\t$start\t$end\t$score\t$strand\t$frame\tgene_id \"$gene_id\"; transcript_id \"$gene_id\";\n";
  }
}

=head1 NAME
exonGTF_2_geneGTF.pl -  input an exon gtf file and it returns a gene gtf file. The left most exon start is gonna be the gene start. The right most exon end is gonna be the gene end

=head1 SYNOPSIS

perl exonGTF_2_geneGTF.pl -gtf file [-help]

=head1 OPTIONS

=over 4

=item B<-gtf> I<file name>

=back 

=head1 DESCRIPTION

The input exon gtf file must have proper transcript and gene ids.
The gene coordinates are assigned by looking at the smallest exon start and the biggest exon end

=cut
