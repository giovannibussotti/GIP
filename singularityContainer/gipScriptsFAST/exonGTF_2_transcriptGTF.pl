#!/usr/bin/perl -w

use strict;
use warnings;
push (@INC,"/bin/");
use Pod::Usage;
use Getopt::Long;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
require "customPerlLib.pl";

my $tag;
my $help                   = 0;
# parsing parameters
my $result = GetOptions ("tag=s"      => \$tag,
                         "help"       => \$help);
$help = 1 unless ($result);
pod2usage(-exitval => 0, -verbose => 2, -noperldoc => 1) if ($help);

$tag = 'transcript_id' unless (defined $tag);


my (@sortedTags , %allTags);
open (GTF,">gtf$$") or die "cannot open gtf$$";
while (<STDIN>){
  chomp $_;
  unless ($_=~/$tag/){system "rm gtf$$";die "ERROR! Line:\n$_\ndo not contain the tag $tag\nPlease set the proper -tag parameter\n";};
  $_=~s/$tag/transcript_id/ unless ($tag eq 'transcript_id');
  $_ .= " gene_id \"FAKE\";" if ($_!~/gene_id\s+\"[^\"]+\"/);

  if ($_=~/transcript_id\s+\"([^\"]+)\"/){
    my $id = $1;
    push (@sortedTags , "$id") unless (defined $allTags{$id});
    $allTags{$id} = 1;
    print GTF "$_\n";
  }
}
close GTF;


my %infoExons       = readingGTFexons("gtf$$");
my %infoTranscripts = gtfExons2Transcripts(%infoExons);
system "rm gtf$$";


foreach my $id (@sortedTags){
  my $outLine;
  $outLine .= $infoTranscripts{$id}->{'chr'}     . "\t";
  $outLine .= $infoTranscripts{$id}->{'source'}  . "\t";
  $outLine .= "transcript" . "\t";
  $outLine .= $infoTranscripts{$id}->{'start'}   . "\t";
  $outLine .= $infoTranscripts{$id}->{'end'}     . "\t";
  $outLine .= $infoTranscripts{$id}->{'score'}   . "\t";
  $outLine .= $infoTranscripts{$id}->{'strand'}  . "\t";
  $outLine .= $infoTranscripts{$id}->{'frame'}   . "\t";
  $outLine .= $infoTranscripts{$id}->{'group'};

  print "$outLine\n";
}



=head1 NAME

exonGTF_2_transcriptGTF.pl - given an exon GTF STDIN it returns the transcripts in GTF format

=head1 SYNOPSIS

perl exonGTF_2_transcriptGTF.pl  [options]

Options:
-tag name
-help 

=head1 OPTIONS

=over 4

=item B<-tag> I<tagName>

default[transcript_id]. Choose the tag to use to indicate which exons belong to what transcripts. Typically you can have tag called like "hitName"

=back

=head1 DESCRIPTION

given GTF STDIN it considers just lines having "exon" has third field
then it assemble the transcripts relying on the transcript_id field. This is the default. If there is any other tag to consider instead like "hitName", use the parameter -tag to specify it.
in any case the output will show the transcript_id field (instead of any specified tag) and it will also include a fake gene_id if this is not found in the input.

=cut
