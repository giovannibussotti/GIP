#!/usr/bin/perl -w

use strict;
use warnings;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);


#given an input array containing multiFasta sequences it parse the array and returns an hash containing the header as keys, and the sequences as values
#It doesn't matter if the multifasta contains spatiation withe lines
sub loadMultifastaIntoHash{
  my (@multifasta) = @_;
  my (%allTheSequences , $sequence , $header);
  my $spy = 1;
  foreach my $line (@multifasta){
    chomp $line;
    next if ($line=~/^\s*$/);
    if ($line=~/(>\S+)/){
      if ($spy == 2){
	$spy = 1;
	$allTheSequences{$header} = $sequence;
	$header   = '';
	$sequence = '';
      }
      $header = $1;
      next;
    }
    if (($line!~/>/) and ($line=~/\w+/)){
      $sequence .= $line;
      $spy = 2;
      next;
    }
  }
  $allTheSequences{$header} = $sequence;
  return %allTheSequences;
}


# it parses the mformat=2 blast output lines. I included so far just some fields. The full field description you might wanna include is here:
# http://blast.advbiocomp.com/doc/tabular.html
sub ABblastMformat2parser{
  my ($line) = @_;
  my %currentHit;
  chomp $line;
  my @fields = split (/\s+/,$line);
  $currentHit{'queryName'}   = $fields[0];
  $currentHit{'targetName'}  = $fields[1];
  $currentHit{'evalue'}      = $fields[2];
  $currentHit{'bitscore'}    = $fields[4];
  $currentHit{'score'}       = $fields[5];
  $currentHit{'alignlen'}    = $fields[6];
  $currentHit{'pcident'}     = $fields[10];
  $currentHit{'pcpos'}       = $fields[11];	
  $currentHit{'strand'}      = $fields[19];
  $currentHit{'targetStart'} = $fields[20];
  $currentHit{'targetEnd'}   = $fields[21];

  return %currentHit;
}





#give an array of hash as input to this function. This is a list of X and Y values (that is points coordinates)
#the function will integrate the area
#WARNING: Better using the R function to do it. You have an R script that does it
sub auc {
  my (@points) = (@_);
  my $auc;
  foreach my $p (1..$#points){
    my $y2 = $points[$p]->{'y'};
    my $x2 = $points[$p]->{'x'};
    my $y1 = $points[$p-1]->{'y'};
    my $x1 = $points[$p-1]->{'x'};
    my $y;
    if ($y2 >= $y1){
      $y = $y1;
    }
    else{
      $y = $y2;
    }
    my $rectangle = ($x2 - $x1) * $y;
    my $triangle  = (abs(($x2 - $x1)) * abs(($y2 - $y1))) / 2;
    $auc += $rectangle + $triangle;
    $p++;
  }
  return $auc;
}




#convert chromosome name from roman to arabic numbering
sub chrNameConversion4C_elegans{
  my ($chrToConvert) = @_;
  my %chrConverter = {

'I'     => '1',
'II'    => '2',
'III'   => '3',
'IV'    => '4',
'V'     => '5'
  };
  my $converted = $chrConverter{$chrToConvert};  
  return $converted;
}


#substitute in the list of multifasta, all the space in the header with a given other symbol
sub substituteSpacesOnHeader {
  my ($caracter) = @_;	
  my (@array) = @_;
  foreach my $element (@array) {
	if ($element=~/^>/) {
		$element=~s/\s/$caracter/g;	
	}
  }
  return (@array);
}


#substitute in the list of multifasta, all the headers with a sequencial fake ID
sub substituteHeaderWithID {
  my (@array) = @_;
  my $i	      = 0;
  my %hash;
  foreach my $element (@array) {
	if ($element=~/^>/) {
		chomp ($element);
		push (@{$hash{'ID'}},"$element***ID$i\n");
		$element=">ID$i";
		$i++;
	}
  }
  foreach my $element (@array) {
		push (@{$hash{'FASTA'}},"$element");
  }
  return (%hash);
}


#given a sequence of a certain length it converts it in a sequence random of DNA ###WARNING! THIS FUNCTION DON'T RESPECT THE COMPOSITION!
#THIS JUST PRODUCE A SEQUENCE RANDOM OF THE SAME LENGTH!!!!
#IF YOU WANNA RESPECT THE COMPOSITION YOU MUST MAKE A FUNCTION THAT READ A FASTA AND FOREACH NUCLEOTIDE IT ASIGN A RANDOM NUMBER
#THEN MAKE A SORT OF THESE NUMBERS AND REPRINT THE NUCLEOTIDES. IN THIS WAY THE NUCLEOTIDES WILL BE SHUFFLED PRESERVING THE SAME COMPOSITION!!
sub FASTArandom {
	my ($FASTA) = @_;
	my $size    = length $FASTA; 
	my $i       = 0;
	my $newFASTA;my $old;
	until ($i == $size){
		my $letter;
		my $num = rand(4);	
		if    ($num < 1)    {$letter = 'A'; $newFASTA .= $letter; $i++; next;}
		elsif (($num>1) and ($num < 2)){$letter = 'C'; $newFASTA .= $letter; $i++; next;}
		elsif (($num>2) and ($num<3)){$letter = 'G'; $newFASTA .= $letter; $i++; next;}
		elsif (($num>3) and ($num<4)){$letter = 'T'; $newFASTA .= $letter; $i++; next;}
		elsif ($num == 4){next;}
		else  {die "error, the rand function produced a number different from 0 1 2 3 4...the value is $num old $old:$!\n";}
	}
	return ($newFASTA);
}



## utilizza la funzione di perl printf o sprintf che e' molto meglio!!!!!!!!!!!!!!!!!!!!

#roundOff a number (approximate)
# you must give the number to approximate and the precision, that is the numberof decimal to consider
sub roundOff {
  my ($number)    = @_;
  my ($precision) = @_;	
	if (int ($number) == $number){
	     return ($number);	
	}
	else{
		my $decimal = $number - int ($number);
		my $i    = 0;
		my $half = 0.5;
		while ($i < $precision){
	    		$half .= '0';	
	    		$i++;
		}
		if ($decimal > $half){
			my $value =  int ($number) + 1;
			return ($value);
		}
		if ($decimal < $half){
			my $value =  int ($number);
			return ($value);	
		}
		if ($decimal == $half){
			return ($number);
		}
        }
}





#non funziona
sub roundOff2 {
  my ($number)    = @_;
  my ($precision) = @_;	
  if (int ($number) == $number){
  return ($number);	
  }
  else{
	my $half = 0.5;
	my $decimal = $number - int ($number);
	my @fake = split(//,$decimal);
	my $value       = int($number) .",";
	my $partial     = $precision - 2; 
	my $last        = $precision -1;
	foreach my $i (0..$partial){
		$value .= $fake[$i];
	}

	if ($fake[$precision] > $half){	
		my $lastValue   = $fake[$last] + 1 ;
		$value .= $lastValue;
                return ($value);
	}
	if ($fake[$precision] < $half){
		my $lastValue   = $fake[$last] - 1 ;
		$value .= $lastValue;
                return ($value);
	}
	if ($fake[$precision] == $half){
		my $lastValue   = $fake[$last];
		$value .= $lastValue . '5';
                return ($value);
	}
	
    }
}


		
#convert a string in FASTA format
sub string2FASTA {
    my ($string) = @_;
    my @fake = split (//, "$string");
    my $count = 0;
    my @copy;
    foreach my $letter (@fake){
	$count++;
	push (@copy,$letter);
	if ($count == 60){
	    push (@copy, "\n");
	    $count = 0;
	}
    }
    $string       = join( "", @copy );
    $string =~ s/^\s+//;
    $string =~ s/\s+$//;
    return ($string);
}	 



#it put a certain sequence of an MSA at the bottom
sub orderingMSA {
        my ($clustalFile , $topID) = @_;	
	open (IN, "<$clustalFile") or die "cannot open the un-ordered MSA\n";	
	my @notOrderedAln = <IN>;
	close IN;
	my @order;
	push (@order, $topID);
	my $maxNameLength = length ($topID);
	my %infoSeq;
	my $spy    = 0;
	my $blocks = 0;
	foreach my $line (@notOrderedAln){
		chomp ($line);
		next if ($line=~/CLUSTAL/);
		if ($line=~/^(\S+)\s+(\S+)/){
			my $name = $1;
			my $seq  = $2;
			$maxNameLength = max (length ($name), $maxNameLength);
			$blocks++ if ($name eq $topID); 

			push (@order, $name) unless (($name eq $topID) or ($spy == 2));        
			push (@{$infoSeq{$name}}, $seq);
			if ($spy == 0){$spy = 1;}
		}
		if (($line =~/^\s*$/) and ($spy == 1)){
			$spy = 2;
			next;	
		}
	}
	my $maxSpacing = $maxNameLength + 4;
	open (OUT,">ordered_MSA") or die "cannot create the ordered MSA";
	print OUT "CLUSTAL\n\n"; 

	my $i = 0;
	for ($i = 0; $i <= $blocks-1; $i++){
		foreach my $ID (@order){  
			my $currentLength = length ($ID);
			my $spacesToAdd    = $maxSpacing - $currentLength;
			my $s = 0;
			my $space;
			until ($s == $spacesToAdd){
				$s++;
				$space .= ' ';
			}
			print OUT "$ID" . "$space" . "$infoSeq{$ID}[$i]" . "\n";
		}
		print OUT "\n\n";
	}
}



#given a list and a number this funcion extract a random subset (as specified by the number)
sub randSubset {
my ($subset , @allNameList)=@_;

my $i = 0;
my $size = scalar(@allNameList);
my (%recall , @extractedNumbers , %outList);

#sanity check
if ($size < $subset){
  die "ERROR: the subset $subset specified is bigger than the list size $size. Impossible to extract more elements than the ones submitted! Aborting\n";
}

#to iterate the loop "$subset" times
for ($i=0;$i<$subset;$i++){

  #generate a random number
  my $rNum = rand($size);
  $rNum = sprintf("%.0f", $rNum);
  
  #skipp already extracted numbers and the last number since the foreach below is from 0 to $size-1
  if (($recall{$rNum}) || ($rNum == $size)){
    $i--;
    next;
  }
  $recall{$rNum} = 1;

  push(@extractedNumbers , $rNum);

}

#generating the subset
foreach my $index (@extractedNumbers){
  my $element = $allNameList[$index];
  $outList{$element} = 1;
  #push (@outList,$element)
}
return(%outList)
}



#given a folder containing multiple alignments it produce two folders, one containing the vienna format secondary structures and one the plots
sub bulkRNAalifold {
  my ($dirName) = @_;
  opendir (D, "$dirName") or die "cannot open $dirName $!\n";
  my @allMSA = readdir (D);
  closedir D;

  my $plotDirectory = cwd() . "/plotDir" . $$ . "/";
  mkdir $plotDirectory or die "Cannot create the plotDir $plotDirectory $!\n";
  my $viennaFormatDirectory = cwd() . "/viennaFormatDir" . $$ . "/";
  mkdir $viennaFormatDirectory or die "Cannot create the viennaFormatDirectory $viennaFormatDirectory $!\n";

  foreach my $msa (@allMSA){
    next if (($msa eq '..') or ($msa eq '.') or ($msa=~/~$/));
    my $command = "RNAalifold183 -r -color -noLp -cv 0.6 -nc 0.5 $dirName"."/$msa > $viennaFormatDirectory"."/$msa".".secStruc";
    (system $command) == 0 or die "RNAalifold Error with command $command \n$!\n";
    (system "mv alirna.ps $plotDirectory"."/$msa".".ps") == 0 or die "cannot move the plot file alirna.ps for $msa in $plotDirectory $!\n";
  }
}


#randomize a sequence preserving the same composition
sub fastaRandomKeepingComposition {
  my ($inputSequence) = @_;
  my $size = length ($inputSequence);
  my %hashSequence;
  my @arraySequence   = split(//, $inputSequence);

  #creating an hash that links the symbols to a random number
  foreach my $symbol (@arraySequence){
    my $rNum = rand($size);
    $hashSequence{$rNum} = $symbol;
  }

  #sorting
  my @sortedRandomArraySequence;
  my @sortedRandomNumbers = sort (keys %hashSequence);
  foreach my $sortedRandomNumber (@sortedRandomNumbers){
    push (@sortedRandomArraySequence , $hashSequence{$sortedRandomNumber});

  }
 #  foreach my $sortedRandomSymbol (sort { $hashSequence{$b}->{"rNum"} <=> $hashSequence{$a}->{"rNum"} } keys %hashSequence){
#     push (@sortedRandomArraySequence , $sortedRandomSymbol);
#   }
  my $outputSequence = join("",@sortedRandomArraySequence);
  return ($outputSequence);
}



#shufflle the columns of an alignment (pairwise or multiple)
#it works perfectly, otherwise you can try "t_coffee -other_pg seq_reformat -in Or9a.aln -action +reorder_column random" .. JiaMing told you!
sub suffleAlnColumns {
   my (@sequences) = @_;
   my $size = length ($sequences[0]);

   #creating the hash that link a random number to the alignment positions	
   my %hash;
   foreach my $position (0..$size){
       my $rNum = rand($size);
       $hash{$rNum} = $position;
   }

   #sorting
   my @sortedRandomNumbers = sort (keys %hash);
   my @randomPositions;	
   foreach my $sortedRandomNumber (@sortedRandomNumbers){
       push (@randomPositions , $hash{$sortedRandomNumber});
   }

   #regenerating the sequences with the random column positions
   my @newSequences;
   foreach my $seq (@sequences){
	my $newSeq;
	foreach my $pos (@randomPositions){
            $newSeq .= substr($seq,$pos,1);     
	}
	push(@newSequences,$newSeq);
   }
   return (@newSequences);
}




#given a directory containing chromosome files (one file one chromosome!) it returns an hash containing a for each chromosome header its size
sub takeAssemblyInfo {
  my ($assembly_dir) = @_;
  opendir ( DIR, $assembly_dir ) || die "Error in opening dir $assembly_dir\n";
  my @allFiles = readdir(DIR);
  closedir(DIR);
  my %info;
  foreach my $file (@allFiles){
    next if (($file eq '.') or ($file eq '..'));
    open (F,"<${assembly_dir}/$file") or die "cannot read ${assembly_dir}/$file\n";
    my ($name , $fullSeq);
    foreach my $line (<F>){
      chomp $line;
      next if ($line=~/^\s*$/);
      $line =~ s/ $//g;

      if ($line=~/>(\S+)/){
	$name = $1;
      }
      $fullSeq .= $line;
    }
    close F;
    $info{$name}= length($fullSeq);
  }
  return %info;
}



















#GTF format converter: give an array containing the gtf file and select if you want listed the transcript intervals or just the exon intervals
#the output is:   chr start end strand geneId transcriptId
sub gtfReformat{
  my ($intervalType,@gtfArray)=@_;
  my @out;
  die 'function gtfReformat require as interval type \'exon\' or \'transcript\'\n' if (($intervalType ne 'exon') and ($intervalType ne 'transcript'));

  if ($intervalType eq 'exon'){
    foreach my $line (@gtfArray){
      chomp $line;
      if ($line=~/chr(\S+)\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+\S+\s+(\S+)\s+\S+\s+gene_id\s+\"([^\"]+)\";\s+transcript_id\s+\"([^\"]+)\"/) {
	my $chr          = $1;
	my $start        = min ($2 , $3);
	my $end          = max ($2 , $3);
	my $strand       = $4;
	my $geneId       = $5;
	my $transcriptId = $6;
	#print "chr$chr $start $end $strand $geneId $transcriptId\n";
	push (@out,"chr$chr $start $end $strand $geneId $transcriptId\n");
      }
      else{
	die "parsing error on line\n$line\n";
      }
    }
  }

  if ($intervalType eq 'transcript'){
    my $oldTranscriptId    = 'entering';
    my $oldRealStart = 900000000000000000000000000000000000000000000000000000000000;
    my $oldRealEnd   = 0;
    my ($oldChr , $oldStrand , $oldGeneId);
    foreach my $line (@gtfArray){
      chomp $line;
      if ($line=~/chr(\S+)\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+\S+\s+(\S+)\s+\S+\s+gene_id\s+\"([^\"]+)\";\s+transcript_id\s+\"([^\"]+)\"/) {
	my $chr          = $1;
	my $start        = min ($2 , $3);
	my $end          = max ($2 , $3);
	my $strand       = $4;
	my $geneId       = $5;
	my $transcriptId = $6;

	if (($transcriptId ne $oldTranscriptId) and ($oldTranscriptId ne 'entering')){
	  #print "chr$oldChr $oldRealStart $oldRealEnd $oldStrand $oldGeneId $oldTranscriptId\n";
	  push (@out , "chr$oldChr $oldRealStart $oldRealEnd $oldStrand $oldGeneId $oldTranscriptId\n" );
	  $oldTranscriptId = $transcriptId;
	  $oldGeneId       = $geneId;
	  $oldChr          = $chr;
	  $oldStrand       = $strand;
	  $oldRealStart    = $start;
	  $oldRealEnd      = $end;
	  next;
	}

	my $realStart = min ($start , $oldRealStart);
	my $realEnd   = max ($end , $oldRealEnd);
	$oldTranscriptId = $transcriptId;
	$oldGeneId       = $geneId;
	$oldChr          = $chr;
	$oldStrand       = $strand;
	$oldRealStart    = $realStart;
	$oldRealEnd      = $realEnd;
      }
    }
    #print "chr$oldChr $oldRealStart $oldRealEnd $oldStrand $oldGeneId $oldTranscriptId\n";
    push (@out,"chr$oldChr $oldRealStart $oldRealEnd $oldStrand $oldGeneId $oldTranscriptId\n");
  }
  return @out;
}







#This function accepts 2 array. As you know in the script you must give the array as references.
#The two lists must contain the structure:
#chrTOT start end allYouWant
#The function looks if the intervals of the second list map inside or overlap the intervals of the first list
#It returns an hash of hash of array like this: @{%comparison{contained|overlapping}->{intelval1}}

#the key "all_VS_all" contains all the intervals of ANY list that overlap or is contained into ANY interval

#not overlapping are all the elements of the list2 that are not contained, contains or overlap with any other interval in list 1
#used in the script positiveSelectionOn3019ncRNA.pl
sub compareIntervals{
  my ($referenceIntervalList1 ,$referenceIntervalList2) = @_;
  my @intervalList1 = @{$referenceIntervalList1};
  my @intervalList2 = @{$referenceIntervalList2};
  my (%firstList , %comparison);
  my $sanityCheck1 = 1;
  my $sanityCheck2 = 1;
  #putting the first list in hash
  foreach my $line (@intervalList1){
    chomp $line;
    if ($line =~/^(chr\S+)\s+\S+\s+\S+/){
      my $chr_list1   = $1;
      push (@{$firstList{$chr_list1}},$line);
      $sanityCheck2++;
    }
  }
  die "ERROR!the second list submitted do not contains the format chrTOT start end whateverYouWant\n" if ($sanityCheck2 == 1);

  foreach my $lineList2 (@intervalList2){
    chomp $lineList2;
    if ($lineList2 =~/^(chr\S+)\s+(\S+)\s+(\S+)/){
      my $chr_list2   = $1;
      my $start_list2 = $2;
      my $end_list2   = $3;
      $sanityCheck1++;

      my $nonOverlap_spy = 0;
      foreach my $lineList1 (@{$firstList{$chr_list2}}){
	if ($lineList1 =~/^(chr\S+)\s+(\S+)\s+(\S+)/){
	  my $chr_list1   = $1;
	  my $start_list1 = $2;
	  my $end_list1   = $3;
	  my $tag1 = "$lineList1".'@@1';
	  my $tag2 = "$lineList2".'@@2';
	  #list2 element included in list1 element
	  if (($start_list2 >= $start_list1) and ($end_list2 <= $end_list1)) {
	    push (@{$comparison{'contained'}->{$lineList1}}, $lineList2);

	    push (@{$comparison{'all_VS_all'}->{$tag1}}, $tag2);
	    $nonOverlap_spy = 1;
	    next;
	  }
	  #list2 element left overlap list1 element
	  if (($start_list2 <= $start_list1) and ($end_list2 <= $end_list1) and ($end_list2 >= $start_list1)){
	    push (@{$comparison{'overlapping'}->{$lineList1}}, $lineList2);
	    push (@{$comparison{'all_VS_all'}->{$tag1}}, $tag2);
	    $nonOverlap_spy = 1;
	    next;
	  }
	  #list3 element right overlap list1 element
	  if (($start_list2 >= $start_list1) and ($end_list2 >= $end_list1) and ($start_list2 <= $end_list1)){
	    push (@{$comparison{'overlapping'}->{$lineList1}}, $lineList2);
	    push (@{$comparison{'all_VS_all'}->{$tag1}}, $tag2);
	    $nonOverlap_spy = 1;
	    next;
	  }
	  #list1 element included in list2 element
	  if (($start_list1 >= $start_list2) and ($end_list1 <= $end_list2)) {
	    push (@{$comparison{'all_VS_all'}->{$tag2}}, $tag1);
	    $nonOverlap_spy = 1;
	  }
	}
      }
     push (@{$comparison{'no_overlapping'}->{$lineList2}},'noOverlap') unless ($nonOverlap_spy == 1);

    }
  }
  die "ERROR!the first list submitted do not contains the format chrTOT start end whateverYouWant\n" if ($sanityCheck1 == 1);
  return %comparison;
}





#read a GTF file into an hash. You have to specify the fileName and the label.
#The label can be either gene_id or transcript_id. 
#Accordingly with your choice you will have back an hash like this $infoGtf{$gene_id}->{'keys'}  or like this$infoGtf{$transcript_id}->{'keys'}
#You can also specify a third variable, called skip. By selecting a feature like "gene" or "transcript" you will take just the selected feature
sub readingGTF {
    my ($gtf , $label, $skip ) = @_;
    open (IN,"<$gtf") or die "cannot open the GTFfile $!";
    my %infoGtf;
    foreach my $line (<IN>){
	chomp $line;
	my ($chr , $group , $frame , $strand , $score , $end , $start , $feature , $source , $gene_id , $transcript_id);
	if ($line=~/^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(.*)$/){
	    $chr     = $1;
	    $source  = $2;
	    $feature = $3;
	    $start   = min ($4, $5);
	    $end     = max ($4, $5);
	    $score   = $6;
	    $strand  = $7;
	    $frame   = $8;
	    $group   = $9;
	    if ($group =~/gene_id \"([^\"]+)\"/){
		$gene_id = $1;
	    }
	    else {die "gtf file in wrong format. Impossible to find the gene_id in line:\n$line\n when it is mandatory for this format\n";}
	    if ($group =~/transcript_id \"([^\"]+)\"/){
		$transcript_id = $1;
	    }
	    else {die "gtf file in wrong format. Impossible to find the transcript_id in line:\n$line\n when it is mandatory for this format\n";}
	}
	else {die "cannot read the line $line\n";}

	if (defined $skip){
	    next unless ($skip eq $feature);
	}

	if ($label eq 'transcript_id'){
	    $infoGtf{$transcript_id}->{'chr'}     = $chr;
	    $infoGtf{$transcript_id}->{'source'}  = $source;
	    $infoGtf{$transcript_id}->{'feature'} = $feature;
	    $infoGtf{$transcript_id}->{'start'}   = $start;
	    $infoGtf{$transcript_id}->{'end'}     = $end;
	    $infoGtf{$transcript_id}->{'score'}   = $score;
	    $infoGtf{$transcript_id}->{'strand'}  = $strand;
	    $infoGtf{$transcript_id}->{'frame'}   = $frame;
	    $infoGtf{$transcript_id}->{'gene_id'} = $gene_id;
	}
	elsif ($label eq 'gene_id'){
	    $infoGtf{$gene_id}->{'chr'}           = $chr;
	    $infoGtf{$gene_id}->{'source'}        = $source;
	    $infoGtf{$gene_id}->{'feature'}       = $feature;
	    $infoGtf{$gene_id}->{'start'}         = $start;
	    $infoGtf{$gene_id}->{'end'}           = $end;
	    $infoGtf{$gene_id}->{'score'}         = $score;
	    $infoGtf{$gene_id}->{'strand'}        = $strand;
	    $infoGtf{$gene_id}->{'frame'}         = $frame;
	    $infoGtf{$gene_id}->{'transcript_id'} = $transcript_id;
	}
	else {
	    die "$label is not accepted. use either gene_id or transcript_id in readingGTF function\n";
	}
    }
    close IN;
    return %infoGtf;
}








#read the exons of a gtf file and return an hash containing for each transcript id all the exons informations
sub readingGTFexons {
  my ($GTFfile) = @_;
  open (IN,"<$GTFfile") or die "cannot open the GTFfile $!";
  my %infoGtf;
  foreach my $line (<IN>){
	chomp $line;
	my ($chr , $group , $frame , $strand , $score , $end , $start , $feature , $source , $gene_id , $transcript_id);
	if ($line=~/^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(.*)$/){
	    $chr     = $1;
	    $source  = $2;
	    $feature = $3;
	    $start   = min ($4, $5);
	    $end     = max ($4, $5);
	    $score   = $6;
	    $strand  = $7;
	    $frame   = $8;
	    $group   = $9;
	    if ($group =~/gene_id \"([^\"]+)\"/){
		$gene_id = $1;
	    }
	    else {die "gtf file in wrong format. Impossible to find the gene_id in line:\n$line\n when it is mandatory for this format\n";}
	    if ($group =~/transcript_id \"([^\"]+)\"/){
		$transcript_id = $1;
	    }
	    else {die "gtf file in wrong format. Impossible to find the transcript_id in line:\n$line\n when it is mandatory for this format\n";}

	    next unless ($feature eq 'exon');

	    my %hash = (
			"chr"           => $chr,
			"source"        => $source,
			"feature"       => $feature,
			"start"         => $start,
			"end"           => $end,
			"score"         => $score,
			"strand"        => $strand,
			"frame"         => $frame,
			"group"         => $group,
			"gene_id"       => $gene_id,
		       );
	    push (@{$infoGtf{$transcript_id}}, \%hash);
	}
	else {die "cannot read the line $line\n";}
      }
  close IN;
  return %infoGtf;
}


#read the hash returned by the function readingGTFexons() and returns another hash in which the exons are assembled into transcripts
#the transcript position are correct.
#Bear in mind that the transcripts will mantain for the other fields (like the score, the group..) just the ones related to a single exon.
sub gtfExons2Transcripts{
  my (%infoGtf) = @_;
  my %infoGtfTranscript;
  foreach my $rna (keys %infoGtf){
    my $minExonStart = 1000000000000000000000000000000000000000000000000000000000000000000000000000;
    my $maxExonEnd   = 0;
    my ($chr , $group , $frame , $strand , $score , $feature , $source , $gene_id);
    foreach my $index (0..$#{$infoGtf{$rna}}){
      $minExonStart = min ($minExonStart, $infoGtf{$rna}->[$index]->{"start"});
      $maxExonEnd   = max ($maxExonEnd, $infoGtf{$rna}->[$index]->{"end"});

      $chr           = $infoGtf{$rna}->[$index]->{"chr"};
      $group         = $infoGtf{$rna}->[$index]->{"group"};
      $frame         = $infoGtf{$rna}->[$index]->{"frame"};
      $strand        = $infoGtf{$rna}->[$index]->{"strand"};
      $score         = $infoGtf{$rna}->[$index]->{"score"};
      $feature       = $infoGtf{$rna}->[$index]->{"feature"};
      $source        = $infoGtf{$rna}->[$index]->{"source"};
      $gene_id       = $infoGtf{$rna}->[$index]->{"gene_id"};
    }
    $infoGtfTranscript{$rna}->{'start'}        = $minExonStart;
    $infoGtfTranscript{$rna}->{'end'}          = $maxExonEnd;
    $infoGtfTranscript{$rna}->{'chr'}          = $chr;
    $infoGtfTranscript{$rna}->{'group'}        = $group;
    $infoGtfTranscript{$rna}->{'frame'}        = $frame;
    $infoGtfTranscript{$rna}->{'strand'}       = $strand;
    $infoGtfTranscript{$rna}->{'score'}        = $score;
    $infoGtfTranscript{$rna}->{'feature'}      = $feature;
    $infoGtfTranscript{$rna}->{'source'}       = $source;
    $infoGtfTranscript{$rna}->{'gene_id'}      = $gene_id;
  }
  return %infoGtfTranscript;
}





#read the hash returned by the function readingGTFexons() and returns two hash references, one having the transcripts gtf fields (the same as gtfExons2Transcripts) and one the exons (the same as readingGTFexons) of the bigger transcript for each gene
sub gtfExons2BiggerTranscript {
  my %infoExons          =  @_;
  #TAKING THE TRANSCRIPTS
  my %infoTranscripts   = gtfExons2Transcripts(%infoExons);
  #initializing
  my (%geneInfo , %allGeneIds);
  foreach my $transcript_id (keys %infoExons){
    my $gene_id = $infoTranscripts{$transcript_id}{'gene_id'};
    $allGeneIds{$gene_id} = 1;
    $geneInfo{$gene_id}{'maxSize'} = 0;
    $geneInfo{$gene_id}{'transcript_id'} = 'none';
  }
  #Taking the biggest transcript for each gene id
  foreach my $transcript_id (keys %infoExons){
    my $size    = abs ($infoTranscripts{$transcript_id}{'end'} - $infoTranscripts{$transcript_id}{'start'});
    my $gene_id = $infoTranscripts{$transcript_id}{'gene_id'};
    if ($size > $geneInfo{$gene_id}{'maxSize'}){
      $geneInfo{$gene_id}{'maxSize'}       = $size;
      $geneInfo{$gene_id}{'transcript_id'} = $transcript_id;
    }
  }
  #Output
  my (%biggerTranscript_ex , %biggerTranscript_tx);
  foreach my $gene_id (keys %allGeneIds){
    my $transcript_id = $geneInfo{$gene_id}{'transcript_id'};
    $biggerTranscript_ex{$transcript_id} = \@{$infoExons{$transcript_id}};
    $biggerTranscript_tx{$transcript_id} = \%{$infoTranscripts{$transcript_id}};
  }
  return (\%biggerTranscript_ex , \%biggerTranscript_tx);
}



#read the hash returned by the function readingGTFexons() and returns 3 hash references, one having the gene sizes (genomic loci), one having mature transcript sizes and one having exon sizes
sub gtfExons2sizes {
  my %infoExons          =  @_;
  my (%infoExonsPerGene , %transcriptSizes) ;
  my (@transcriptNames, @geneNames);

  foreach my $transcript (keys %infoExons){
    foreach my $exon (0..$#{$infoExons{$transcript}}){
      my $chr      = $infoExons{$transcript}[$exon]->{'chr'};
      my $start    = $infoExons{$transcript}[$exon]->{'start'};
      my $end      = $infoExons{$transcript}[$exon]->{'end'};
      my $gene     = $infoExons{$transcript}[$exon]->{'gene_id'};
      my $rna      = $infoExons{$transcript}[$exon]->{'transcript_id'};
      my $exonSize = abs ($end - $start);
      $infoExons{$transcript}[$exon]->{'size'} = $exonSize;

      push (@transcriptNames , $transcript) if (! defined $transcriptSizes{$transcript});
      push (@geneNames, $gene)       if (! defined $infoExonsPerGene{$gene});
      $transcriptSizes{$transcript} += $exonSize;

      my %hash = ("chr" =>  $chr, "start" => $start, "end" => $end);
      push (@{$infoExonsPerGene{$gene}}, \%hash);
    }
  }

  #TAKING GENE SIZES
  my %geneSizes;
  foreach my $gene (keys %infoExonsPerGene){
    my $minExonStart = 1000000000000000000000000000000000000000000000000000000000000000000000000000;
    my $maxExonEnd   = 0;
    my $chrKeeper;
    foreach my $index (0..$#{$infoExonsPerGene{$gene}}){
      $minExonStart = min ($minExonStart, $infoExonsPerGene{$gene}->[$index]->{"start"});
      $maxExonEnd   = max ($maxExonEnd, $infoExonsPerGene{$gene}->[$index]->{"end"});
      $chrKeeper    = $infoExonsPerGene{$gene}->[$index]->{"chr"};
    }
    my $geneSize = abs ($maxExonEnd - $minExonStart);
    $geneSizes{$gene} = $geneSize;
  }

  return (\%geneSizes , \%transcriptSizes , \%infoExons);

}






#read the hash returned by the function readingGTFexons() and returns another hash with the same structure containing introns insted than exons
#used in the script: exonGTF_2_intronGTF.pl
sub gtfExons2gtfIntrons {
    my (%infoGtf) = @_;
    my %infoIntrons;
    foreach my $txID (keys %infoGtf){
      @{$infoGtf{$txID}} = sort {$a->{"start"}  <=>  $b->{"start"}} @{$infoGtf{$txID}};
      foreach my $exon (0..$#{$infoGtf{$txID}}){
	next if (! defined $infoGtf{$txID}[1+$exon]);
	my $ex_start = $infoGtf{$txID}[$exon]{'start'};
	my $ex_end   = $infoGtf{$txID}[$exon]{'end'};
	my $next_ex_start = $infoGtf{$txID}[1+$exon]{'start'};

	my %hash = (
			"chr"           => $infoGtf{$txID}[$exon]{'chr'},
			"source"        => $infoGtf{$txID}[$exon]{'source'},
			"feature"       => "intron",
			"start"         => $ex_end +1 ,
			"end"           => $next_ex_start -1,
			"score"         => $infoGtf{$txID}[$exon]{'score'},
			"strand"        => $infoGtf{$txID}[$exon]{'strand'},
			"frame"         => $infoGtf{$txID}[$exon]{'frame'},
			"group"         => $infoGtf{$txID}[$exon]{'group'},
			"gene_id"       => $infoGtf{$txID}[$exon]{'group'}
		       );
	push (@{$infoIntrons{$txID}}, \%hash);
      }
    }
    return %infoIntrons;
}


#read the hash returned by the function readingGTFexons() and returns another hash giving for each transcript ID all the promoter GTF fields
#by defaults it consider the strand of each transcript, it take the start exon (which after the sorting can be the first or the last one depending on the transcript strand) and take the coordinates of the promoter before the start exon
#by default the span of the promoter is 1kb and the spatiation is 0, this means that it is stuck to the start exon (you can edit these values if needed)
#For appropriate promoter sizes see:
#Assigning roles to DNA regulatory motifs using comparative genomics
#RSAT: regulatory sequence analysis tools
#To me 1kb is the safest span
sub gtfExons2gtfPromoters {
    my (%infoGtf) = @_;
    my %infoPromoters;
    foreach my $txID (keys %infoGtf){
      @{$infoGtf{$txID}} = sort {$a->{"start"}  <=>  $b->{"start"}} @{$infoGtf{$txID}};
      my $spatiation = 0;
      my $span       = 1000;
      my ($beginning , $promoterStart , $promoterEnd);
      if ($infoGtf{$txID}[0]{'strand'} eq "+"){
	$beginning = $infoGtf{$txID}[0]{'start'};
	$promoterStart = ($beginning - $spatiation) - $span;
	$promoterEnd   = ($promoterStart + $span) -1;
      }
      elsif ($infoGtf{$txID}[0]{'strand'} eq "-"){
	my $lastExon = scalar(@{$infoGtf{$txID}}) -1;
	$beginning = $infoGtf{$txID}[$lastExon]{'end'};
	$promoterStart = ($beginning + $spatiation) +1;
	$promoterEnd   = ($promoterStart + $span) -1;
      }
      else {die "Error! $infoGtf{$txID}[0]{'strand'} strand is not accepted. Choose either + or -\n";}

      my %hash = (
		  "chr"           => $infoGtf{$txID}[0]{'chr'},
		  "source"        => $infoGtf{$txID}[0]{'source'},
		  "feature"       => "promoter",
		  "start"         => $promoterStart ,
		  "end"           => $promoterEnd ,
		  "score"         => $infoGtf{$txID}[0]{'score'},
		  "strand"        => $infoGtf{$txID}[0]{'strand'},
		  "frame"         => $infoGtf{$txID}[0]{'frame'},
		  "group"         => $infoGtf{$txID}[0]{'group'},
		  "gene_id"       => $infoGtf{$txID}[0]{'group'}
		 );
      $infoPromoters{$txID} = \%hash;
    }
    return %infoPromoters;
}



#read the hash returned by the function readingGTFexons() and returns another hash giving for each transcript ID all the TSS GTF fields
#by defaults it consider the strand of each transcript, it take the start exon (which after the sorting can be the first or the last one depending on the transcript strand) and take the coordinates of the TSS around the start exon
#by default the span of the TSS area is 100
sub gtfExons2gtfTSS {
    my ($span , $ref_infoGtf) = @_;
    my %infoGtf = %{$ref_infoGtf};
    my %infoTSS;
    foreach my $txID (keys %infoGtf){
      @{$infoGtf{$txID}} = sort {$a->{"start"}  <=>  $b->{"start"}} @{$infoGtf{$txID}};
      my ($beginning , $TSSStart , $TSSEnd);
      if ($infoGtf{$txID}[0]{'strand'} eq "+"){
	$beginning = $infoGtf{$txID}[0]{'start'};
      }
      elsif ($infoGtf{$txID}[0]{'strand'} eq "-"){
	my $lastExon = scalar(@{$infoGtf{$txID}}) -1;
	$beginning = $infoGtf{$txID}[$lastExon]{'end'};
      }
      else {die "Error! $infoGtf{$txID}[0]{'strand'} strand is not accepted. Choose either + or -\n";}
      $TSSStart = $beginning -  $span;
      $TSSEnd   = $beginning + $span ;

      my %hash = (
		  "chr"           => $infoGtf{$txID}[0]{'chr'},
		  "source"        => $infoGtf{$txID}[0]{'source'},
		  "feature"       => "TSSarea",
		  "start"         => $TSSStart ,
		  "end"           => $TSSEnd ,
		  "score"         => $infoGtf{$txID}[0]{'score'},
		  "strand"        => $infoGtf{$txID}[0]{'strand'},
		  "frame"         => $infoGtf{$txID}[0]{'frame'},
		  "group"         => $infoGtf{$txID}[0]{'group'},
		  "gene_id"       => $infoGtf{$txID}[0]{'group'}
		 );
      $infoTSS{$txID} = \%hash;
    }
    return %infoTSS;
}




#read the hash returned by the function readingGTFexons() reading the output of overlap (from Sara Djerbaly)
#The normal gtf output include a last field indicating how many overlap there was among two compared gtf, and how much
#it looks like this:
##nb_ov_feat2: 0 list_feat2: .
#nb_ov_feat2: 1 list_feat2: chr7_121968523_121968619_.,
#nb_ov_feat2: 2 list_feat2: chr1_113541552_113542020_.,chr1_113541552_113542060_.,
#the function so far returns a value which indicates the total ammount of overlapping nt and the number of overlapping exons for a given set
sub overlapping_nts{
  my (%infoExons) = @_;
  my (%hash , @exonIdList);
  my $ntOverlap4set   = 0;
  my $exonOverlap4set = 0;
  foreach my $transcript_id (keys %infoExons){
    my $overlap4transcript          = 0;
    my $numberOfOverlappingFeatures = 0;
    my $numberOfOverlappingExons    = 0;
    foreach my $index (0..$#{$infoExons{$transcript_id}}){
      my $exonGroup = $infoExons{$transcript_id}->[$index]->{'group'};
      if ($exonGroup=~/nb_ov_feat2:\s+(\d+)\s+list_feat2:\s+(.*)/){
	my $numberOfFeatures = $1;
	my $rest             = $2;
	if ($numberOfFeatures > 0){
	  while ($rest){
	    if ($rest=~/^.+_([^_]+)_([^_]+)_\.,(.*)$/){
	      my $start = min($1,$2);
	      my $end   = max($1,$2);
	      my $currentOverlappingNucleotides = abs ($end - $start);
	      $rest  = $3;
	      $overlap4transcript += $currentOverlappingNucleotides;
	    }
	    else{
	      $rest='';
	    }
	  }
	  $numberOfOverlappingFeatures += $numberOfFeatures;
	  $numberOfOverlappingExons++;
          push(@exonIdList,$transcript_id);
	}
      }
    }
    $hash{$transcript_id}{'overlappingNucleotides'}      = $overlap4transcript;
    $hash{$transcript_id}{'numberOfOverlappingFeatures'} = $numberOfOverlappingFeatures;
    $hash{$transcript_id}{'numberOfOverlappingExons'}    = $numberOfOverlappingExons;
    $ntOverlap4set   += $overlap4transcript;
    $exonOverlap4set += $numberOfOverlappingExons;
  }
  return ($ntOverlap4set , $exonOverlap4set , \@exonIdList);
}




#it also verify inclusion, as inclusion is a specific kind of overlap
sub verifyOverlap {
  my ($traSta , $geneEnd , $geneSta , $traEnd) = @_;
  my $overlap = 0;
  if ((($traSta <= $geneEnd) && ($traSta >= $geneSta)) || (($traEnd >= $geneSta) && ($traEnd <= $geneEnd)) || (($traSta <= $geneSta) && ($traEnd >= $geneEnd)) ||   (($traSta >= $geneSta) && ($traEnd <= $geneEnd))){
    $overlap = 1;
  }
  return $overlap;
}





sub read_axtNet {
  my ($axtAlignmentFile)= @_;
  my %infoAxt;
  open (AXT,"<$axtAlignmentFile") or die "cannot open the axtAlignmentFile : $!";
  my $spy = 0;
  my $humanChr;
  my $ID;  print "reading the axt file\n";
  foreach my $line (<AXT>){
    next if ($line =~/^#/);

    if (($line =~/^\d+\s+(chr\S+)\s+(\S+)\s+(\S+)\s+(chr\S+)\s+(\S+)\s+(\S+)\s+[+|-]/) and ($spy == 0)) {
      $humanChr     = $1;
      my $humanStart   = min ($2,$3);
      my $humanEnd     = max ($3,$2);
      my $mouseChr     = $4;
      my $mouseStart   = min ($5,$6);
      my $mouseEnd     = max ($6,$5);
      my $mouseStrand  = $7;
      $ID = "$humanChr;"."$humanStart;"."$humanEnd";

      $infoAxt{$humanChr}->{$ID}->{'humanChr'}   = $humanChr;
      $infoAxt{$humanChr}->{$ID}->{'humanStart'} = $humanStart;
      $infoAxt{$humanChr}->{$ID}->{'humanEnd'}   = $humanEnd;
      $infoAxt{$humanChr}->{$ID}->{'mouseChr'}   = $mouseChr;
      $infoAxt{$humanChr}->{$ID}->{'mouseStart'} = $mouseStart;
      $infoAxt{$humanChr}->{$ID}->{'mouseEnd'}   = $mouseEnd;
      $spy = 1;
      next;
    }
    if (($spy == 1) && ($line=~/\S+/)){
      chomp $line;
      $line = trim ($line);
      $infoAxt{$humanChr}->{$ID}->{'humanSequence'}   = $line;
      $spy = 2;
      next;
    }
    if (($spy == 2) && ($line=~/\S+/)){
      chomp $line;
      $line = trim ($line);
      $infoAxt{$humanChr}->{$ID}->{'mouseSequence'}   = $line;
      next;
    }
    if ($line=~/^\s*$/){
      $spy = 0;
      next;
    }
  }
  return %infoAxt;
}
















1;
