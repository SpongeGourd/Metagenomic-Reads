#! /usr/bin/perl -w
## Copyright 2012.05.28 by lina, my first Perl script.
use Bio::SeqIO;

my $usage= "usage: perl cut_genome_randomly.pl input_file_name output_name output_desc frag_length frag_number\n";
#input : input_file_name|output_name.fasta|output_desc.txt|50 or 400|seq_number
my $input = shift @ARGV or die $usage;
my $output_fasta = shift @ARGV or die $usage;
my $output_desc = shift @ARGV or die $usage;
my $len_set = shift @ARGV or die $usage;
my $number = shift @ARGV or die $usage;

#read the Genbank file first and creat the out file
$in = Bio::SeqIO ->new (-file => $input, -format => "Genbank");
$out = Bio::SeqIO ->new (-file => ">$output_fasta", -format => "fasta");
open OUT, ">$output_desc";

while (my $seq = $in ->next_seq() ){
   my $i = 0;

  #if seq is empty, then jump this seq.
	my $jump =0;
	my $seq_content=$seq ->seq;
	if(! $seq_content){
		 $jump++;
		 print "jump $jump!\n";
		 next;	
	} 

	#set the inner cycle times
	while($i<$number){
	#generate two numbers randomly
	my $frag_length = int ($len_set +rand 25);
	my $geno_length = $seq ->length;
	if ($geno_length < ($len_set+25) ){
		last;
	}
	my $start = int (1 + rand ($geno_length-$len_set-25) );
	my $end = $start + $frag_length;
	if ($end > $geno_length){
		$end = $geno_length
		} 

	#get the title of each seq
	my $acc = $seq ->accession_number;

	#get the organism description
	my @desc = $seq ->species->classification;
	my $desc = join "; ", @desc;
	
	#get the cut seq
	$frag = $seq ->subseq($start, $end);

	#remove the "NNN" frag
	if ($frag=~ /NNN/){
	next;
	}

	#write the result into file
	$seq_obj = Bio::Seq ->new (
                           -seq => $frag,
                           -display_id => $acc."_(".$start.",".$end.")",
                           );
	$out -> write_seq($seq_obj);

	#add the seq inf into another file
	print OUT $acc."_(".$start.",".$end.")\n";
	print OUT $desc."\n";
   $i++;
	}
}
close OUT;
