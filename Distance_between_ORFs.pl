#! /usr/bin/perl -w
##2012.5.29 by lina.
use Bio::SeqIO;

#type in the input_file_name and the output_file_name
my $input = shift @ARGV;
my $output = shift @ARGV;

#read the input file and creat the output
$in = Bio::SeqIO ->new (-file =>$input, -format => "Genbank");
open OUT, ">$output";
print OUT "Length_of_genome\\Ratio_of_coding_seq\\Average_value_of_distances\\Classification_information\\Distances_between_ORFs\n";

while (my $seq = $in ->next_seq() ){
  #get the length of the genome
	$geno_len = $seq ->length;

	#get the ratio of coding seq
	my $CDS_len = 0;
	for my $feature ($seq->top_SeqFeatures){
		if ( $feature->location->isa('Bio::Location::SplitLocationI')               #CDS   join(5516..5591,7925..8199)的情况
    	           && $feature->primary_tag eq 'CDS' )  {
			for my $location ( $feature->location->sub_Location ) {
				my $each = $location->end - $location->start +1;
				$CDS_len += $each;
			}
		}
		elsif ($feature->primary_tag eq "CDS") {                                    #CDS   4587..5165的情况
			my $each = ($feature->location->end) - ($feature->location->start) +1;
			$CDS_len += $each;
		}
	}
	my $ratio = $CDS_len/$geno_len;


	#get each distance between ORFs
	my @feat = $seq ->top_SeqFeatures;
	my @start;
	my @end;
	for (my $i=0; $i<=$#feat; $i++){
		if ($feat[$i]->location->isa('Bio::Location::SplitLocationI')               #CDS   join(5516..5591,7925..8199)的情况
    	           && $feat[$i]->primary_tag eq 'CDS' ){
			for my $location ( $feat[$i]->location->sub_Location ) {                 #Take join as separate ORFs
				push @start, $location->start;
				push @end, $location->end;
			}
		}
		elsif ($feat[$i]->primary_tag eq "CDS"){                                    #CDS   4587..5165的情况
				push @start, $feat[$i]->location->start;
				push @end, $feat[$i]->location->end;
		}
	}
	my @dis;
	for (my $n=0; $n<$#start; $n++){
		my $distance = $start[$n+1] - $end[$n] -1;
		push @dis, $distance;
	}
	my $dis = join " ", @dis;

	#calculate the average value of distances betw ORFs
	my $dis_sum=0;	
	for (@dis){
		$dis_sum += $_;
	}
	if ($#dis == -1){
	next;
	}
	my $aver = $dis_sum/($#dis+1);

	#get the classification information
	my @classi = $seq ->species->classification;
	@classi = reverse @classi;
	my $classi = join "; ", @classi;

	#output
	print OUT $geno_len."\\".$ratio."\\".$aver."\\".$classi."\\".$dis."\n";
}
close OUT;
