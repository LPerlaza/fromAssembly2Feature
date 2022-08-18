#!/usr/bin/perl
use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use File::Basename;
use Cwd;
use File::Spec;
use Data::Dumper;
use File::Copy;

 
### TO DO
#add parameters TRUE/FLASE to skip steps
#check for dependencies before starting
#fork with a set number of cores-DONE
#remove duplication in final table -no necesary
#fusion gff and fasta -post processing
#annotate gff with table as reference -post processing
#agregar coordenadas a los outputs - DONE
#alignear con peptido -DONE
#agregar porcentaje de variabilidad con peptidos -DONE
#agregar posicion de las copias -Done necesito que los procentajes de similaridad correspondan a cada copia diferente
#add a spining wheel
#my $waiting=();
#while( $waiting or 1 ) {
#    spin();
#}

#{
#    my $c=0;  # closure to remember spin state
 #   sub spin {
  #      local $| = 1;
  #      print "\r", qw( | / - \ )[$c++%4];
  #      select undef, undef, undef, 0.25;  # sleep 250 msec
   # }
#}

#el enter del final de los genes todavia esta molestando 
#arreglar stadisticas de la segunda copia
#arreglar classification para que se base en el porcentaje de similaridad con la proteina menos del 60 es truncada



## run like ./fromAssembly2gene.pl -g gene/*fasta -a genome/*

my $start = time;


my @seq=();
my @genes=();
my $out=();
my $R_script=();
my $R_file=();
my $cores;
my $version;
my $help;
my $Kleb;
my $Esch;
my $Ent;
my $percentage;
my $plasmid;
my $species;
my $genes_type;



GetOptions(
    'assemblies|a=s{,}' => \@seq,
    'genes|g=s{,}' => \@genes,
    'out|o=s' =>\$out,
    'cores|c=s' =>\$cores,
    'percentage|i=s' =>\$percentage,
     'plasmid|p' =>\$plasmid,
     'genestype|t=s' =>\$genes_type,
     'version|v' =>\$version,
      'help|h' =>\$help,

    );


if ($help){
&usage();
 exit(1); 
 }

if ($version){
print STDERR "\nVersion 2.21. Program Last updated 22th June 2021 \n\n";
 exit(1); 
 }

if(!$seq[0]){
  &usage();
  die "\nERROR: No assemblies or file empty. Can not continue.\n\n";

exit(1);
}

if(!$genes[0]){
  &usage();
 die "\nERROR: No genes or file empty. Can not continue.\n\n";
exit(1);
}


if (!$out){
&usage();
 print STDERR "\nERROR: No output prefix given.This is obligatory option\n\n";
 exit(1); 
 }


if (!$cores){
$cores=4;
 }
if (!$percentage){
$percentage=80;
 }

if (!$genes_type){
$genes_type ='nucl';
 }

if ($genes_type eq "protein"){
$genes_type='prot';
 }


sub usage{
 print STDERR <<EOF;
 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
																																	
 From Assemblies to genes (Version 22 June 2021)																									
 																																	
Detect genes in an assembly, and get their detailed alignments and description using fromAssembly2gene. 		
It generates descriptive tables of presence, absence and truncated  genes, curated alignments and peptides predictions.													
Dependencies: Blast, prodigal and samtools. R packages "msa", "reshape2", "Biostrings", "seqinr".
the program PlasClass is used to predict pasmids.
																																	
run like: 

./fromAssembly2gene.pl -g gene/*fasta -a genome/* -o GenesInterest -c 10 												


OPTIONS
	--assemblies  -a put all your files (assemblies) in a folder and write here the path with * at the end.							
	--genes 	  -g put all your files (genes in nucleotide sequences) in a folder and write here the path with *fasta at the end. extension ".fasta".	
	--out 		  -o prefix for output folders																						
	--cores		  -c number of cores to use (4 default).
	--percentage  -i number that indicates the minimum percentage identity considered. defaults is 80%. if your genes are not not highly conserved consider lowering this value
	--plasmid	  -p prediction of plasmids. (if not specify it won't predict plasmids)
	--genes_type  -t type of sequence of the genes input (nucleotypes or aminoacid sequences). default is nucl. if you have protein type: protein 
    --help        -h print this help
    --version     -v version

NOTE 1: When you write * at the end of the path the program will take every file in the folder. 
NOTE 2: The genes must have the extension ".fasta" as it is used as a tag for handling file in the program. 
NOTE 3: The blast indexes are going to be created inside the genes folder. Have this in mind when runing several times, use *fasta to avoid using the index files as inputs.	
NOTE 4: Make sure your input genes have a starting and c codon, and are not in reverse. Unexpected mistakes on the identification of the gene can happen when these are not met. We could account for this but the classification of possibly truncated matches is based on the movement of the stop codon in comparison with the reference.
NOTE 5: If you failed to note 4, the program will add a ATG at the start of your sequence and a TAG at the end.
NOTE 6: If you are completely SURE you dont have plasmids in your assemblies you don't need to set the option and the program will assumed all your assembly is chromosome. 																														
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	
EOF
}

  my $rel_path_genes = $genes[0];
  my $abs_path_genes= File::Spec->rel2abs( $rel_path_genes ) ;
  my $dir_genes=dirname($abs_path_genes);
  my $rel_path_genomes = $seq[0];
  my $abs_path_genomes= File::Spec->rel2abs( $rel_path_genomes ) ;
  my $dir_genomes=dirname($abs_path_genomes);

my $cmd_fai="rm $dir_genomes/*fai";
#print $cmd_fai;
system($cmd_fai);


@seq = grep ! /fai$/, @seq;
#print Dumper @seq;



my $output_folder=$out."_output";
my $results_folder=$out."_results";

my $dir= getcwd();
my $diroutput= $dir."/".$output_folder;
my $dirresults= $dir."/".$results_folder;
my $dirresultspep= $dir."/".$results_folder."/Peptides";


#foreach my $files(@genes){ copy($files,$diroriginalgenes)}


if (-d $output_folder) {
die "ERROR: $output_folder folder Exists! Please rename or delete folder called: $output_folder";
}

if (-d $results_folder) {
die "ERROR: $results_folder folder Exists! Please rename or delete folder called: $results_folder";
}

if (! -d  $dir_genes) {
die "ERROR: $dir_genes folder Does NOT Exist! Please direct the program to the right folder path";
}

if (! -d $dir_genomes) {
die "ERROR: $dir_genomes older Does NOT Exist! Please direct the program to the right folder path";
}



system("mkdir $output_folder");
system("mkdir $results_folder");
system("mkdir $dirresultspep");
system("mkdir $results_folder/scripts");

my $diroriginalgenes= $diroutput."/originalGenes";
system("mkdir $diroriginalgenes");

print "
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";

print "folder for intermedia files: $output_folder\n";
print "folder for results files: $results_folder\n";
print "folder for results files: originalGenes\n";


#delete especial characters from windows

if (defined $plasmid ){
	if($plasmid){ $species="SI"};
}else{ 
$species="NO";
}


my $diroriginalgenomes= $diroutput."/originalGenomes";

system("mkdir $diroriginalgenomes");
print "folder for results files: originalGenomes\n";
print "
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";

foreach my $files(@seq){ copy($files,$diroriginalgenomes)}

print "\n\nPLEASE NOTE: Copies of the original files are stored in *output/Originals folders and some format modification are going to be performed on them.\n\n";


     my $out_alignment_description=$dirresults."/out.alignments.description.txt";
	open OUT, ">$out_alignment_description" or die "Cannot open $out_alignment_description for writing\n";
	print OUT "Gene\tRef\tRef_length_nt\tRef_length_AA\tRef_stopCodon\tGenomeFile\tQuery\tStart\tEnd\tStrand\tQuery_length_nt\tQuery_length_AA\tnt_gaps\tstopCodon\tSNPs\tGAPs\tINSERTION\tlocation\tclassification\tPercentageSimilarityAA\tPercentageSimilarityDNA\tN.copies\n";
	close OUT;


my $gene_out_alignment_nt_fasta=();

sub subgene{

my $genes = shift;

unless (-e $genes | -s $genes){ die "ERROR: $genes is empty or doesn't exist";}

my $typefile=`file $genes`;
my $cmd_gene;
my $linegene=$genes."tmp";


#print "\n\n>>>>>>>> $typefile\n\n";

if ( $typefile=~m/CRLF/g){
				my $cmd_gene="awk '/^>/{gsub(\"\\r\",\"\\n\");printf(\"\%s\",\$0);next; }{gsub(\"\\r\",\"\");printf(\"\%s\",\$0);}  END {printf(\"\");} ' < $genes >$linegene";
						system($cmd_gene);
						#print $cmd_gene;
					}elsif( $typefile=~m/CR/g){
					    my $cmd_gene="awk '/^>/{gsub(\"\\r\",\"\\n\");printf(\"\%s\",\$0);next;}{gsub(\"\\\$\",\"\");printf(\"\%s\",\$0);}  END {printf(\"\");} ' < $genes | awk '/^>/{printf(\"\%s\\n\",\$0);next; }{printf(\"\%s\",\$0);}  END {printf(\"\");}'  >$linegene";
						system($cmd_gene);
						#print $cmd_gene;
						}else{ 
					    my $cmd_gene="awk '/^>/{printf(\"\%s\\n\",\$0);next; }{printf(\"\%s\",\$0);}  END {printf(\"\");} ' < $genes >$linegene";
						system($cmd_gene);
						#print $cmd_gene;
						}
						
	
	my $lastcodon=`grep -o '..[^ ]\$' $linegene|tail -1`;
	chomp($lastcodon);
   	$lastcodon =~ s/(.*?)\s?[\n?|\r?]/$1/;
   #	print "\n\n>>$lastcodon<<";
   	
   	my $dec;
   	my @array=("TAG","TAA","TGA","tag","taa","tga");
   	
   	
   	my %params = map { $_ => 1 } @array;


	if(exists($params{$lastcodon}) | $genes_type eq 'prot') { system("echo '\n' >> $linegene")}elsif ( $genes_type eq 'nucl'){system("echo 'TAG\n' >> $linegene")} 

   	my $gene_tag=basename($genes);
   	$gene_tag=~s/.fasta$//g;
 		
	$gene_out_alignment_nt_fasta=$dirresults."/".$gene_tag.".nt_alignment.fasta";
	open OUT1, ">$gene_out_alignment_nt_fasta" or die "Cannot open $gene_out_alignment_nt_fasta for writing\n";
	#sed 's/.*\(...\)/\1/'

	system("cat $linegene  >$gene_out_alignment_nt_fasta.tmp");
	close OUT1;
	

	my $genenew=basename($genes);
	system("mv $linegene $diroutput/originalGenes/$genenew");

	print "\t\tRunning makeblastdb on $genes\n";
	my $gene_location=$diroutput."/originalGenes/".$genenew;
 	my $gene_name=$gene_location;
	$gene_name=~s/.fasta//g;
 	my $cmd_makeblastdb= "makeblastdb -in ".$gene_location." -dbtype ".$genes_type." -out ".$gene_name."\n";
  	print "You are asking to run\n $cmd_makeblastdb\n";
  
  	system($cmd_makeblastdb);
  	

 
}




sub sub0{
my $seq = shift;
my $genome=basename($seq);
my $dirproplas=$diroutput."/".$genome.".plasmid";
my $cmd_prod_plasmid="mkdir $dirproplas";
system("$cmd_prod_plasmid");	
my $dirprochr=$diroutput."/".$genome.".chromosome";
my $cmd_prod_chr="mkdir $dirprochr";
system("$cmd_prod_chr");

my $plaspred=$diroutput."/".$genome."_plas.prob.out";


#print "Rscript plasmid_chromosome_Detection.r $seq  $dirprochr $dirproplas \"$species\" --slave\n";
		system("classify_fasta.py -f  $seq  -o $plaspred");
      	

 }
 

sub sub1{
my $seq = shift;
my $genome=basename($seq);
my $dirproplas=$diroutput."/".$genome.".plasmid";
my $cmd_prod_plasmid="mkdir $dirproplas";
my $dirprochr=$diroutput."/".$genome.".chromosome";
my $cmd_prod_chr="mkdir $dirprochr";

my $plaspred=$diroutput."/".$genome."_plas.prob.out";
      	
      	system("Rscript $results_folder/scripts/plasmid_chromosome_Detection.r $plaspred $dirprochr $dirproplas --slave");
 		#print "done with $seq";
 }
 
 
 
 
 sub sub2{
my $seq = shift;
my $genome=basename($seq);

my $dirproplas=$diroutput."/".$genome.".plasmid";	
my $dirprochr=$diroutput."/".$genome.".chromosome";
my $input_plasmid=$dirproplas."/".$genome.".plasmid_contigslist.txt";
my $input_chr=$dirprochr."/".$genome.".chromosome_contigslist.txt";
my $output_plasmid=$input_plasmid;
my $output_chr=$input_chr;
	
	$output_plasmid=~s/_contigslist.txt/.fasta/g;
	$output_chr=~s/_contigslist.txt/.fasta/g;

	#print "      processing    $output_plasmid \n";	
	#print "      processing    $output_chr \n";

if(-s $input_plasmid){
#print "enter Plasmid \n";
	my $temp_plasmid=$output_plasmid.".tmp";
	my $temp_plasmid_1=$output_plasmid.".tmp1";
    my $cmd_cut="cut -d \" \" -f1 $input_plasmid >$temp_plasmid";
    system($cmd_cut);
	system("xargs samtools faidx $seq< $temp_plasmid > $temp_plasmid_1");
	my $cmd_awk="awk '/^>/{gsub(\" \",\"_\"); \$0=\$0\"\\:$genome"."_plasmid\"}1' $temp_plasmid_1 >$output_plasmid";
	#my $cmd_awk="awk '/^>/{\$0=\$0\"\\:$genome"."_plasm\"}1' $temp_plasmid_1 >$output_plasmid";
	#print ">>>>>>$cmd_awk\n\n";
	system($cmd_awk);
	system( "rm $temp_plasmid");
	system( "rm $temp_plasmid_1");
	}
	
if( -s $input_chr){
#print "enter chromosome \n";
	my $temp_chr=$output_chr.".tmp";
	my $temp_chr_1=$output_chr.".tmp1";
    my $cmd_cut1="cut -d \" \" -f1 $input_chr >$temp_chr";
	system($cmd_cut1);
	system("xargs samtools faidx $seq< $temp_chr > $temp_chr_1");
	my $cmd_awk1="awk '/^>/{gsub(\" \",\"_\"); \$0=\$0\"\\:$genome"."_chr\"}1' $temp_chr_1 >$output_chr";
	#my $cmd_awk1="awk '/^>/{\$0=\$0\"\\:$genome"."_chr\"}1' $temp_chr_1 >$output_chr";	
#	print ">>>>>$cmd_awk1\n\n";
	system($cmd_awk1);
	system( "rm $temp_chr");
	system( "rm $temp_chr_1");

	}
	

}


sub sub3{
my $seq = shift;
my $genome=basename($seq);
my $dirproplas=$diroutput."/".$genome.".plasmid";
my $dirprochr=$diroutput."/".$genome.".chromosome";
my  $input_plasmid=$dirproplas."/".$genome.".plasmid_contigslist.txt";
my  $input_chr=$dirprochr."/".$genome.".chromosome_contigslist.txt";
my $output_plasmid=$input_plasmid;
my $output_chr=$input_chr;
	$output_plasmid=~s/_contigslist.txt/.fasta/g;
	$output_chr=~s/_contigslist.txt/.fasta/g;

my $name_output_plasmid=$output_plasmid;
 $name_output_plasmid=~s/.fasta$//;

my $name_output_chr=$output_chr;
 $name_output_chr=~s/.fasta$//;

if(-s $output_plasmid){

my $genes_gff= $dirproplas."/".$genome.".plasmid.genes.gff";
my $genes_gfftmp= $dirproplas."/".$genome.".plasmid.genes.gfftmp";
my $genes_faa= $dirproplas."/".$genome.".plasmid.genes.faa";
my $genes_fasta= $dirproplas."/".$genome.".plasmid.genes.fasta";
my $temp=$genes_fasta.".tmp";

#print "\n\n\nprodigal -i $output_plasmid -o $gene_gbk -a $genes_faa -d $genes_fasta\n\n\n";
system("prodigal -i $output_plasmid -o $genes_gfftmp -a $genes_faa -d $genes_fasta -f gff -q -p meta");
system("\(cat $genes_gfftmp; echo '##FASTA'; cat $output_plasmid\) >$genes_gff");              
system("rm $genes_gfftmp");       
}

if(-s $output_chr ){

my $genes_gff1=$dirprochr."/".$genome.".chr.genes.gff";
my $genes_gff1tmp=$dirprochr."/".$genome.".chr.genes.gfftmp";
my $genes_faa1=$dirprochr."/".$genome.".chr.genes.faa";
my $genes_fasta1=$dirprochr."/".$genome.".chr.genes.fasta";
my $temp1=$genes_fasta1.".tmp";

#print "\n\n\nprodigal -i $output_chr -o $gene_gbk1 -a $genes_faa1 -d $genes_fasta1\n\n\n";
#system("prodigal -i $output_chr -o $gene_gbk1 -a $genes_faa1 -d $genes_fasta1 -p meta -q");
system("prodigal -i $output_chr -o $genes_gff1tmp -a $genes_faa1 -d $genes_fasta1 -f gff -q  -p meta");
system("\(cat $genes_gff1tmp; echo '##FASTA'; cat  $output_chr\) >$genes_gff1");              
system("rm $genes_gff1tmp");    

	}
}



sub sub4{
my $seq = shift;
my $genome=basename($seq);
my $dirproplas=$diroutput."/".$genome.".plasmid";
my $dirprochr=$diroutput."/".$genome.".chromosome";

 my  $input_plasmid=$dirproplas."/".$genome.".plasmid_contigslist.txt";
 my  $input_chr=$dirprochr."/".$genome.".chromosome_contigslist.txt";
 

	my $output_plasmid=$input_plasmid;
    my $output_chr=$input_chr;
	
	$output_plasmid=~s/_contigslist.txt/.fasta/g;
	$output_chr=~s/_contigslist.txt/.fasta/g;


my $name_output_plasmid=$output_plasmid;
 $name_output_plasmid=~s/.fasta$//;

my $name_output_chr=$output_chr;
 $name_output_chr=~s/.fasta$//;


my $genes_gff= $dirproplas."/".$genome.".plasmid.genes.gff";
my $genes_faa= $dirproplas."/".$genome.".plasmid.genes.faa";
my $genes_fasta= $dirproplas."/".$genome.".plasmid.genes.fasta";
my $temp=$genes_fasta.".tmp";
my $temp_1=$genes_faa.".tmp";

if(-s $genes_fasta){
system("cp $genes_fasta $temp");
my $cmd="awk '/^>/ {printf(\"\\n%s\\n\",\$0);next; } { printf(\"\%s\",\$0);}  END {printf(\"\\n\");}' < $temp  > $genes_fasta";
system( $cmd);
system("sed 's/ # /###/g' $genes_fasta >$temp");
system("mv $temp $genes_fasta");

system("cp $genes_faa $temp_1");
my $cmd_prt="awk '/^>/ {printf(\"\\n%s\\n\",\$0);next; } { printf(\"\%s\",\$0);}  END {printf(\"\\n\");}' < $temp_1  > $genes_faa";
system( $cmd_prt);
system("sed 's/ # /###/g' $genes_faa >$temp_1");
system("mv $temp_1 $genes_faa");


}

my $genes_gff1=$dirprochr."/".$genome.".chr.genes.gff";
my $genes_faa1=$dirprochr."/".$genome.".chr.genes.faa";
my $genes_fasta1=$dirprochr."/".$genome.".chr.genes.fasta";
my $temp1=$genes_fasta1.".tmp";
my $temp1_1=$genes_faa1.".tmp";


if(-s $genes_fasta1){
system("cp $genes_fasta1 $temp1");
my $cmd1="awk '/^>/ {printf(\"\\n%s\\n\",\$0);next; } { printf(\"\%s\",\$0);}  END {printf(\"\\n\");}' < $temp1  > $genes_fasta1";
#print $cmd1."\n";
system( $cmd1);
system("sed 's/ # /###/g' $genes_fasta1 >$temp1");
system("mv $temp1 $genes_fasta1");


system("cp $genes_faa1 $temp1_1");
my $cmd_prt1="awk '/^>/ {printf(\"\\n%s\\n\",\$0);next; } { printf(\"\%s\",\$0);}  END {printf(\"\\n\");}' < $temp1_1  > $genes_faa1";
system( $cmd_prt1);
system("sed 's/ # /###/g' $genes_faa1 >$temp1_1");
system("mv $temp1_1 $genes_faa1");
}
}


sub sub5{
my $seq = shift;
my $genome=basename($seq);
my $dirproplas=$diroutput."/".$genome.".plasmid";	
my $dirprochr=$diroutput."/".$genome.".chromosome";
my $input_plasmid=$dirproplas."/".$genome.".plasmid_contigslist.txt";
my $input_chr=$dirprochr."/".$genome.".chromosome_contigslist.txt";
my $output_plasmid=$input_plasmid;
my $output_chr=$input_chr;
	$output_plasmid=~s/_contigslist.txt/.fasta/g;
	$output_chr=~s/_contigslist.txt/.fasta/g;
	 

my $name_output_plasmid=$output_plasmid;
 $name_output_plasmid=~s/.fasta$//;

my $name_output_chr=$output_chr;
 $name_output_chr=~s/.fasta$//;


my $genes_gff= $dirproplas."/".$genome.".plasmid.genes.gff";
my $genes_faa= $dirproplas."/".$genome.".plasmid.genes.faa";
my $genes_fasta= $dirproplas."/".$genome.".plasmid.genes.fasta";
my $temp=$genes_fasta.".tmp";

my $genes_gff1=$dirprochr."/".$genome.".chr.genes.gff";
my $genes_faa1=$dirprochr."/".$genome.".chr.genes.faa";
my $genes_fasta1=$dirprochr."/".$genome.".chr.genes.fasta";
my $temp1=$genes_fasta1.".tmp";

	foreach my $gene(@genes){
	
	my $genename=basename($gene);
	my $gene_id=$diroriginalgenes."/".$genename;
	$gene_id=~s/.fasta$//g;
	
	my $gene_tag=basename($gene);
	$gene_tag=~s/.fasta$//g;
	
	
	my $outputblast= $name_output_plasmid.".".$gene_tag.".Blast.txt";
	my $outputblast1= $name_output_chr.".".$gene_tag.".Blast.txt";
	
	my $tempblast=$outputblast.".tmp";
	my $tempblast1=$outputblast1.".tmp";
	
if ($genes_type eq"nucl"){	
	if(-s $genes_fasta){
	my $cmd_blast="blastn -query $genes_fasta -db $gene_id -dust no -outfmt \"6 qseqid sseqid pident qlen slen length qstart qend sstart send mismatch gapope evalue bitscore\" -perc_identity $percentage -out $outputblast";

#print "$cmd_blast\n";
	  system ($cmd_blast);

		
	  }
	  
	if(-s $genes_fasta1){
	my $cmd_blast1="blastn -query $genes_fasta1 -db $gene_id -dust no -outfmt \"6 qseqid sseqid pident qlen slen length qstart qend sstart send mismatch gapope evalue bitscore\" -perc_identity $percentage -out $outputblast1";
#print "$cmd_blast1\n";
	system ($cmd_blast1);
	
	}
}elsif ($genes_type eq "prot"){

	if(-s $genes_fasta){
	my $cmd_blast="blastp -query $genes_faa -db $gene_id  -outfmt \"6 qseqid sseqid pident qlen slen length qstart qend sstart send evalue bitscore qcovs qcovhsp scovhsp\"  -out $outputblast";
#print "$cmd_blast\n";
	  system ($cmd_blast);
	  
	  my $cmd_highper="awk '\$3>=".$percentage." && \$13>=".$percentage." {print \$0}' < $outputblast >$tempblast";
		system ($cmd_highper);
		system (" mv $tempblast $outputblast")
		
	  }
	  
	if(-s $genes_fasta1){
	my $cmd_blast1="blastp -query $genes_faa1 -db $gene_id -outfmt \"6 qseqid sseqid pident qlen slen length qstart qend sstart send evalue bitscore qcovs qcovhsp scovhsp\"  -out $outputblast1";
#print "$cmd_blast1\n";
	system ($cmd_blast1);
	my $cmd_highper="awk '\$3>=".$percentage." && \$13>=".$percentage." {print \$0}' < $outputblast1 >$tempblast1";
		system ($cmd_highper);
		system (" mv $tempblast1 $outputblast1");
			}
		}
	}
}



sub sub5_1{
my $seq = shift;
my $genome=basename($seq);
my $dirproplas=$diroutput."/".$genome.".plasmid";	
my $dirprochr=$diroutput."/".$genome.".chromosome";
my $input_plasmid=$dirproplas."/".$genome.".plasmid_contigslist.txt";
my $input_chr=$dirprochr."/".$genome.".chromosome_contigslist.txt";
my $output_plasmid=$input_plasmid;
my $output_chr=$input_chr;
	$output_plasmid=~s/_contigslist.txt/.fasta/g;
	$output_chr=~s/_contigslist.txt/.fasta/g;
	 

my $name_output_plasmid=$output_plasmid;
 $name_output_plasmid=~s/.fasta$//;

my $name_output_chr=$output_chr;
 $name_output_chr=~s/.fasta$//;


my $genes_gff= $dirproplas."/".$genome.".plasmid.genes.gff";
my $genes_faa= $dirproplas."/".$genome.".plasmid.genes.faa";
my $genes_fasta= $dirproplas."/".$genome.".plasmid.genes.fasta";
my $temp=$genes_fasta.".tmp";

my $genes_gff1=$dirprochr."/".$genome.".chr.genes.gff";
my $genes_faa1=$dirprochr."/".$genome.".chr.genes.faa";
my $genes_fasta1=$dirprochr."/".$genome.".chr.genes.fasta";
my $temp1=$genes_fasta1.".tmp";

	foreach my $gene(@genes){
	
	my $genename=basename($gene);
	my $gene_id=$diroriginalgenes."/".$genename;
	$gene_id=~s/.fasta$//g;
	
	my $gene_tag=basename($gene);
	$gene_tag=~s/.fasta$//g;
	
	
	my $outputblast= $name_output_plasmid.".".$gene_tag.".Blast.txt";
	my $outputblast1= $name_output_chr.".".$gene_tag.".Blast.txt";
	
	my $tempblast=$outputblast.".tmp";
	my $tempblast1=$outputblast1.".tmp";
	
if ($genes_type eq"nucl"){	
	if(-s $genes_fasta){

	  my $cmd_highper="awk '\$3>=".$percentage." && \$13>=".$percentage." {print \$0}' < $outputblast >$tempblast";
		system ($cmd_highper);
		system (" mv $tempblast $outputblast")
		
	  }
	  
	if(-s $genes_fasta1){

	my $cmd_highper="awk '\$3>=".$percentage." && \$13>=".$percentage." {print \$0}' < $outputblast1 >$tempblast1";
		system ($cmd_highper);
		system (" mv $tempblast1 $outputblast1");
		
	}
}elsif ($genes_type eq "prot"){

	if(-s $genes_fasta){
	  
	  my $cmd_highper="awk '\$3>=".$percentage." && \$13>=".$percentage." {print \$0}' < $outputblast >$tempblast";
		system ($cmd_highper);
		system (" mv $tempblast $outputblast")
		
	  }
	  
	if(-s $genes_fasta1){

	my $cmd_highper="awk '\$3>=".$percentage." && \$13>=".$percentage." {print \$0}' < $outputblast1 >$tempblast1";
		system ($cmd_highper);
		system (" mv $tempblast1 $outputblast1");
		
			}
		}
	}
}




sub sub6{
my $seq = shift;
my $genome=basename($seq);
my $dirproplas=$diroutput."/".$genome.".plasmid";	
my $dirprochr=$diroutput."/".$genome.".chromosome";
my $input_plasmid=$dirproplas."/".$genome.".plasmid_contigslist.txt";
my $input_chr=$dirprochr."/".$genome.".chromosome_contigslist.txt";
my $output_plasmid=$input_plasmid;
my $output_chr=$input_chr;
	$output_plasmid=~s/_contigslist.txt/.fasta/g;
	$output_chr=~s/_contigslist.txt/.fasta/g;
	 

my $name_output_plasmid=$output_plasmid;
 $name_output_plasmid=~s/.fasta$//;

my $name_output_chr=$output_chr;
 $name_output_chr=~s/.fasta$//;


my $genes_gff= $dirproplas."/".$genome.".plasmid.genes.gff";
my $genes_faa= $dirproplas."/".$genome.".plasmid.genes.faa";
my $genes_fasta= $dirproplas."/".$genome.".plasmid.genes.fasta";
my $temp=$genes_fasta.".tmp";

my $genes_gff1=$dirprochr."/".$genome.".chr.genes.gff";
my $genes_faa1=$dirprochr."/".$genome.".chr.genes.faa";
my $genes_fasta1=$dirprochr."/".$genome.".chr.genes.fasta";
my $temp1=$genes_fasta1.".tmp";

	foreach my $gene(@genes){
	
	my $genename=basename($gene);
	my $gene_id=$diroriginalgenes."/".$genename;
	$gene_id=~s/.fasta$//g;
	
	my $gene_tag=basename($gene);
	$gene_tag=~s/.fasta$//g;
	

	my $outputblast= $name_output_plasmid.".".$gene_tag.".Blast.txt";
	my $outputblast1= $name_output_chr.".".$gene_tag.".Blast.txt";
	
	my $gene_ref=$gene_id.".fasta";

if(-s $outputblast){
 my $temp_plasmid=$name_output_plasmid.".".$gene_tag.".tmp";
	system("cut -f 1 $outputblast >$temp_plasmid");
	my $output_gene= $name_output_plasmid.".".$gene_tag.".fasta";


if( $genes_type eq 'nucl'){
	system("xargs samtools faidx $genes_fasta< $temp_plasmid >$output_gene");
}elsif($genes_type eq 'prot'){
system("xargs samtools faidx $genes_faa< $temp_plasmid >$output_gene");}

	system("rm $temp_plasmid");
	my $output_ref_plasmid=$output_gene.".plusRef.fasta";
	system("cat $gene_ref $output_gene >$output_ref_plasmid");

my $out_alignment_description=$dirresults."/".$genome.".".$gene_tag.".plasmid.out.alignments.description";
	open OUT, ">$out_alignment_description" or die "Cannot open $out_alignment_description for writing\n";

my $gene_genome_plas_out_alignment_AA_fasta=$dirresultspep."/".$gene_tag."_".$genome.".plasm.AAsequence.fasta";
		open OUT1, ">$gene_genome_plas_out_alignment_AA_fasta" or die "Cannot open $gene_genome_plas_out_alignment_AA_fasta for writing\n";

my $gene_genome_plas_out_alignment_nt_fasta=$dirresults."/".$gene_tag."_".$genome.".plasm.nt_alignment.fasta.tmp";
		open OUT2, ">$gene_genome_plas_out_alignment_nt_fasta" or die "Cannot open  $gene_genome_plas_out_alignment_nt_fasta for writing\n";

	
		#print "Rscript curatingAlignments.r $output_ref_plasmid $gene_tag $gene_genome_plas_out_alignment_nt_fasta $gene_genome_plas_out_alignment_AA_fasta $out_alignment_description --slave\n\n";

  		system(" Rscript $results_folder/scripts/curatingAlignments.r $output_ref_plasmid $gene_tag $gene_genome_plas_out_alignment_nt_fasta $gene_genome_plas_out_alignment_AA_fasta $out_alignment_description --slave");

}else { 
 my $out_alignment_description=$dirresults."/".$genome.".".$gene_tag.".plasmid.out.alignments.description";
	open OUT, ">$out_alignment_description" or die "Cannot open $out_alignment_description for writing\n";

	my $firstline_Ref=`grep ">" $gene_id.fasta`;
	chomp($firstline_Ref);
	$firstline_Ref=~s/[\n\r]//g;
	$firstline_Ref=~s/>//g;
	print OUT "$genename\t\"$firstline_Ref\"\tNA\tNA\tNA\t\"$genome\"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tplasmid\tabsent\t0\t0\t0\n";
	close OUT;
 }

if(-s $outputblast1){
my $temp_chr=$name_output_chr.".".$gene_tag.".tmp";
	system("cut -f 1 $outputblast1 >$temp_chr");
	my $output_gene1= $name_output_chr.".".$gene_tag.".fasta";
	
if( $genes_type eq 'nucl'){
	system("xargs samtools faidx $genes_fasta1< $temp_chr >$output_gene1");
}elsif($genes_type eq 'prot'){
system("xargs samtools faidx $genes_faa1< $temp_chr >$output_gene1");}
	
	system("rm $temp_chr");
	my $output_ref_chr=$output_gene1.".plusRef.fasta";
	system("cat $gene_ref $output_gene1 >$output_ref_chr");

my $out_alignment_description=$dirresults."/".$genome.".".$gene_tag.".chr.out.alignments.description";
	open OUT, ">$out_alignment_description" or die "Cannot open $out_alignment_description for writing\n";

my $gene_genome_chr_out_alignment_AA_fasta=$dirresultspep."/".$gene_tag."_".$genome.".chr.AAsequence.fasta";
		open OUT1, ">$gene_genome_chr_out_alignment_AA_fasta" or die "Cannot open $gene_genome_chr_out_alignment_AA_fasta for writing\n";

my $gene_genome_chr_out_alignment_nt_fasta=$dirresults."/".$gene_tag."_".$genome.".chr.nt_alignment.fasta.tmp";
		open OUT2, ">$gene_genome_chr_out_alignment_nt_fasta" or die "Cannot open  $gene_genome_chr_out_alignment_nt_fasta for writing\n";
	  	
	  	#print " Rscript curatingAlignments.r  $output_ref_chr $gene_tag $gene_genome_chr_out_alignment_nt_fasta $gene_genome_chr_out_alignment_AA_fasta $out_alignment_description --slave\n\n";

  		system(" Rscript $results_folder/scripts/curatingAlignments.r  $output_ref_chr $gene_tag $gene_genome_chr_out_alignment_nt_fasta $gene_genome_chr_out_alignment_AA_fasta $out_alignment_description --slave");

		}else { 
 my $out_alignment_description=$dirresults."/".$genome.".".$gene_tag.".chr.out.alignments.description";
	open OUT, ">$out_alignment_description" or die "Cannot open $out_alignment_description for writing\n";
	

	my $firstline_Ref=`grep ">" $gene_id.fasta`;
	chomp($firstline_Ref);
		$firstline_Ref=~s/[\n\r]//g;
		$firstline_Ref=~s/>//g;
	print OUT "$genename\t\"$firstline_Ref\"\tNA\tNA\tNA\t\"$genome\"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tchromosome\tabsent\t0\t0\t0\n";
	close OUT;
	 }
	}

}




##############################  No Plasmid #######################################


sub sub1_chr{

my $seq = shift;
my $genome=basename($seq);
	
my $dirprochr=$diroutput."/".$genome.".chromosome";
my $cmd_prod_chr="mkdir $dirprochr";
system("$cmd_prod_chr");

 }
 
sub sub2_chr{
my $seq = shift;
$seq=$diroriginalgenomes."/".$seq;
my $genome=basename($seq);
my $dirprochr=$diroutput."/".$genome.".chromosome";

my $genes_gff1=$dirprochr."/".$genome.".chr.genes.gff";
my $genes_gff1tmp=$dirprochr."/".$genome.".chr.genes.gfftmp";
my $genes_faa1=$dirprochr."/".$genome.".chr.genes.faa";
my $genes_fasta1=$dirprochr."/".$genome.".chr.genes.fasta";
my $temp1=$genes_fasta1.".tmp";

#print "\n\n\nprodigal -i $output_chr -o $gene_gbk1 -a $genes_faa1 -d $genes_fasta1\n\n\n";
#system("prodigal -i $output_chr -o $gene_gbk1 -a $genes_faa1 -d $genes_fasta1 -p meta -q");
system("prodigal -i $seq -o $genes_gff1tmp -a $genes_faa1 -d $genes_fasta1 -f gff -q  -p meta");
system("\(cat $genes_gff1tmp; echo '##FASTA'; cat  $seq\) >$genes_gff1");              
system("rm $genes_gff1tmp");    

}


sub sub3_chr{
my $seq = shift;
my $genome=basename($seq);
$seq=$diroriginalgenomes."/".$genome;
my $dirprochr=$diroutput."/".$genome.".chromosome";

my $genes_gff1=$dirprochr."/".$genome.".chr.genes.gff";
my $genes_gff1tmp=$dirprochr."/".$genome.".chr.genes.gfftmp";
my $genes_faa1=$dirprochr."/".$genome.".chr.genes.faa";
my $genes_faa1tmp=$dirprochr."/".$genome.".chr.genes.faatmp";
my $genes_fasta1=$dirprochr."/".$genome.".chr.genes.fasta";
my $genes_fasta1tmp=$dirprochr."/".$genome.".chr.genes.fastatmp";

	my $cmd_awk1="awk '/^>/{\$0=\$0\"\\:$genome\"}1'  $genes_gff1 > $genes_gff1tmp";
	my $cmd_awk2="awk '/^>/{\$0=\$0\"\\:$genome\"}1'  $genes_faa1 > $genes_faa1tmp";
	my $cmd_awk3="awk '/^>/{\$0=\$0\"\\:$genome\"}1'  $genes_fasta1 > $genes_fasta1tmp";
	
	#print ">>>>> 1 $cmd_awk\n\n";
	system($cmd_awk1);
	system($cmd_awk2);
	system($cmd_awk3);
	
	system("mv $genes_gff1tmp $genes_gff1");
	system("mv $genes_faa1tmp $genes_faa1");
	system("mv $genes_fasta1tmp $genes_fasta1");
	
}





sub sub4_chr{
my $seq = shift;
$seq=$diroriginalgenomes."/".$seq;
my $genome=basename($seq);

my $dirprochr=$diroutput."/".$genome.".chromosome";
my $genes_gff1=$dirprochr."/".$genome.".chr.genes.gff";
my $genes_faa1=$dirprochr."/".$genome.".chr.genes.faa";
my $genes_fasta1=$dirprochr."/".$genome.".chr.genes.fasta";
my $temp1=$genes_fasta1.".tmp";
my $temp2=$genes_faa1.".tmp";

if(-s $genes_fasta1){
system("cp $genes_fasta1 $temp1");
my $cmd1="awk '/^>/ {printf(\"\\n%s\\n\",\$0);next; } { printf(\"\%s\",\$0);}  END {printf(\"\\n\");}' < $temp1  > $genes_fasta1";
#print $cmd1."\n";
system( $cmd1);

system("cp $genes_fasta1 $temp2");
my $cmd2="awk '/^>/ {printf(\"\\n%s\\n\",\$0);next; } { printf(\"\%s\",\$0);}  END {printf(\"\\n\");}' < $temp2  > $genes_faa1";
#print $cmd1."\n";
system( $cmd2);


system("sed 's/ # /###/g' $genes_fasta1 >$temp1");
system("mv $temp1 $genes_fasta1");

system("sed 's/ # /###/g' $genes_faa1 >$temp2");
system("mv $temp2 $genes_faa1");
}
}



sub sub5_chr{
my $seq = shift;
$seq=$diroriginalgenomes."/".$seq;
my $genome=basename($seq);

my $dirprochr=$diroutput."/".$genome.".chromosome";
my $genes_gff1=$dirprochr."/".$genome.".chr.genes.gff";
my $genes_faa1=$dirprochr."/".$genome.".chr.genes.faa";
my $genes_fasta1=$dirprochr."/".$genome.".chr.genes.fasta";
my $temp1=$genes_fasta1.".tmp";



	foreach my $gene(@genes){
	
	my $genename=basename($gene);
	my $gene_id=$diroriginalgenes."/".$genename;
	$gene_id=~s/.fasta$//g;
	
	my $gene_tag=basename($gene);
	$gene_tag=~s/.fasta$//g;
	
	
my $outputblast1= $dirprochr."/".$genome.".".$gene_tag.".Blast.txt";
my $tempblast=$outputblast1.".tmp";

	 
  if ($genes_type eq "nucl"){
		if(-s $genes_fasta1){
		my $cmd_blast1="blastn -query $genes_fasta1 -db $gene_id -dust no -outfmt \"6 qseqid sseqid pident qlen slen length qstart qend sstart send mismatch gapope evalue bitscore\" -perc_identity $percentage -out $outputblast1";
		#print "$cmd_blast1\n";
		system ($cmd_blast1);
		
		}
		#print "\n\n\n\n find $diroutput/*/ -size 0 |  xargs rm \n\n\n\n ";
		#system("find $diroutput/*/ -size 0 |  xargs rm");
}elsif ($genes_type eq "prot"){
	if(-s $genes_fasta1){
		my $cmd_blast1="blastp -query $genes_faa1 -db $gene_id -outfmt \"6 qseqid sseqid pident qlen slen length qstart qend sstart send evalue bitscore qcovs qcovhsp scovhsp\" -out $outputblast1";
		print "$cmd_blast1\n";
		system ($cmd_blast1);
		my $cmd_highper="awk '\$3>=".$percentage." && \$13>=".$percentage."  && \$13>=".$percentage." {print \$0}' < $outputblast1 >$tempblast";
		system ($cmd_highper);
		system (" mv $tempblast $outputblast1")
		}
		#print "\n\n\n\n find $diroutput/*/ -size 0 |  xargs rm \n\n\n\n ";
		#system("find $diroutput/*/ -size 0 |  xargs rm");
		}
		
	}
}


sub sub5_1chr{
my $seq = shift;
$seq=$diroriginalgenomes."/".$seq;
my $genome=basename($seq);

my $dirprochr=$diroutput."/".$genome.".chromosome";
my $genes_gff1=$dirprochr."/".$genome.".chr.genes.gff";
my $genes_faa1=$dirprochr."/".$genome.".chr.genes.faa";
my $genes_fasta1=$dirprochr."/".$genome.".chr.genes.fasta";
my $temp1=$genes_fasta1.".tmp";



	foreach my $gene(@genes){
	
	my $genename=basename($gene);
	my $gene_id=$diroriginalgenes."/".$genename;
	$gene_id=~s/.fasta$//g;
	
	my $gene_tag=basename($gene);
	$gene_tag=~s/.fasta$//g;
	
	
my $outputblast1= $dirprochr."/".$genome.".".$gene_tag.".Blast.txt";
my $tempblast=$outputblast1.".tmp";

	 
  if ($genes_type eq "nucl"){
		if(-s $genes_fasta1){

		my $cmd_highper="awk '\$3>=".$percentage." && \$13>=".$percentage." {print \$0}' < $outputblast1 >$tempblast";
		system ($cmd_highper);
		system (" mv $tempblast $outputblast1")
		
		}
		#print "\n\n\n\n find $diroutput/*/ -size 0 |  xargs rm \n\n\n\n ";
		#system("find $diroutput/*/ -size 0 |  xargs rm");
}elsif ($genes_type eq "prot"){
	if(-s $genes_fasta1){

		my $cmd_highper="awk '\$3>=".$percentage." && \$13>=".$percentage." {print \$0}' < $outputblast1 >$tempblast";
		system ($cmd_highper);
		system (" mv $tempblast $outputblast1")
		}
		#print "\n\n\n\n find $diroutput/*/ -size 0 |  xargs rm \n\n\n\n ";
		#system("find $diroutput/*/ -size 0 |  xargs rm");
		}
		
	}

}




sub sub6_chr{
my $seq = shift;
$seq=$diroriginalgenomes."/".$seq;
my $genome=basename($seq);
my $dirprochr=$diroutput."/".$genome.".chromosome";
my $genes_gff1=$dirprochr."/".$genome.".chr.genes.gff";
my $genes_faa1=$dirprochr."/".$genome.".chr.genes.faa";
my $genes_fasta1=$dirprochr."/".$genome.".chr.genes.fasta";
my $temp1=$genes_fasta1.".tmp";

	foreach my $gene(@genes){
	
	my $genename=basename($gene);
	my $gene_id=$diroriginalgenes."/".$genename;
	$gene_id=~s/.fasta$//g;
	
	my $gene_tag=basename($gene);
	$gene_tag=~s/.fasta$//g;

	my $outputblast1= $dirprochr."/".$genome.".".$gene_tag.".Blast.txt";
	my $gene_ref=$gene_id.".fasta";

print($outputblast1);

if(-s $outputblast1){
my $temp_chr=$dirprochr.".".$gene_tag.".tmp";
	system("cut -f 1 $outputblast1 >$temp_chr");

	my $output_gene1= $dirprochr."/".$gene_tag.".fasta";
	
	if( $genes_type eq 'nucl'){
	system("xargs samtools faidx $genes_fasta1< $temp_chr >$output_gene1");
}elsif($genes_type eq 'prot'){
system("xargs samtools faidx $genes_faa1< $temp_chr >$output_gene1");}
	
	
	
	system("rm $temp_chr");
	my $output_ref_chr=$dirprochr."/".$genome.".".$gene_tag.".plusRef.fasta";
	system("cat $gene_ref $output_gene1 >$output_ref_chr");

my $out_alignment_description=$dirresults."/".$genome.".".$gene_tag.".chr.out.alignments.description";
	open OUT, ">$out_alignment_description" or die "Cannot open $out_alignment_description for writing\n";

my $gene_genome_chr_out_alignment_AA_fasta=$dirresultspep."/".$gene_tag."_".$genome.".AAsequence.fasta";
		open OUT1, ">$gene_genome_chr_out_alignment_AA_fasta" or die "Cannot open $gene_genome_chr_out_alignment_AA_fasta for writing\n";

my $gene_genome_chr_out_alignment_nt_fasta=$dirresults."/".$gene_tag."_".$genome.".nt_alignment.fasta.tmp";
		open OUT2, ">$gene_genome_chr_out_alignment_nt_fasta" or die "Cannot open  $gene_genome_chr_out_alignment_nt_fasta for writing\n";

	  	
	  	#print " Rscript curatingAlignments.r  $output_ref_chr $gene_tag $gene_genome_chr_out_alignment_nt_fasta $gene_genome_chr_out_alignment_AA_fasta $out_alignment_description --slave\n\n";

  		system(" Rscript $results_folder/scripts/curatingAlignments.r  $output_ref_chr $gene_tag $gene_genome_chr_out_alignment_nt_fasta $gene_genome_chr_out_alignment_AA_fasta $out_alignment_description --slave");

		}else { 
 my $out_alignment_description=$dirresults."/".$genome.".".$gene_tag.".chr.out.alignments.description";
	open OUT, ">$out_alignment_description" or die "Cannot open $out_alignment_description for writing\n";
	

	my $firstline_Ref=`grep ">" $gene_id.fasta`;
	chomp($firstline_Ref);
		$firstline_Ref=~s/[\n\r]//g;
		$firstline_Ref=~s/>//g;
	print OUT "$genename\t\"$firstline_Ref\"\tNA\tNA\tNA\t\"$genome\"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tchromosome\tabsent\t0\t0\t0\n";
	close OUT;


 }
	}
}




#########################################################################################

#afork, the first argument is an array - a child will be 
#forked for each array element. The second argument indicates the maximum 
#number of children that may be alive at one time. The third argument is a 
#code reference; this is the code that will be executed by the child. One 
#argument will be given to this code fragment; for mfork it will be an increasing number,
#starting at one. Each next child gets the next number. For afork, the array element is 
#passed. Note that this code will assume no other children will be spawned, 
#and that $SIG {CHLD} hasn't been set to IGNORE. 

sub afork (\@$&) {
my ($data, $max, $code) = @_;
my $c = 0;
foreach my $data (@$data) {
wait unless ++ $c <= $max;
die "Fork failed: $!\n" unless defined (my $pid = fork);
exit $code -> ($data) unless $pid;

}

my $child_pid = wait();
my $status = $?;

1 until -1 == wait;
if    ( $status == -1  ) { die "wait failed: $!\n"; }
elsif ( $status & 127) { warn "Child $child_pid killed by signal ".(  $? & 127)."\n"; }
elsif ( $status & 128) { warn "Child $child_pid. There was a core dump!\n ".(  $? & 128)."\n"; }
elsif ( $status & 0x7F ) { warn "Child $child_pid killed by signal ".( $? & 0x7F )."\n"; }
elsif ( $status >> 8   ) { warn "Child $child_pid exited with error $data".( $? >> 8 )."\n"; }
else                     { print "Child $child_pid exited successfully\n"; }

}
#######
##usage afork(@array,$numberofcores,\&subroutine);

#sub normal(\@&){

#my ($data, $code) = @_;
# foreach my $sequence(@$data){
#  $code -> ($sequence);
# }
#}
#normal(@seq,\&sub);




###################################CREATING R SCRIPTS
 if($species ne "NO"){
 $R_file="$results_folder/scripts/plasmid_chromosome_Detection.r";
	
		open R_SCRIPT,">$R_file" or die "Cannot write $R_file script\n";

        	 $R_script= "
    rm(list=ls()); 
    suppressMessages(suppressWarnings(library(mlplasmids)))
    options(warn=-1)
    options(stringsAsFactors = FALSE)
    args = commandArgs(trailingOnly=TRUE)

    filepath=args[1]
    dirproplas=args[3]
    dirprochr=args[2]

	all=read.table(filepath)
	colnames(all)=c('Contig_name','Probability')

	all\$Prediction='Chromosome'
	all\$Prediction[all\$Probability>0.75]='Plasmid'

x=all[ all['Prediction']=='Plasmid',]
	if( dim(x)[1]>0){
		namefull=basename(filepath)
		name=gsub('_plas.prob.out','',namefull)
		x['name']=name	
		nameplasmid=paste(name,'plasmidsummary.txt',sep='.')
		namelist=paste(name,'plasmid_contigslist.txt',sep='.')


		namefile=paste(dirproplas,nameplasmid,sep='/')
		namefilelist=paste(dirproplas,namelist,sep='/')
		#print (namefile)
		#print (namefilelist)
 
    	write.table(x,namefile,quote=F,row.names=F,sep='\\t')
		write.table(x['Contig_name'],namefilelist,quote=F,row.names=F,sep='\\t',col.names = F)
 	}


x=all[ all['Prediction']=='Chromosome',]
	if( dim(x)[1]>0){
		namefull=basename(filepath)
		name=gsub('_plas.prob.out','',namefull)
		x['name']=name
		nameplasmid=paste(name,'chromosomesummary.txt',sep='.')
		namelist=paste(name,'chromosome_contigslist.txt',sep='.')
	

		namefile=paste(dirprochr,nameplasmid,sep='/')
		namefilelist=paste(dirprochr,namelist,sep='/')
		#print (namefile)
		#print (namefilelist)
 
		write.table(x,namefile,quote=F,row.names=F,sep='\\t')
		write.table(x['Contig_name'],namefilelist,quote=F,row.names=F,sep='\\t',col.names = F)
}
";

       print R_SCRIPT $R_script;
        close R_SCRIPT;
        system("chmod a+x $R_file"); }else{ print "No predicting Plasmids\n\n"}
  		
	
	
	
	
if ( $genes_type eq 'nucl'){

$R_file="$results_folder/scripts/curatingAlignments.r";

		open R_SCRIPT,">$R_file" or die "Cannot write  $R_file script\n";

        	 $R_script= "
        	
        	 rm(list=ls()); 
        	 
        	options(stringsAsFactors = FALSE)
        	args = commandArgs(trailingOnly=TRUE)
			file=args[1]
			genename=args[2]
			out1=args[3]
			out2=args[4]
			out3=args[5]
			setwd(dirname(file))
			#print(getwd())
				Sys.setenv('R_MAX_VSIZE'=32000000000)

				suppressMessages(suppressWarnings(library(msa)))
				suppressMessages(suppressWarnings(library(reshape2)))
				suppressMessages(suppressWarnings(library(seqinr)))
				 suppressMessages(suppressWarnings(library( Biostrings)))
				options(warn=-1)
				    library(stringr)
					##Functions
	
	

findPotentialStartsAndStops <- function(DNA_string)
  { positions =c()
   	types =c()
     if(length(positions)==0){
     codons            <- c('atg', 'taa', 'tag', 'tga','ATG', 'TAA', 'TAG', 'TGA')
     for (i in  1:length(codons))
     {
        codon <- codons[i]
        occurrences <- matchPattern(codon, DNA_string)
      codonpositions <- occurrences\@ranges\@start  
        numoccurrences <- length(codonpositions) 	
   	positions   <- append(positions,codonpositions, after=length(positions))
     types       <- append(types,rep(codon, numoccurrences), after=length(types))
 }
    indices <- order(positions)
     positions <- positions[indices]
     types <- types[indices]
     mylist <- list(positions,types)
 return(mylist)
 
     }
   }     

 
  findORFsinSeq <- function(sequence)
  {
     mylist <- findPotentialStartsAndStops(sequence)
     positions <- mylist[[1]]
     types <- mylist[[2]]
     orfstarts <- numeric()
     orfstops <- numeric()
     orflengths <- numeric()
     numpositions <- length(positions)
 
     if (numpositions >= 2)
     {
        for (i in 1:(numpositions-1))
        {
           posi <- positions[i]
           typei <- types[i]
           found <- 0
           while (found == 0)
           {
              for (j in (i+1):numpositions)
              {
                 posj  <- positions[j]
                 typej <- types[j]
                 posdiff <- posj - posi
                 posdiffmod3 <- posdiff %% 3
           
                 orflength <- posj - posi + 3
                 if ((typei == 'atg' || typei == 'ATG') && (typej == 'taa' || typej == 'tag' || typej == 'tga'||typej == 'TAA' || typej == 'TAG' || typej == 'TGA') && posdiffmod3 == 0)
                 {
                    numorfs <- length(orfstops)
                    usedstop <- -1
                    if (numorfs > 0)
                    {
                      for (k in 1:numorfs)
                      {
                          orfstopk <- orfstops[k]
                          if (orfstopk == (posj + 2)) { usedstop <- 1 }
                      }
                    }
                    if (usedstop == -1)
                    {
                       orfstarts <- append(orfstarts, posi, after=length(orfstarts))
                       orfstops <- append(orfstops, posj+2, after=length(orfstops)) 
                       orflengths <- append(orflengths, orflength, after=length(orflengths))
                    }
                    found <- 1
                    break
                 }
                 if (j == numpositions) { found <- 1 }
              }
           }
        }
     }
     indices <- order(orfstarts)
     orfstarts <- orfstarts[indices]
     orfstops <- orfstops[indices]
     orflengths <- numeric()
     numorfs <- length(orfstarts)
     for (i in 1:numorfs)
     {
        orfstart <- orfstarts[i]
        orfstop <- orfstops[i]
        orflength <- orfstop - orfstart + 1
        orflengths <- append(orflengths,orflength,after=length(orflengths))
     }
     mylist <- data.frame(orfstarts=orfstarts,orfstops= orfstops,orflengths= orflengths)
     return(mylist)
  }
  
			seq2Fasta <- function(seq,name,filename) 
				{
  						sink(filename,append=TRUE)
    					cat(paste0('>', name),file=filename,append=TRUE)
    					cat('\\n',file=filename,append=TRUE)
    					the.sequence <- toString(seq)
    					the.sequence=gsub('N','-',the.sequence)
    					cat(the.sequence,file=filename,append=TRUE)
    					cat('\\n',file=filename,append=TRUE)  
  						sink(NULL)
				}
				
						seqAA2Fasta <- function(seq,name,filename) 
				{
  						sink(filename,append=TRUE)
    					cat(paste0('>', name),file=filename,append=TRUE)
    					cat('\\n',file=filename,append=TRUE)
    					the.sequence <- toString(seq)
    					cat(the.sequence,file=filename,append=TRUE)
    					cat('\\n',file=filename,append=TRUE)  
  						sink(NULL)
				}

		Alignment <- function(x){
    	tryCatch(
        expr = {
 		return(invisible(msa(x,gapExtension =150,gap=1500,method='Muscle')))
            message('Successfully executed')
       		 },
        error = function(e){
            message('Caught an error!')
		return(invisible(msa(x)))

        },
        warning = function(w){
           message('Caught an warning!')
            #print(w)
        },
        finally = {
            message('All done, quitting.')
      		  }
    		)    
		}


		mySequences <- readDNAStringSet(file)
		namesunique=unique(mySequences\@ranges\@NAMES)
		mySequences <- mySequences[namesunique]
		
		myFirstAlignment=Alignment(mySequences)

		myFirstAlignment_reverse=Alignment(c(reverseComplement(mySequences[1]),mySequences[2]))

 		FA=msaConsensusSequence(myFirstAlignment)
 		FAR=msaConsensusSequence(myFirstAlignment_reverse)

		Consensus=str_count(FA, pattern = '\\\\?')
 		Consensus_reverse=str_count(FAR, pattern ='\\\\?')
	
		if(Consensus<Consensus_reverse){myFirstAlignment=myFirstAlignment}else{myFirstAlignment=myFirstAlignment_reverse}		

		#check for several genes
		location=ifelse(grepl('plasmid',file),'plasmid','chromosome')


		nSeq=length(rownames(myFirstAlignment))
		#print(paste('number found',nSeq,sep=' '))
	
  if(nSeq==2){
  		#print('2 sequences')
      		
     X_Alignment=myFirstAlignment

query=rownames(X_Alignment)[2]
      		Seq=toString(unmasked(X_Alignment)[[query]] )
      		Seq_string=strsplit(toString(Seq),'')
  			x=table(strsplit(Seq,''))
  			gaps=x['-']
  			if(is.na(gaps)){gaps=0}
  			Seq_nogaps=gsub('-','',Seq)
  			x=table(strsplit(Seq_nogaps,''))
  			xORF=findORFsinSeq(Seq_nogaps)
  			bestORF=xORF[ xORF\$orflengths==max(xORF\$orflengths),]
			query_trans=Biostrings::translate(DNAString(substring(Seq_nogaps, bestORF\$orfstart,bestORF\$orfstops)),if.fuzzy.codon='solve')
      		Seq_trl=query_trans
  			length_Seqtrl=length(Seq_trl)
      		Seq_trl_string=strsplit(toString(Seq_trl),'')	
      		stops=which( Seq_trl_string[[1]]=='*')
  		

			######
			ref=rownames(X_Alignment)[1]
      		Seqref=toString(unmasked(X_Alignment)[[ref]] )

      		Seqref_string=strsplit(toString(Seqref),'')
      		xref=table(Seqref_string)
      		gapsref=xref['-']
      		if(is.na(gapsref)){gapsref=0}
      		Seq_nogapsref=gsub('-','',Seqref)
      		xref=table(strsplit(Seq_nogapsref,''))
      		xORFref=findORFsinSeq(Seq_nogapsref)
			bestORF=xORFref[ xORFref\$orflengths==max(xORFref\$orflengths),]
			ref_trans=Biostrings::translate(DNAString(substring(Seq_nogapsref, bestORF\$orfstart,bestORF\$orfstops)),if.fuzzy.codon='solve')
      		Seq_trlref=ref_trans
      		length_Seqtrlref=length(ref_trans)
      		Seq_trl_stringref=strsplit(toString(Seq_trlref),'')	
      		stopsref=which( Seq_trl_stringref[[1]]=='*')

######	
			Align=data.frame(Ref=c(unlist(Seqref_string)),Query=c(unlist(Seq_string)))
				Align\$position=c(1:dim(Align)[1])

	
			Align\$Ref=as.character(Align\$Ref)
			Align\$Query=as.character(Align\$Query)
			Align\$mutation=paste(Align\$Ref,Align\$position,Align\$Query,sep='')
			Align\$Change='identical'
			Align\$Change[ Align\$Query!=Align\$Ref ]='SNP'
			Align\$Change[ Align\$Query=='-' | Align\$Query=='N' ]='GAP'
			Align\$Change[ Align\$Ref=='-' | Align\$Ref=='N' ]='INSERTION'
			
			
			changesSNPs=Align[Align\$Change=='SNP',]
			changesGAPs=Align[Align\$Change=='GAP',]
			changesINSERTION=Align[Align\$Change=='INSERTION',]
			SNPs=paste(changesSNPs\$mutation,collapse=',')
			GAPs=paste(changesGAPs\$mutation,collapse=',')
			INSERTION=paste(changesINSERTION\$mutation,collapse=',')
			
			myAlignment=myFirstAlignment
			Ncopies=1
    }
    
      k=0
 	 if(nSeq==3){
 		#print('3 sequences')
 	 sequences=length(mySequences\@ranges\@NAMES)
 	  xx <- as.data.frame(matrix(ncol=sequences, nrow=myFirstAlignment\@unmasked\@ranges\@width[1]))
 
	  names=matrix(ncol=sequences, nrow=1)
	 
    		for( i in mySequences\@ranges\@NAMES){ 
        	x=c(unlist(strsplit(toString(unmasked(myFirstAlignment)[[i]]),'')))
        	k=k+1
    		xx[,k]=x
        	names[,k]=i
    		}
  
    		
    		xx\$consensus=apply(xx, 1, function(x) {
						y=unique(x[-1])	
    					y=y[y!='-']	
						if(length(y)==0){y='-'}
						if(length(y)>1){if(x[1]==x[2]){y=x[1]}}
						if(length(y)>1){if(x[1]==x[3]){y=x[1]}}
						if(length(y)>1){y='N'}
						return(y)})

		conSeq=paste(xx\$consensus,collapse='')
		conSeq=gsub('-','',conSeq)
  		consensusSeq=DNAStringSet(conSeq)
  		
  		consensusSeq\@ranges\@NAMES=mySequences\@ranges\@NAMES[2]
  		ref=DNAStringSet(mySequences[[1]])
  		ref\@ranges\@NAMES=mySequences\@ranges\@NAMES[1]
  		
  		myAlignment=invisible(Alignment(c(ref,consensusSeq)))

		myAlignment_reverse=Alignment(c(reverseComplement(ref),consensusSeq))

 		FA=msaConsensusSequence(myAlignment)
 		FAR=msaConsensusSequence(myAlignment_reverse)

		Consensus=str_count(FA, pattern = '\\\\?')
 		Consensus_reverse=str_count(FAR, pattern ='\\\\?')
	
		if(Consensus<Consensus_reverse){Alignment=myAlignment}else{myAlignment=myAlignment_reverse}
	
     X_Alignment=myAlignment

query=rownames(X_Alignment)[2]
      		Seq=toString(unmasked(X_Alignment)[[query]] )
      		Seq_string=strsplit(toString(Seq),'')
  			x=table(strsplit(Seq,''))
  			gaps=x['-']
  			if(is.na(gaps)){gaps=0}
  			Seq_nogaps=gsub('-','',Seq)
  			x=table(strsplit(Seq_nogaps,''))
  			xORF=findORFsinSeq(Seq_nogaps)
  			bestORF=xORF[ xORF\$orflengths==max(xORF\$orflengths),]
			query_trans=Biostrings::translate(DNAString(substring(Seq_nogaps, bestORF\$orfstart,bestORF\$orfstops)),if.fuzzy.codon='solve')
      		Seq_trl=query_trans
  			length_Seqtrl=length(Seq_trl)
      		Seq_trl_string=strsplit(toString(Seq_trl),'')	
      		stops=which( Seq_trl_string[[1]]=='*')
  		

			######
			ref=rownames(X_Alignment)[1]
      		Seqref=toString(unmasked(X_Alignment)[[ref]] )

      		Seqref_string=strsplit(toString(Seqref),'')
      		xref=table(Seqref_string)
      		gapsref=xref['-']
      		if(is.na(gapsref)){gapsref=0}
      		Seq_nogapsref=gsub('-','',Seqref)
      		xref=table(strsplit(Seq_nogapsref,''))
      		xORFref=findORFsinSeq(Seq_nogapsref)
			bestORF=xORFref[ xORFref\$orflengths==max(xORFref\$orflengths),]
			ref_trans=Biostrings::translate(DNAString(substring(Seq_nogapsref, bestORF\$orfstart,bestORF\$orfstops)),if.fuzzy.codon='solve')
      		Seq_trlref=ref_trans
      		length_Seqtrlref=length(ref_trans)
      		Seq_trl_stringref=strsplit(toString(Seq_trlref),'')	
      		stopsref=which( Seq_trl_stringref[[1]]=='*')

######	
			Align=data.frame(Ref=c(unlist(Seqref_string)),Query=c(unlist(Seq_string)))
				Align\$position=c(1:dim(Align)[1])

	
			Align\$Ref=as.character(Align\$Ref)
			Align\$Query=as.character(Align\$Query)
			Align\$mutation=paste(Align\$Ref,Align\$position,Align\$Query,sep='')
			Align\$Change='identical'
			Align\$Change[ Align\$Query!=Align\$Ref ]='SNP'
			Align\$Change[ Align\$Query=='-' | Align\$Query=='N' ]='GAP'
			Align\$Change[ Align\$Ref=='-' | Align\$Ref=='N' ]='INSERTION'
			
			
			changesSNPs=Align[Align\$Change=='SNP',]
			changesGAPs=Align[Align\$Change=='GAP',]
			changesINSERTION=Align[Align\$Change=='INSERTION',]
			SNPs=paste(changesSNPs\$mutation,collapse=',')
			GAPs=paste(changesGAPs\$mutation,collapse=',')
			INSERTION=paste(changesINSERTION\$mutation,collapse=',')

			if(any(xx[,2]=='-') && any(xx[,3]=='-')){Ncopies=1}else{Ncopies=2}		
   
}




k=0
 if(nSeq==4){
 		#print('4 sequences')
 
  sequences=length(mySequences\@ranges\@NAMES)
 	  xx <- as.data.frame(matrix(ncol=sequences, nrow=myFirstAlignment\@unmasked\@ranges\@width[1]))
 
	  names=matrix(ncol=sequences, nrow=1)
	 
    		for( i in mySequences\@ranges\@NAMES){ 
        	x=c(unlist(strsplit(toString(unmasked(myFirstAlignment)[[i]]),'')))
        	k=k+1
    		xx[,k]=x
        	names[,k]=i
    		}

		xx\$consensus=apply(xx, 1, function(x) {
						y=unique(x[-1])	
    					y=y[y!='-']	
						if(length(y)==0){y='-'}
						if(length(y)>1){if(x[1]==x[2]){y=x[1]}}
						if(length(y)>1){if(x[1]==x[3]){y=x[1]}}
						if(length(y)>1){if(x[1]==x[4]){y=x[1]}}
						if(length(y)>1){y='N'}
						return(y)})

  		
  		conSeq=paste(xx\$consensus,collapse='')
		conSeq=gsub('-','',conSeq)
  		consensusSeq=DNAStringSet(conSeq)
  		consensusSeq\@ranges\@NAMES=mySequences\@ranges\@NAMES[2]
  		ref=DNAStringSet(mySequences[[1]])
  		ref\@ranges\@NAMES=mySequences\@ranges\@NAMES[1]
  		
  		myAlignment=invisible(Alignment(c(ref,consensusSeq)))
		myAlignment_reverse=Alignment(c(reverseComplement(ref),consensusSeq))

 		FA=msaConsensusSequence(myAlignment)
 		FAR=msaConsensusSequence(myAlignment_reverse)

		Consensus=str_count(FA, pattern = '\\\\?')
 		Consensus_reverse=str_count(FAR, pattern ='\\\\?')
	
		if(Consensus<Consensus_reverse){Alignment=myAlignment}else{myAlignment=myAlignment_reverse}
  

     X_Alignment=myAlignment

query=rownames(X_Alignment)[2]
      		Seq=toString(unmasked(X_Alignment)[[query]] )
      		Seq_string=strsplit(toString(Seq),'')
  			x=table(strsplit(Seq,''))
  			gaps=x['-']
  			if(is.na(gaps)){gaps=0}
  			Seq_nogaps=gsub('-','',Seq)
  			x=table(strsplit(Seq_nogaps,''))
  			xORF=findORFsinSeq(Seq_nogaps)
  			bestORF=xORF[ xORF\$orflengths==max(xORF\$orflengths),]
			query_trans=Biostrings::translate(DNAString(substring(Seq_nogaps, bestORF\$orfstart,bestORF\$orfstops)),if.fuzzy.codon='solve')
      		Seq_trl=query_trans
  			length_Seqtrl=length(Seq_trl)
      		Seq_trl_string=strsplit(toString(Seq_trl),'')	
      		stops=which( Seq_trl_string[[1]]=='*')
  		

			######
			ref=rownames(X_Alignment)[1]
      		Seqref=toString(unmasked(X_Alignment)[[ref]] )

      		Seqref_string=strsplit(toString(Seqref),'')
      		xref=table(Seqref_string)
      		gapsref=xref['-']
      		if(is.na(gapsref)){gapsref=0}
      		Seq_nogapsref=gsub('-','',Seqref)
      		xref=table(strsplit(Seq_nogapsref,''))
      		xORFref=findORFsinSeq(Seq_nogapsref)
			bestORF=xORFref[ xORFref\$orflengths==max(xORFref\$orflengths),]
			ref_trans=Biostrings::translate(DNAString(substring(Seq_nogapsref, bestORF\$orfstart,bestORF\$orfstops)),if.fuzzy.codon='solve')
      		Seq_trlref=ref_trans
      		length_Seqtrlref=length(ref_trans)
      		Seq_trl_stringref=strsplit(toString(Seq_trlref),'')	
      		stopsref=which( Seq_trl_stringref[[1]]=='*')

######	
			Align=data.frame(Ref=c(unlist(Seqref_string)),Query=c(unlist(Seq_string)))
				Align\$position=c(1:dim(Align)[1])

	
			Align\$Ref=as.character(Align\$Ref)
			Align\$Query=as.character(Align\$Query)
			Align\$mutation=paste(Align\$Ref,Align\$position,Align\$Query,sep='')
			Align\$Change='identical'
			Align\$Change[ Align\$Query!=Align\$Ref ]='SNP'
			Align\$Change[ Align\$Query=='-' | Align\$Query=='N' ]='GAP'
			Align\$Change[ Align\$Ref=='-' | Align\$Ref=='N' ]='INSERTION'
			
			
			changesSNPs=Align[Align\$Change=='SNP',]
			changesGAPs=Align[Align\$Change=='GAP',]
			changesINSERTION=Align[Align\$Change=='INSERTION',]
			SNPs=paste(changesSNPs\$mutation,collapse=',')
			GAPs=paste(changesGAPs\$mutation,collapse=',')
			INSERTION=paste(changesINSERTION\$mutation,collapse=',')

			if(any(xx[,2]=='-') && any(xx[,3]=='-') && any(xx[,4]=='-')){Ncopies=1}else{Ncopies=3}		
			if(any(xx[,2]=='-') && any(xx[,3]=='-') ){Ncopies=1}else{Ncopies=2}		
			
}


k=0
 if(nSeq==5){
 		#print('5 sequences')
 
  sequences=length(mySequences\@ranges\@NAMES)
 	  xx <- as.data.frame(matrix(ncol=sequences, nrow=myFirstAlignment\@unmasked\@ranges\@width[1]))
 
	  names=matrix(ncol=sequences, nrow=1)
	 
    		for( i in mySequences\@ranges\@NAMES){ 
        	x=c(unlist(strsplit(toString(unmasked(myFirstAlignment)[[i]]),'')))
        	k=k+1
    		xx[,k]=x
        	names[,k]=i
    		}
    				xx\$consensus=apply(xx, 1, function(x) {
						y=unique(x[-1])	
    					y=y[y!='-']	
						if(length(y)==0){y='-'}
						if(length(y)>1){if(x[1]==x[2]){y=x[1]}}
						if(length(y)>1){if(x[1]==x[3]){y=x[1]}}
						if(length(y)>1){if(x[1]==x[4]){y=x[1]}}
						if(length(y)>1){if(x[1]==x[5]){y=x[1]}}
						if(length(y)>1){y='N'}
						return(y)})

  		
  		conSeq=paste(xx\$consensus,collapse='')
		conSeq=gsub('-','',conSeq)
  		consensusSeq=DNAStringSet(conSeq)
  		consensusSeq\@ranges\@NAMES=mySequences\@ranges\@NAMES[2]
  		
  		ref=DNAStringSet(mySequences[[1]])
  		ref\@ranges\@NAMES=mySequences\@ranges\@NAMES[1]
  		
  		myAlignment=invisible(Alignment(c(ref,consensusSeq)))
		myAlignment_reverse=Alignment(c(reverseComplement(ref),consensusSeq))

 		FA=msaConsensusSequence(myAlignment)
 		FAR=msaConsensusSequence(myAlignment_reverse)

		Consensus=str_count(FA, pattern = '\\\\?')
 		Consensus_reverse=str_count(FAR, pattern ='\\\\?')
	
		if(Consensus<Consensus_reverse){Alignment=myAlignment}else{myAlignment=myAlignment_reverse}
       
       X_Alignment=myAlignment

			query=rownames(X_Alignment)[2]
      		Seq=toString(unmasked(X_Alignment)[[query]] )
      		Seq_string=strsplit(toString(Seq),'')
  			x=table(strsplit(Seq,''))
  			gaps=x['-']
  			if(is.na(gaps)){gaps=0}
  			Seq_nogaps=gsub('-','',Seq)
  			x=table(strsplit(Seq_nogaps,''))
  			xORF=findORFsinSeq(Seq_nogaps)
  			bestORF=xORF[ xORF\$orflengths==max(xORF\$orflengths),]
			query_trans=Biostrings::translate(DNAString(substring(Seq_nogaps, bestORF\$orfstart,bestORF\$orfstops)),if.fuzzy.codon='solve')
      		Seq_trl=query_trans
  			length_Seqtrl=length(Seq_trl)
      		Seq_trl_string=strsplit(toString(Seq_trl),'')	
      		stops=which( Seq_trl_string[[1]]=='*')
  		

			######
			ref=rownames(X_Alignment)[1]
      		Seqref=toString(unmasked(X_Alignment)[[ref]] )

      		Seqref_string=strsplit(toString(Seqref),'')
      		xref=table(Seqref_string)
      		gapsref=xref['-']
      		if(is.na(gapsref)){gapsref=0}
      		Seq_nogapsref=gsub('-','',Seqref)
      		xref=table(strsplit(Seq_nogapsref,''))
      		xORFref=findORFsinSeq(Seq_nogapsref)
			bestORF=xORFref[ xORFref\$orflengths==max(xORFref\$orflengths),]
			ref_trans=Biostrings::translate(DNAString(substring(Seq_nogapsref, bestORF\$orfstart,bestORF\$orfstops)),if.fuzzy.codon='solve')
      		Seq_trlref=ref_trans
      		length_Seqtrlref=length(ref_trans)
      		Seq_trl_stringref=strsplit(toString(Seq_trlref),'')	
      		stopsref=which( Seq_trl_stringref[[1]]=='*')

######	
			Align=data.frame(Ref=c(unlist(Seqref_string)),Query=c(unlist(Seq_string)))
				Align\$position=c(1:dim(Align)[1])

	
			Align\$Ref=as.character(Align\$Ref)
			Align\$Query=as.character(Align\$Query)
			Align\$mutation=paste(Align\$Ref,Align\$position,Align\$Query,sep='')
			Align\$Change='identical'
			Align\$Change[ Align\$Query!=Align\$Ref ]='SNP'
			Align\$Change[ Align\$Query=='-' | Align\$Query=='N' ]='GAP'
			Align\$Change[ Align\$Ref=='-' | Align\$Ref=='N' ]='INSERTION'
			
			
			changesSNPs=Align[Align\$Change=='SNP',]
			changesGAPs=Align[Align\$Change=='GAP',]
			changesINSERTION=Align[Align\$Change=='INSERTION',]
			SNPs=paste(changesSNPs\$mutation,collapse=',')
			GAPs=paste(changesGAPs\$mutation,collapse=',')
			INSERTION=paste(changesINSERTION\$mutation,collapse=',')
			
			if(any(xx[,2]=='-') && any(xx[,3]=='-') && any(xx[,4]=='-') && any(xx[,5]=='-')){Ncopies=1}else{Ncopies=4}	
			if(any(xx[,2]=='-') && any(xx[,3]=='-') && any(xx[,4]=='-' )){Ncopies=1}else{Ncopies=3}		
			if(any(xx[,2]=='-') && any(xx[,3]=='-') ){Ncopies=1}else{Ncopies=2}		

}

k=0
 if(nSeq==6){
 		#print('6 sequences')
 
  sequences=length(mySequences\@ranges\@NAMES)
 	  xx <- as.data.frame(matrix(ncol=sequences, nrow=myFirstAlignment\@unmasked\@ranges\@width[1]))
 
	  names=matrix(ncol=sequences, nrow=1)
	 
    		for( i in mySequences\@ranges\@NAMES){ 
        	x=c(unlist(strsplit(toString(unmasked(myFirstAlignment)[[i]]),'')))
        	k=k+1
    		xx[,k]=x
        	names[,k]=i
    		}
    				xx\$consensus=apply(xx, 1, function(x) {
						y=unique(x[-1])	
    					y=y[y!='-']	
						if(length(y)==0){y='-'}
						if(length(y)>1){if(x[1]==x[2]){y=x[1]}}
						if(length(y)>1){if(x[1]==x[3]){y=x[1]}}
						if(length(y)>1){if(x[1]==x[4]){y=x[1]}}
						if(length(y)>1){if(x[1]==x[5]){y=x[1]}}
						if(length(y)>1){if(x[1]==x[6]){y=x[1]}}
						if(length(y)>1){y='N'}
						return(y)})

  		
  		conSeq=paste(xx\$consensus,collapse='')
		conSeq=gsub('-','',conSeq)
  		consensusSeq=DNAStringSet(conSeq)
  		consensusSeq\@ranges\@NAMES=mySequences\@ranges\@NAMES[2]
  		ref=DNAStringSet(mySequences[[1]])
  		ref\@ranges\@NAMES=mySequences\@ranges\@NAMES[1]
  		
  		myAlignment=invisible(Alignment(c(ref,consensusSeq)))
		myAlignment_reverse=Alignment(c(reverseComplement(ref),consensusSeq))

 		FA=msaConsensusSequence(myAlignment)
 		FAR=msaConsensusSequence(myAlignment_reverse)

		Consensus=str_count(FA, pattern = '\\\\?')
 		Consensus_reverse=str_count(FAR, pattern ='\\\\?')
	
		if(Consensus<Consensus_reverse){Alignment=myAlignment}else{myAlignment=myAlignment_reverse}
       
       X_Alignment=myAlignment

			query=rownames(X_Alignment)[2]
      		Seq=toString(unmasked(X_Alignment)[[query]] )
      		Seq_string=strsplit(toString(Seq),'')
  			x=table(strsplit(Seq,''))
  			gaps=x['-']
  			if(is.na(gaps)){gaps=0}
  			Seq_nogaps=gsub('-','',Seq)
  			x=table(strsplit(Seq_nogaps,''))
  			xORF=findORFsinSeq(Seq_nogaps)
  			bestORF=xORF[ xORF\$orflengths==max(xORF\$orflengths),]
			query_trans=Biostrings::translate(DNAString(substring(Seq_nogaps, bestORF\$orfstart,bestORF\$orfstops)),if.fuzzy.codon='solve')
      		Seq_trl=query_trans
  			length_Seqtrl=length(Seq_trl)
      		Seq_trl_string=strsplit(toString(Seq_trl),'')	
      		stops=which( Seq_trl_string[[1]]=='*')
  		

			######
			ref=rownames(X_Alignment)[1]
      		Seqref=toString(unmasked(X_Alignment)[[ref]] )

      		Seqref_string=strsplit(toString(Seqref),'')
      		xref=table(Seqref_string)
      		gapsref=xref['-']
      		if(is.na(gapsref)){gapsref=0}
      		Seq_nogapsref=gsub('-','',Seqref)
      		xref=table(strsplit(Seq_nogapsref,''))
      		xORFref=findORFsinSeq(Seq_nogapsref)
			bestORF=xORFref[ xORFref\$orflengths==max(xORFref\$orflengths),]
			ref_trans=Biostrings::translate(DNAString(substring(Seq_nogapsref, bestORF\$orfstart,bestORF\$orfstops)),if.fuzzy.codon='solve')
      		Seq_trlref=ref_trans
      		length_Seqtrlref=length(ref_trans)
      		Seq_trl_stringref=strsplit(toString(Seq_trlref),'')	
      		stopsref=which( Seq_trl_stringref[[1]]=='*')

######	
			Align=data.frame(Ref=c(unlist(Seqref_string)),Query=c(unlist(Seq_string)))
				Align\$position=c(1:dim(Align)[1])

	
			Align\$Ref=as.character(Align\$Ref)
			Align\$Query=as.character(Align\$Query)
			Align\$mutation=paste(Align\$Ref,Align\$position,Align\$Query,sep='')
			Align\$Change='identical'
			Align\$Change[ Align\$Query!=Align\$Ref ]='SNP'
			Align\$Change[ Align\$Query=='-' | Align\$Query=='N' ]='GAP'
			Align\$Change[ Align\$Ref=='-' | Align\$Ref=='N' ]='INSERTION'
			
			
			changesSNPs=Align[Align\$Change=='SNP',]
			changesGAPs=Align[Align\$Change=='GAP',]
			changesINSERTION=Align[Align\$Change=='INSERTION',]
			SNPs=paste(changesSNPs\$mutation,collapse=',')
			GAPs=paste(changesGAPs\$mutation,collapse=',')
			INSERTION=paste(changesINSERTION\$mutation,collapse=',')
			
			if(any(xx[,2]=='-') && any(xx[,3]=='-') && any(xx[,4]=='-') && any(xx[,5]=='-') && any(xx[,6]=='-')){Ncopies=1}else{Ncopies=5}
			if(any(xx[,2]=='-') && any(xx[,3]=='-') && any(xx[,4]=='-') && any(xx[,5]=='-')){Ncopies=1}else{Ncopies=4}	
			if(any(xx[,2]=='-') && any(xx[,3]=='-') && any(xx[,4]=='-' )){Ncopies=1}else{Ncopies=3}		
			if(any(xx[,2]=='-') && any(xx[,3]=='-') ){Ncopies=1}else{Ncopies=2}		

}

k=0
 if(nSeq==7){
 		#print('7 sequences')
 
  sequences=length(mySequences\@ranges\@NAMES)
 	  xx <- as.data.frame(matrix(ncol=sequences, nrow=myFirstAlignment\@unmasked\@ranges\@width[1]))
 
	  names=matrix(ncol=sequences, nrow=1)
	 
    		for( i in mySequences\@ranges\@NAMES){ 
        	x=c(unlist(strsplit(toString(unmasked(myFirstAlignment)[[i]]),'')))
        	k=k+1
    		xx[,k]=x
        	names[,k]=i
    		}
    				xx\$consensus=apply(xx, 1, function(x) {
						y=unique(x[-1])	
    					y=y[y!='-']	
						if(length(y)==0){y='-'}
						if(length(y)>1){if(x[1]==x[2]){y=x[1]}}
						if(length(y)>1){if(x[1]==x[3]){y=x[1]}}
						if(length(y)>1){if(x[1]==x[4]){y=x[1]}}
						if(length(y)>1){if(x[1]==x[5]){y=x[1]}}
						if(length(y)>1){if(x[1]==x[6]){y=x[1]}}
						if(length(y)>1){if(x[1]==x[7]){y=x[1]}}
						if(length(y)>1){y='N'}
						return(y)})

  		
  		conSeq=paste(xx\$consensus,collapse='')
		conSeq=gsub('-','',conSeq)
  		consensusSeq=DNAStringSet(conSeq)
  		consensusSeq\@ranges\@NAMES=mySequences\@ranges\@NAMES[2]
  		ref=DNAStringSet(mySequences[[1]])
  		ref\@ranges\@NAMES=mySequences\@ranges\@NAMES[1]
  		
  		myAlignment=invisible(Alignment(c(ref,consensusSeq)))
		myAlignment_reverse=Alignment(c(reverseComplement(ref),consensusSeq))

 		FA=msaConsensusSequence(myAlignment)
 		FAR=msaConsensusSequence(myAlignment_reverse)

		Consensus=str_count(FA, pattern = '\\\\?')
 		Consensus_reverse=str_count(FAR, pattern ='\\\\?')
	
		if(Consensus<Consensus_reverse){Alignment=myAlignment}else{myAlignment=myAlignment_reverse}
       
      	 X_Alignment=myAlignment

			query=rownames(X_Alignment)[2]
      		Seq=toString(unmasked(X_Alignment)[[query]] )
      		Seq_string=strsplit(toString(Seq),'')
  			x=table(strsplit(Seq,''))
  			gaps=x['-']
  			if(is.na(gaps)){gaps=0}
  			Seq_nogaps=gsub('-','',Seq)
  			x=table(strsplit(Seq_nogaps,''))
  			xORF=findORFsinSeq(Seq_nogaps)
  			bestORF=xORF[ xORF\$orflengths==max(xORF\$orflengths),]
			query_trans=Biostrings::translate(DNAString(substring(Seq_nogaps, bestORF\$orfstart,bestORF\$orfstops)),if.fuzzy.codon='solve')
      		Seq_trl=query_trans
  			length_Seqtrl=length(Seq_trl)
      		Seq_trl_string=strsplit(toString(Seq_trl),'')	
      		stops=which( Seq_trl_string[[1]]=='*')
  		

			######
			ref=rownames(X_Alignment)[1]
      		Seqref=toString(unmasked(X_Alignment)[[ref]] )

      		Seqref_string=strsplit(toString(Seqref),'')
      		xref=table(Seqref_string)
      		gapsref=xref['-']
      		if(is.na(gapsref)){gapsref=0}
      		Seq_nogapsref=gsub('-','',Seqref)
      		xref=table(strsplit(Seq_nogapsref,''))
      		xORFref=findORFsinSeq(Seq_nogapsref)
			bestORF=xORFref[ xORFref\$orflengths==max(xORFref\$orflengths),]
			ref_trans=Biostrings::translate(DNAString(substring(Seq_nogapsref, bestORF\$orfstart,bestORF\$orfstops)),if.fuzzy.codon='solve')
      		Seq_trlref=ref_trans
      		length_Seqtrlref=length(ref_trans)
      		Seq_trl_stringref=strsplit(toString(Seq_trlref),'')	
      		stopsref=which( Seq_trl_stringref[[1]]=='*')

######	
			Align=data.frame(Ref=c(unlist(Seqref_string)),Query=c(unlist(Seq_string)))
				Align\$position=c(1:dim(Align)[1])

	
			Align\$Ref=as.character(Align\$Ref)
			Align\$Query=as.character(Align\$Query)
			Align\$mutation=paste(Align\$Ref,Align\$position,Align\$Query,sep='')
			Align\$Change='identical'
			Align\$Change[ Align\$Query!=Align\$Ref ]='SNP'
			Align\$Change[ Align\$Query=='-' | Align\$Query=='N' ]='GAP'
			Align\$Change[ Align\$Ref=='-' | Align\$Ref=='N' ]='INSERTION'
			
			
			changesSNPs=Align[Align\$Change=='SNP',]
			changesGAPs=Align[Align\$Change=='GAP',]
			changesINSERTION=Align[Align\$Change=='INSERTION',]
			SNPs=paste(changesSNPs\$mutation,collapse=',')
			GAPs=paste(changesGAPs\$mutation,collapse=',')
			INSERTION=paste(changesINSERTION\$mutation,collapse=',')
		
			if(any(xx[,2]=='-') && any(xx[,3]=='-') && any(xx[,4]=='-') && any(xx[,5]=='-') && any(xx[,6]=='-') && any(xx[,6]=='-')){Ncopies=1}else{Ncopies=6}
			if(any(xx[,2]=='-') && any(xx[,3]=='-') && any(xx[,4]=='-') && any(xx[,5]=='-') && any(xx[,6]=='-')){Ncopies=1}else{Ncopies=5}
			if(any(xx[,2]=='-') && any(xx[,3]=='-') && any(xx[,4]=='-') && any(xx[,5]=='-')){Ncopies=1}else{Ncopies=4}	
			if(any(xx[,2]=='-') && any(xx[,3]=='-') && any(xx[,4]=='-' )){Ncopies=1}else{Ncopies=3}		
			if(any(xx[,2]=='-') && any(xx[,3]=='-') ){Ncopies=1}else{Ncopies=2}		

}




names=mySequences\@ranges\@NAMES

if(Ncopies==1){

myAlignmentAA=invisible(msa(c(AAStringSet(Seq_trlref),AAStringSet(Seq_trl))))
RefAA=strsplit(toString(  toString(unmasked(myAlignmentAA)[[1]] ) ),'')
SeqAA=strsplit(toString(  toString(unmasked(myAlignmentAA)[[2]] ) ),'')
PepAlign=data.frame(Ref=unlist(RefAA),Query=unlist(SeqAA))
 
 
PepAlign\$match[PepAlign\$Ref== PepAlign\$Query]=1
PepAlign\$match[PepAlign\$Ref!= PepAlign\$Query]=0

myAlignmentDNA=invisible(msa(c(DNAStringSet(Seqref),DNAStringSet(Seq))))
RefDNA=strsplit(toString(  toString(unmasked(myAlignmentDNA)[[1]] ) ),'')
SeqDNA=strsplit(toString(  toString(unmasked(myAlignmentDNA)[[2]] ) ),'')
DNAAlign=data.frame(Ref=unlist(RefDNA),Query=unlist(SeqDNA))
 
 
PepAlign\$match[PepAlign\$Ref== PepAlign\$Query]=1
PepAlign\$match[PepAlign\$Ref!= PepAlign\$Query]=0

DNAAlign\$match[DNAAlign\$Ref== DNAAlign\$Query]=1
DNAAlign\$match[DNAAlign\$Ref!= DNAAlign\$Query]=0

PepAlignSimilarity=sum(PepAlign\$match/length(PepAlign\$match))*100	
DNAAlignSimilarity=sum(DNAAlign\$match/length(DNAAlign\$match))*100

query=strsplit(names[2],'###')[[1]][1]
start=strsplit(names[2],'###')[[1]][2]
end=strsplit(names[2],'###')[[1]][3]
strand=strsplit(names[2],'###')[[1]][4]
GenomeFile=strsplit(names[2],':')[[1]][2]


startingcodon= min(which(PepAlign\$match==1))

if (stopsref != stops[1]){classification='possibly_truncated'}else{classification='present'}
if (startingcodon != 1){classification='late_startcodon'}else{classification='possibly_truncated'}


   		metadatasummary=data.frame(
      		gene=genename,
      			Ref=paste('\"',ref,'\"',sep=''),
 					Ref_length_nt=sum(xref),
 					Ref_length_AA=length_Seqtrlref,
 					Ref_stopCodon=stopsref,
 					GenomeFile=GenomeFile,
 					Query=paste('\"',query,'\"',sep=''),
 					Start=start,
 					End=end,
 					Strand= strand,
 					Query_length_nt=sum(x),
 					Query_length_AA=length_Seqtrl,
 					nt_gaps=gaps,
 					stopCodon= stops[1],
 					SNPs=SNPs,
 					GAPs=GAPs,
 					INSERTION=INSERTION,
 					location=location,
 					classification=classification,
 					PepAlignSimilarity=PepAlignSimilarity,
 					DNAAlignSimilarity=DNAAlignSimilarity,
 						Ncopies=Ncopies)

SeqQueryGenomeName=paste(GenomeFile,query,sep='::')
      		seq2Fasta(Seq,SeqQueryGenomeName,out1)
       		seqAA2Fasta(Seq_trl,SeqQueryGenomeName, out2)
  write.table(metadatasummary,file=out3,append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t')
}

count=1

if(Ncopies>1){

names=names[-1]

for (i in names){

count=count+1

query=strsplit(i,'###')[[1]][1]
start=strsplit(i,'###')[[1]][2]
end=strsplit(i,'###')[[1]][3]
strand=strsplit(i,'###')[[1]][4]
GenomeFile=strsplit(names[1],':')[[1]][2]

SeqName=xx[count]
ff=gsub('\"','',SeqName)
ff=gsub('\\n','',ff)
ff=gsub(', ','',ff)
ff=gsub('c\\\\(','',ff)
ff=gsub('\\\\)','',ff)

SeqName=ff
x=table(strsplit(SeqName,''))
xORF=findORFsinSeq(SeqName)
bestORF=xORF[ xORF\$orflengths==max(xORF\$orflengths),]
query_trans=Biostrings::translate(DNAString(substring(SeqName, bestORF\$orfstart,bestORF\$orfstops)),if.fuzzy.codon='solve')
Seq_trl=query_trans
length_Seqtrl=length(Seq_trl)
Seq_trl_string=strsplit(toString(Seq_trl),'')	
stops=which( Seq_trl_string[[1]]=='*')


myAlignmentAA=invisible(msa(c(AAStringSet(Seq_trlref),AAStringSet(Seq_trl))))
RefAA=strsplit(toString(  toString(unmasked(myAlignmentAA)[[1]] ) ),'')
SeqAA=strsplit(toString(  toString(unmasked(myAlignmentAA)[[2]] ) ),'')
PepAlign=data.frame(Ref=unlist(RefAA),Query=unlist(SeqAA))
 
 
PepAlign\$match[PepAlign\$Ref== PepAlign\$Query]=1
PepAlign\$match[PepAlign\$Ref!= PepAlign\$Query]=0

myAlignmentDNA=invisible(msa(c(DNAStringSet(Seqref),DNAStringSet(SeqName))))
RefDNA=strsplit(toString(  toString(unmasked(myAlignmentDNA)[[1]] ) ),'')
SeqDNA=strsplit(toString(  toString(unmasked(myAlignmentDNA)[[2]] ) ),'')
DNAAlign=data.frame(Ref=unlist(RefDNA),Query=unlist(SeqDNA))
 
 
PepAlign\$match[PepAlign\$Ref== PepAlign\$Query]=1
PepAlign\$match[PepAlign\$Ref!= PepAlign\$Query]=0

DNAAlign\$match[DNAAlign\$Ref== DNAAlign\$Query]=1
DNAAlign\$match[DNAAlign\$Ref!= DNAAlign\$Query]=0

PepAlignSimilarity=sum(PepAlign\$match/length(PepAlign\$match))*100	
DNAAlignSimilarity=sum(DNAAlign\$match/length(DNAAlign\$match))*100

startingcodon= min(which(PepAlign\$match==1))

if (stopsref != stops[1]){classification='possibly_truncated'}else{classification='present'}
if (startingcodon != 1){classification='late_startcodon'}else{classification='present'}

  		metadatasummary=data.frame(
      		gene=genename,
      			Ref=paste('\"',ref,'\"',sep=''),
 					Ref_length_nt=sum(xref),
 					Ref_length_AA=length_Seqtrlref,
 					Ref_stopCodon=stopsref,
 					GenomeFile=GenomeFile,
 					Query=paste('\"',query,'\"',sep=''),
 					Start=start,
 					End=end,
 					Strand= strand,
 					Query_length_nt=sum(x),
 					Query_length_AA=length_Seqtrl,
 					nt_gaps=gaps,
 					stopCodon= stops[1],
 					SNPs=SNPs,
 					GAPs=GAPs,
 					INSERTION=INSERTION,
 					location=location,
 					classification=classification,
 					PepAlignSimilarity=PepAlignSimilarity,
 					DNAAlignSimilarity=DNAAlignSimilarity,
 						Ncopies=Ncopies)

SeqQueryGenomeName=paste(GenomeFile,query,sep='::')
      		seq2Fasta(SeqName,SeqQueryGenomeName,out1)
       		seqAA2Fasta(Seq_trl,SeqQueryGenomeName, out2)
  write.table(metadatasummary,file=out3,append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t')
}


#print (genename)
#print (out3)
}
";
 print R_SCRIPT $R_script;
        close R_SCRIPT;
        system("chmod a+x $R_file"); 

}


if ( $genes_type eq 'prot'){

$R_file="$results_folder/scripts/curatingAlignments.r";

		open R_SCRIPT,">$R_file" or die "Cannot write  $R_file script\n";

        	 $R_script= "
        	
        	 rm(list=ls()); 
        	 
        	options(stringsAsFactors = FALSE)
        	args = commandArgs(trailingOnly=TRUE)
			file=args[1]
			genename=args[2]
			out1=args[3]
			out2=args[4]
			out3=args[5]
			setwd(dirname(file))
			#print(getwd())
				Sys.setenv('R_MAX_VSIZE'=32000000000)

				suppressMessages(suppressWarnings(library(msa)))
				suppressMessages(suppressWarnings(library(reshape2)))
				suppressMessages(suppressWarnings(library(seqinr)))
				 suppressMessages(suppressWarnings(library( Biostrings)))
				options(warn=-1)
				    library(stringr)
					##Functions
	
	

findPotentialStartsAndStops <- function(DNA_string)
  { positions =c()
   	types =c()
     if(length(positions)==0){
     codons            <- c('atg', 'taa', 'tag', 'tga','ATG', 'TAA', 'TAG', 'TGA')
     for (i in  1:length(codons))
     {
        codon <- codons[i]
        occurrences <- matchPattern(codon, DNA_string)
      codonpositions <- occurrences\@ranges\@start  
        numoccurrences <- length(codonpositions) 	
   	positions   <- append(positions,codonpositions, after=length(positions))
     types       <- append(types,rep(codon, numoccurrences), after=length(types))
 }
    indices <- order(positions)
     positions <- positions[indices]
     types <- types[indices]
     mylist <- list(positions,types)
 return(mylist)
 
     }
   }     

 
  findORFsinSeq <- function(sequence)
  {
     mylist <- findPotentialStartsAndStops(sequence)
     positions <- mylist[[1]]
     types <- mylist[[2]]
     orfstarts <- numeric()
     orfstops <- numeric()
     orflengths <- numeric()
     numpositions <- length(positions)
 
     if (numpositions >= 2)
     {
        for (i in 1:(numpositions-1))
        {
           posi <- positions[i]
           typei <- types[i]
           found <- 0
           while (found == 0)
           {
              for (j in (i+1):numpositions)
              {
                 posj  <- positions[j]
                 typej <- types[j]
                 posdiff <- posj - posi
                 posdiffmod3 <- posdiff %% 3
           
                 orflength <- posj - posi + 3
                 if ((typei == 'atg' || typei == 'ATG') && (typej == 'taa' || typej == 'tag' || typej == 'tga'||typej == 'TAA' || typej == 'TAG' || typej == 'TGA') && posdiffmod3 == 0)
                 {
                    numorfs <- length(orfstops)
                    usedstop <- -1
                    if (numorfs > 0)
                    {
                      for (k in 1:numorfs)
                      {
                          orfstopk <- orfstops[k]
                          if (orfstopk == (posj + 2)) { usedstop <- 1 }
                      }
                    }
                    if (usedstop == -1)
                    {
                       orfstarts <- append(orfstarts, posi, after=length(orfstarts))
                       orfstops <- append(orfstops, posj+2, after=length(orfstops)) 
                       orflengths <- append(orflengths, orflength, after=length(orflengths))
                    }
                    found <- 1
                    break
                 }
                 if (j == numpositions) { found <- 1 }
              }
           }
        }
     }
     indices <- order(orfstarts)
     orfstarts <- orfstarts[indices]
     orfstops <- orfstops[indices]
     orflengths <- numeric()
     numorfs <- length(orfstarts)
     for (i in 1:numorfs)
     {
        orfstart <- orfstarts[i]
        orfstop <- orfstops[i]
        orflength <- orfstop - orfstart + 1
        orflengths <- append(orflengths,orflength,after=length(orflengths))
     }
     mylist <- data.frame(orfstarts=orfstarts,orfstops= orfstops,orflengths= orflengths)
     return(mylist)
  }
  
			seq2Fasta <- function(seq,name,filename) 
				{
  						sink(filename,append=TRUE)
    					cat(paste0('>', name),file=filename,append=TRUE)
    					cat('\\n',file=filename,append=TRUE)
    					the.sequence <- toString(seq)
    					the.sequence=gsub('N','-',the.sequence)
    					cat(the.sequence,file=filename,append=TRUE)
    					cat('\\n',file=filename,append=TRUE)  
  						sink(NULL)
				}
				
						seqAA2Fasta <- function(seq,name,filename) 
				{
  						sink(filename,append=TRUE)
    					cat(paste0('>', name),file=filename,append=TRUE)
    					cat('\\n',file=filename,append=TRUE)
    					the.sequence <- toString(seq)
    					cat(the.sequence,file=filename,append=TRUE)
    					cat('\\n',file=filename,append=TRUE)  
  						sink(NULL)
				}

	
	
	Alignment <- function(x){
    	tryCatch(
        expr = {
 		return(invisible(msa(x)))
            message('Successfully executed')
       		 },
        error = function(e){
            message('Caught an error!')
		return(invisible(msa(x)))

        },
        warning = function(w){
           message('Caught an warning!')
            #print(w)
        },
        finally = {
            message('All done, quitting.')
      		  }
    		)    
		}


		mySequences <- readAAStringSet(file)
		namesunique=unique(mySequences\@ranges\@NAMES)
		mySequences <- mySequences[namesunique]
		
		myFirstAlignment=Alignment(mySequences)

		location=ifelse(grepl('plasmid',file),'plasmid','chromosome')


		nSeq=length(rownames(myFirstAlignment))
		#print(paste('number found',nSeq,sep=' '))
	
  if(nSeq==2){
  		#print('2 sequences')
      		
     X_Alignment=myFirstAlignment

			query=rownames(X_Alignment)[2]
      		Seq=toString(unmasked(X_Alignment)[[query]] )
      		Seq_string=strsplit(toString(Seq),'')
      		stops=which( Seq_string[[1]]=='-')
  		

			######
			ref=rownames(X_Alignment)[1]
      		Seqref=toString(unmasked(X_Alignment)[[ref]] )
      		Seqref_string=strsplit(toString(Seqref),'')
      		stopsref=which( Seqref_string[[1]]=='-')

######	
			Align=data.frame(Ref=c(unlist(Seqref_string)),Query=c(unlist(Seq_string)))
				Align\$position=c(1:dim(Align)[1])

			Align\$Ref=as.character(Align\$Ref)
			Align\$Query=as.character(Align\$Query)
			Align\$mutation=paste(Align\$Ref,Align\$position,Align\$Query,sep='')
			Align\$Change='identical'
			Align\$Change[ Align\$Query!=Align\$Ref ]='SNP'
			Align\$Change[ Align\$Query=='-' ]='GAP'
			Align\$Change[ Align\$Ref=='-' ]='INSERTION'
			gaps=length(which(Align\$Query=='-'))
				
			
			changesSNPs=Align[Align\$Change=='SNP',]
			changesGAPs=Align[Align\$Change=='GAP',]
			changesINSERTION=Align[Align\$Change=='INSERTION',]
			SNPs=paste(changesSNPs\$mutation,collapse=',')
			GAPs=paste(changesGAPs\$mutation,collapse=',')
			INSERTION=paste(changesINSERTION\$mutation,collapse=',')
			
			myAlignment=myFirstAlignment
			Ncopies=1
    }
    
      k=0
 	 if(nSeq==3){
 		#print('3 sequences')
 	 sequences=length(mySequences\@ranges\@NAMES)
 	  xx <- as.data.frame(matrix(ncol=sequences, nrow=myFirstAlignment\@unmasked\@ranges\@width[1]))
 
	  names=matrix(ncol=sequences, nrow=1)
	 
    		for( i in mySequences\@ranges\@NAMES){ 
        	x=c(unlist(strsplit(toString(unmasked(myFirstAlignment)[[i]]),'')))
        	k=k+1
    		xx[,k]=x
        	names[,k]=i
    		}
  
    		
    		xx\$consensus=apply(xx, 1, function(x) {
						y=unique(x[-1])	
    					y=y[y!='-']	
						if(length(y)==0){y='-'}
						if(length(y)>1){if(x[1]==x[2]){y=x[1]}}
						if(length(y)>1){if(x[1]==x[3]){y=x[1]}}
						if(length(y)>1){y='N'}
						return(y)})

		conSeq=paste(xx\$consensus,collapse='')
		conSeq=gsub('-','',conSeq)
  		consensusSeq=DNAStringSet(conSeq)
  		
  		consensusSeq\@ranges\@NAMES=mySequences\@ranges\@NAMES[2]
  		ref=DNAStringSet(mySequences[[1]])
  		ref\@ranges\@NAMES=mySequences\@ranges\@NAMES[1]
  		
  		myAlignment=invisible(Alignment(c(ref,consensusSeq)))

		
     X_Alignment=myAlignment

			query=rownames(X_Alignment)[2]
      		Seq=toString(unmasked(X_Alignment)[[query]] )
      		Seq_string=strsplit(toString(Seq),'')
      		stops=which( Seq_string[[1]]=='-')
  		

			######
			ref=rownames(X_Alignment)[1]
      		Seqref=toString(unmasked(X_Alignment)[[ref]] )
      		Seqref_string=strsplit(toString(Seqref),'')
      		stopsref=which( Seqref_string[[1]]=='-')

######	
			Align=data.frame(Ref=c(unlist(Seqref_string)),Query=c(unlist(Seq_string)))
				Align\$position=c(1:dim(Align)[1])

			Align\$Ref=as.character(Align\$Ref)
			Align\$Query=as.character(Align\$Query)
			Align\$mutation=paste(Align\$Ref,Align\$position,Align\$Query,sep='')
			Align\$Change='identical'
			Align\$Change[ Align\$Query!=Align\$Ref ]='SNP'
			Align\$Change[ Align\$Query=='-' ]='GAP'
			Align\$Change[ Align\$Ref=='-' ]='INSERTION'
			gaps=length(which(Align\$Query=='-'))
				
			
			changesSNPs=Align[Align\$Change=='SNP',]
			changesGAPs=Align[Align\$Change=='GAP',]
			changesINSERTION=Align[Align\$Change=='INSERTION',]
			SNPs=paste(changesSNPs\$mutation,collapse=',')
			GAPs=paste(changesGAPs\$mutation,collapse=',')
			INSERTION=paste(changesINSERTION\$mutation,collapse=',')
######	
		
			if(any(xx[,2]=='-') && any(xx[,3]=='-')){Ncopies=1}else{Ncopies=2}		
   
}




k=0
 if(nSeq==4){
 		#print('4 sequences')
 
  sequences=length(mySequences\@ranges\@NAMES)
 	  xx <- as.data.frame(matrix(ncol=sequences, nrow=myFirstAlignment\@unmasked\@ranges\@width[1]))
 
	  names=matrix(ncol=sequences, nrow=1)
	 
    		for( i in mySequences\@ranges\@NAMES){ 
        	x=c(unlist(strsplit(toString(unmasked(myFirstAlignment)[[i]]),'')))
        	k=k+1
    		xx[,k]=x
        	names[,k]=i
    		}

		xx\$consensus=apply(xx, 1, function(x) {
						y=unique(x[-1])	
    					y=y[y!='-']	
						if(length(y)==0){y='-'}
						if(length(y)>1){if(x[1]==x[2]){y=x[1]}}
						if(length(y)>1){if(x[1]==x[3]){y=x[1]}}
						if(length(y)>1){if(x[1]==x[4]){y=x[1]}}
						if(length(y)>1){y='N'}
						return(y)})

  		
  		conSeq=paste(xx\$consensus,collapse='')
		conSeq=gsub('-','',conSeq)
  		consensusSeq=DNAStringSet(conSeq)
  		consensusSeq\@ranges\@NAMES=mySequences\@ranges\@NAMES[2]
  		ref=DNAStringSet(mySequences[[1]])
  		ref\@ranges\@NAMES=mySequences\@ranges\@NAMES[1]
  		
  		myAlignment=invisible(Alignment(c(ref,consensusSeq)))
		

     X_Alignment=myAlignment


			query=rownames(X_Alignment)[2]
      		Seq=toString(unmasked(X_Alignment)[[query]] )
      		Seq_string=strsplit(toString(Seq),'')
      		stops=which( Seq_string[[1]]=='-')
  		

			######
			ref=rownames(X_Alignment)[1]
      		Seqref=toString(unmasked(X_Alignment)[[ref]] )
      		Seqref_string=strsplit(toString(Seqref),'')
      		stopsref=which( Seqref_string[[1]]=='-')

######	
			Align=data.frame(Ref=c(unlist(Seqref_string)),Query=c(unlist(Seq_string)))
				Align\$position=c(1:dim(Align)[1])

			Align\$Ref=as.character(Align\$Ref)
			Align\$Query=as.character(Align\$Query)
			Align\$mutation=paste(Align\$Ref,Align\$position,Align\$Query,sep='')
			Align\$Change='identical'
			Align\$Change[ Align\$Query!=Align\$Ref ]='SNP'
			Align\$Change[ Align\$Query=='-' ]='GAP'
			Align\$Change[ Align\$Ref=='-' ]='INSERTION'
			gaps=length(which(Align\$Query=='-'))
			
			
			changesSNPs=Align[Align\$Change=='SNP',]
			changesGAPs=Align[Align\$Change=='GAP',]
			changesINSERTION=Align[Align\$Change=='INSERTION',]
			SNPs=paste(changesSNPs\$mutation,collapse=',')
			GAPs=paste(changesGAPs\$mutation,collapse=',')
			INSERTION=paste(changesINSERTION\$mutation,collapse=',')
######	

			if(any(xx[,2]=='-') && any(xx[,3]=='-') && any(xx[,4]=='-')){Ncopies=1}else{Ncopies=3}		
			if(any(xx[,2]=='-') && any(xx[,3]=='-') ){Ncopies=1}else{Ncopies=2}		
			
}


k=0
 if(nSeq==5){
 		#print('5 sequences')
 
  sequences=length(mySequences\@ranges\@NAMES)
 	  xx <- as.data.frame(matrix(ncol=sequences, nrow=myFirstAlignment\@unmasked\@ranges\@width[1]))
 
	  names=matrix(ncol=sequences, nrow=1)
	 
    		for( i in mySequences\@ranges\@NAMES){ 
        	x=c(unlist(strsplit(toString(unmasked(myFirstAlignment)[[i]]),'')))
        	k=k+1
    		xx[,k]=x
        	names[,k]=i
    		}
    				xx\$consensus=apply(xx, 1, function(x) {
						y=unique(x[-1])	
    					y=y[y!='-']	
						if(length(y)==0){y='-'}
						if(length(y)>1){if(x[1]==x[2]){y=x[1]}}
						if(length(y)>1){if(x[1]==x[3]){y=x[1]}}
						if(length(y)>1){if(x[1]==x[4]){y=x[1]}}
						if(length(y)>1){if(x[1]==x[5]){y=x[1]}}
						if(length(y)>1){y='N'}
						return(y)})

  		
  		conSeq=paste(xx\$consensus,collapse='')
		conSeq=gsub('-','',conSeq)
  		consensusSeq=DNAStringSet(conSeq)
  		consensusSeq\@ranges\@NAMES=mySequences\@ranges\@NAMES[2]
  		
  		ref=DNAStringSet(mySequences[[1]])
  		ref\@ranges\@NAMES=mySequences\@ranges\@NAMES[1]
  		
  		myAlignment=invisible(Alignment(c(ref,consensusSeq)))
		     X_Alignment=myAlignment

			query=rownames(X_Alignment)[2]
      		Seq=toString(unmasked(X_Alignment)[[query]] )
      		Seq_string=strsplit(toString(Seq),'')
      		stops=which( Seq_string[[1]]=='-')
  		

			######
			ref=rownames(X_Alignment)[1]
      		Seqref=toString(unmasked(X_Alignment)[[ref]] )
      		Seqref_string=strsplit(toString(Seqref),'')
      		stopsref=which( Seqref_string[[1]]=='-')

######	
			Align=data.frame(Ref=c(unlist(Seqref_string)),Query=c(unlist(Seq_string)))
				Align\$position=c(1:dim(Align)[1])

			Align\$Ref=as.character(Align\$Ref)
			Align\$Query=as.character(Align\$Query)
			Align\$mutation=paste(Align\$Ref,Align\$position,Align\$Query,sep='')
			Align\$Change='identical'
			Align\$Change[ Align\$Query!=Align\$Ref ]='SNP'
			Align\$Change[ Align\$Query=='-' ]='GAP'
			Align\$Change[ Align\$Ref=='-' ]='INSERTION'
			gaps=length(which(Align\$Query=='-'))
			
			
			changesSNPs=Align[Align\$Change=='SNP',]
			changesGAPs=Align[Align\$Change=='GAP',]
			changesINSERTION=Align[Align\$Change=='INSERTION',]
			SNPs=paste(changesSNPs\$mutation,collapse=',')
			GAPs=paste(changesGAPs\$mutation,collapse=',')
			INSERTION=paste(changesINSERTION\$mutation,collapse=',')
######	
			if(any(xx[,2]=='-') && any(xx[,3]=='-') && any(xx[,4]=='-') && any(xx[,5]=='-')){Ncopies=1}else{Ncopies=4}	
			if(any(xx[,2]=='-') && any(xx[,3]=='-') && any(xx[,4]=='-' )){Ncopies=1}else{Ncopies=3}		
			if(any(xx[,2]=='-') && any(xx[,3]=='-') ){Ncopies=1}else{Ncopies=2}		

}

k=0
 if(nSeq==6){
 		#print('6 sequences')
 
  sequences=length(mySequences\@ranges\@NAMES)
 	  xx <- as.data.frame(matrix(ncol=sequences, nrow=myFirstAlignment\@unmasked\@ranges\@width[1]))
 
	  names=matrix(ncol=sequences, nrow=1)
	 
    		for( i in mySequences\@ranges\@NAMES){ 
        	x=c(unlist(strsplit(toString(unmasked(myFirstAlignment)[[i]]),'')))
        	k=k+1
    		xx[,k]=x
        	names[,k]=i
    		}
    				xx\$consensus=apply(xx, 1, function(x) {
						y=unique(x[-1])	
    					y=y[y!='-']	
						if(length(y)==0){y='-'}
						if(length(y)>1){if(x[1]==x[2]){y=x[1]}}
						if(length(y)>1){if(x[1]==x[3]){y=x[1]}}
						if(length(y)>1){if(x[1]==x[4]){y=x[1]}}
						if(length(y)>1){if(x[1]==x[5]){y=x[1]}}
						if(length(y)>1){if(x[1]==x[6]){y=x[1]}}
						if(length(y)>1){y='N'}
						return(y)})

  		
  		conSeq=paste(xx\$consensus,collapse='')
		conSeq=gsub('-','',conSeq)
  		consensusSeq=DNAStringSet(conSeq)
  		consensusSeq\@ranges\@NAMES=mySequences\@ranges\@NAMES[2]
  		ref=DNAStringSet(mySequences[[1]])
  		ref\@ranges\@NAMES=mySequences\@ranges\@NAMES[1]
  		
  		myAlignment=invisible(Alignment(c(ref,consensusSeq)))
		     X_Alignment=myAlignment

			query=rownames(X_Alignment)[2]
      		Seq=toString(unmasked(X_Alignment)[[query]] )
      		Seq_string=strsplit(toString(Seq),'')
      		stops=which( Seq_string[[1]]=='-')
  		

			######
			ref=rownames(X_Alignment)[1]
      		Seqref=toString(unmasked(X_Alignment)[[ref]] )
      		Seqref_string=strsplit(toString(Seqref),'')
      		stopsref=which( Seqref_string[[1]]=='-')

######	
			Align=data.frame(Ref=c(unlist(Seqref_string)),Query=c(unlist(Seq_string)))
				Align\$position=c(1:dim(Align)[1])

			Align\$Ref=as.character(Align\$Ref)
			Align\$Query=as.character(Align\$Query)
			Align\$mutation=paste(Align\$Ref,Align\$position,Align\$Query,sep='')
			Align\$Change='identical'
			Align\$Change[ Align\$Query!=Align\$Ref ]='SNP'
			Align\$Change[ Align\$Query=='-' ]='GAP'
			Align\$Change[ Align\$Ref=='-' ]='INSERTION'
			gaps=length(which(Align\$Query=='-'))

			
			changesSNPs=Align[Align\$Change=='SNP',]
			changesGAPs=Align[Align\$Change=='GAP',]
			changesINSERTION=Align[Align\$Change=='INSERTION',]
			SNPs=paste(changesSNPs\$mutation,collapse=',')
			GAPs=paste(changesGAPs\$mutation,collapse=',')
			INSERTION=paste(changesINSERTION\$mutation,collapse=',')
######	
			
			if(any(xx[,2]=='-') && any(xx[,3]=='-') && any(xx[,4]=='-') && any(xx[,5]=='-') && any(xx[,6]=='-')){Ncopies=1}else{Ncopies=5}
			if(any(xx[,2]=='-') && any(xx[,3]=='-') && any(xx[,4]=='-') && any(xx[,5]=='-')){Ncopies=1}else{Ncopies=4}	
			if(any(xx[,2]=='-') && any(xx[,3]=='-') && any(xx[,4]=='-' )){Ncopies=1}else{Ncopies=3}		
			if(any(xx[,2]=='-') && any(xx[,3]=='-') ){Ncopies=1}else{Ncopies=2}		

}

k=0
 if(nSeq==7){
 		#print('7 sequences')
 
  sequences=length(mySequences\@ranges\@NAMES)
 	  xx <- as.data.frame(matrix(ncol=sequences, nrow=myFirstAlignment\@unmasked\@ranges\@width[1]))
 
	  names=matrix(ncol=sequences, nrow=1)
	 
    		for( i in mySequences\@ranges\@NAMES){ 
        	x=c(unlist(strsplit(toString(unmasked(myFirstAlignment)[[i]]),'')))
        	k=k+1
    		xx[,k]=x
        	names[,k]=i
    		}
    				xx\$consensus=apply(xx, 1, function(x) {
						y=unique(x[-1])	
    					y=y[y!='-']	
						if(length(y)==0){y='-'}
						if(length(y)>1){if(x[1]==x[2]){y=x[1]}}
						if(length(y)>1){if(x[1]==x[3]){y=x[1]}}
						if(length(y)>1){if(x[1]==x[4]){y=x[1]}}
						if(length(y)>1){if(x[1]==x[5]){y=x[1]}}
						if(length(y)>1){if(x[1]==x[6]){y=x[1]}}
						if(length(y)>1){if(x[1]==x[7]){y=x[1]}}
						if(length(y)>1){y='N'}
						return(y)})

  		
  		conSeq=paste(xx\$consensus,collapse='')
		conSeq=gsub('-','',conSeq)
  		consensusSeq=DNAStringSet(conSeq)
  		consensusSeq\@ranges\@NAMES=mySequences\@ranges\@NAMES[2]
  		ref=DNAStringSet(mySequences[[1]])
  		ref\@ranges\@NAMES=mySequences\@ranges\@NAMES[1]
  		
  		myAlignment=invisible(Alignment(c(ref,consensusSeq)))
		     X_Alignment=myAlignment

			query=rownames(X_Alignment)[2]
      		Seq=toString(unmasked(X_Alignment)[[query]] )
      		Seq_string=strsplit(toString(Seq),'')
      		stops=which( Seq_string[[1]]=='-')
  		

			######
			ref=rownames(X_Alignment)[1]
      		Seqref=toString(unmasked(X_Alignment)[[ref]] )
      		Seqref_string=strsplit(toString(Seqref),'')
      		stopsref=which( Seqref_string[[1]]=='-')
      		
######	
			Align=data.frame(Ref=c(unlist(Seqref_string)),Query=c(unlist(Seq_string)))
				Align\$position=c(1:dim(Align)[1])

			Align\$Ref=as.character(Align\$Ref)
			Align\$Query=as.character(Align\$Query)
			Align\$mutation=paste(Align\$Ref,Align\$position,Align\$Query,sep='')
			Align\$Change='identical'
			Align\$Change[ Align\$Query!=Align\$Ref ]='SNP'
			Align\$Change[ Align\$Query=='-' ]='GAP'
			Align\$Change[ Align\$Ref=='-' ]='INSERTION'
			gaps=length(which(Align\$Query=='-'))
			
			changesSNPs=Align[Align\$Change=='SNP',]
			changesGAPs=Align[Align\$Change=='GAP',]
			changesINSERTION=Align[Align\$Change=='INSERTION',]
			SNPs=paste(changesSNPs\$mutation,collapse=',')
			GAPs=paste(changesGAPs\$mutation,collapse=',')
			INSERTION=paste(changesINSERTION\$mutation,collapse=',')
######	
		
			if(any(xx[,2]=='-') && any(xx[,3]=='-') && any(xx[,4]=='-') && any(xx[,5]=='-') && any(xx[,6]=='-') && any(xx[,6]=='-')){Ncopies=1}else{Ncopies=6}
			if(any(xx[,2]=='-') && any(xx[,3]=='-') && any(xx[,4]=='-') && any(xx[,5]=='-') && any(xx[,6]=='-')){Ncopies=1}else{Ncopies=5}
			if(any(xx[,2]=='-') && any(xx[,3]=='-') && any(xx[,4]=='-') && any(xx[,5]=='-')){Ncopies=1}else{Ncopies=4}	
			if(any(xx[,2]=='-') && any(xx[,3]=='-') && any(xx[,4]=='-' )){Ncopies=1}else{Ncopies=3}		
			if(any(xx[,2]=='-') && any(xx[,3]=='-') ){Ncopies=1}else{Ncopies=2}		

}


names=mySequences\@ranges\@NAMES

if(Ncopies==1){

myAlignmentAA=invisible(msa(c(AAStringSet(Seqref),AAStringSet(Seq))))
RefAA=unlist(strsplit(toString(  toString(unmasked(myAlignmentAA)[[1]] ) ),''))
SeqAA=unlist(strsplit(toString(  toString(unmasked(myAlignmentAA)[[2]] ) ),''))
PepAlign=data.frame(Ref=unlist(RefAA),Query=unlist(SeqAA))
 
PepAlign\$match[PepAlign\$Ref== PepAlign\$Query]=1
PepAlign\$match[PepAlign\$Ref!= PepAlign\$Query]=0

PepAlign\$match[PepAlign\$Ref== PepAlign\$Query]=1
PepAlign\$match[PepAlign\$Ref!= PepAlign\$Query]=0

PepAlignSimilarity=sum(PepAlign\$match/length(PepAlign\$match))*100	

query=strsplit(names[2],'###')[[1]][1]
start=strsplit(names[2],'###')[[1]][2]
end=strsplit(names[2],'###')[[1]][3]
strand=strsplit(names[2],'###')[[1]][4]
GenomeFile=strsplit(names[2],':')[[1]][2]


startingcodon= min(which(PepAlign\$match==1))

if (stopsref != stops[1]){classification='possibly_truncated'}else{classification='present'}
if (startingcodon != 1){classification='late_startcodon'}else{classification='possibly_truncated'}


   		metadatasummary=data.frame(
      		gene=genename,
      			Ref=paste('\"',ref,'\"',sep=''),
 					Ref_length_nt=length(SeqAA)*3,
 					Ref_length_AA=length(RefAA),
 					Ref_stopCodon=stopsref,
 					GenomeFile=GenomeFile,
 					Query=paste('\"',query,'\"',sep=''),
 					Start=start,
 					End=end,
 					Strand= strand,
 					Query_length_nt=NA,
 					Query_length_AA=length(SeqAA),
 					nt_gaps=gaps,
 					stopCodon= stops[1],
 					SNPs=SNPs,
 					GAPs=GAPs,
 					INSERTION=INSERTION,
 					location=location,
 					classification=classification,
 					PepAlignSimilarity=PepAlignSimilarity,
 					DNAAlignSimilarity=NA,
 						Ncopies=Ncopies)

SeqQueryGenomeName=paste(GenomeFile,query,sep='::')
      		seq2Fasta(Seq,SeqQueryGenomeName,out1)
       		seqAA2Fasta(Seq,SeqQueryGenomeName, out2)
  write.table(metadatasummary,file=out3,append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t')
}

count=1

if(Ncopies>1){

names=names[-1]

for (i in names){

count=count+1

query=strsplit(i,'###')[[1]][1]
start=strsplit(i,'###')[[1]][2]
end=strsplit(i,'###')[[1]][3]
strand=strsplit(i,'###')[[1]][4]
GenomeFile=strsplit(names[1],':')[[1]][2]

SeqName=xx[count]
ff=gsub('\"','',SeqName)
ff=gsub('\\n','',ff)
ff=gsub(', ','',ff)
ff=gsub('c\\\\(','',ff)
ff=gsub('\\\\)','',ff)

SeqName=ff
stops=which( Seq_trl_string[[1]]=='-')


myAlignmentAA=invisible(msa(c(AAStringSet(Seqref),AAStringSet(SeqName))))
RefAA=strsplit(toString(  toString(unmasked(myAlignmentAA)[[1]] ) ),'')
SeqAA=strsplit(toString(  toString(unmasked(myAlignmentAA)[[2]] ) ),'')
PepAlign=data.frame(Ref=unlist(RefAA),Query=unlist(SeqAA))
 
 
PepAlign\$match[PepAlign\$Ref== PepAlign\$Query]=1
PepAlign\$match[PepAlign\$Ref!= PepAlign\$Query]=0
 
PepAlign\$match[PepAlign\$Ref== PepAlign\$Query]=1
PepAlign\$match[PepAlign\$Ref!= PepAlign\$Query]=0


PepAlignSimilarity=sum(PepAlign\$match/length(PepAlign\$match))*100	

startingcodon= min(which(PepAlign\$match==1))

if (stopsref != stops[1]){classification='possibly_truncated'}else{classification='present'}
if (startingcodon != 1){classification='late_startcodon'}else{classification='present'}

  		metadatasummary=data.frame(
      		gene=genename,
      			Ref=paste('\"',ref,'\"',sep=''),
 					Ref_length_nt=length(SeqAA)*3,
 					Ref_length_AA=length(RefAA),
 					Ref_stopCodon=stopsref,
 					GenomeFile=GenomeFile,
 					Query=paste('\"',query,'\"',sep=''),
 					Start=start,
 					End=end,
 					Strand= strand,
 					Query_length_nt=NA,
 					Query_length_AA=length(SeqAA),
 					nt_gaps=gaps,
 					stopCodon= stops[1],
 					SNPs=SNPs,
 					GAPs=GAPs,
 					INSERTION=INSERTION,
 					location=location,
 					classification=classification,
 					PepAlignSimilarity=PepAlignSimilarity,
 					DNAAlignSimilarity=NA,
 						Ncopies=Ncopies)

SeqQueryGenomeName=paste(GenomeFile,query,sep='::')
      		seq2Fasta(SeqName,SeqQueryGenomeName,out1)
       		seqAA2Fasta(SeqName,SeqQueryGenomeName, out2)
  write.table(metadatasummary,file=out3,append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t')
}


#print (genename)
#print (out3)
}
";
 print R_SCRIPT $R_script;
        close R_SCRIPT;
        system("chmod a+x $R_file"); 
}



























sub sub7{
my $gene= shift;
	my $genename=basename($gene);
	$genename=~s/.fasta$//g;
	my $outAlign=$dirresults."/".$genename.".nt_alignment.fasta";
	my $outAlign_forR=$dirresults."/".$genename.".nt_alignment.fasta_forR";
	my $outAlign_forRR=$dirresults."/".$genename.".nt_alignment.fasta_forR_forR";


	my $cmd_align="cat ".$dirresults."/".$genename."*nt_alignment.fasta.tmp >$outAlign";
	#print ">>>>>>>>>>>>>$cmd_align\n";
	system($cmd_align);

	
	
	
	my $cmd_sed="sed 's/-//g' $outAlign >$outAlign_forRR";
	my $cmd_sed1="sed 's/:/-/g' $outAlign_forRR >$outAlign_forR";

	system($cmd_sed);
	system($cmd_sed1);
	
	$R_file=$outAlign."temp.R";
	
		open R_SCRIPT,">$R_file" or die "Cannot write $R_file script\n";

         $R_script= "
		rm(list=ls()); 
					
		Sys.setenv('R_MAX_VSIZE'=32000000000)

				suppressMessages(suppressWarnings(library(msa)))
				suppressMessages(suppressWarnings(library(reshape2)))
				options(warn=-1)
				library(stringr)

					##Functions
	
			seq2Fasta <- function(seq,name,filename) 
				{
  						sink(filename,append=TRUE)
    					cat(paste0('>', name),file=filename,append=TRUE)
    					cat('\\n',file=filename,append=TRUE)
    					the.sequence <- toString(seq)
    					the.sequence=gsub('N','-',the.sequence)
    					cat(the.sequence,file=filename,append=TRUE)
    					cat('\\n',file=filename,append=TRUE)  
  						sink(NULL)
				}
				
						seqAA2Fasta <- function(seq,name,filename) 
				{
  						sink(filename,append=TRUE)
    					cat(paste0('>', name),file=filename,append=TRUE)
    					cat('\\n',file=filename,append=TRUE)
    					the.sequence <- toString(seq)
    					cat(the.sequence,file=filename,append=TRUE)
    					cat('\\n',file=filename,append=TRUE)  
  						sink(NULL)
				}

		file=\"$outAlign_forR\"
		#setwd(dirname(file))
		############
		mySequences <- readDNAStringSet(file)
nSeq=length(mySequences\@ranges\@NAMES)

Alignment <- function(x){
    	tryCatch(
        expr = {
 		return(invisible(msa(x,gapExtension =150,gap=1500,method='Muscle')))
            message('Successfully executed')
       		 },
        error = function(e){
          #  message('Caught an error!')
		return(invisible(  msa(x)))

        },
        warning = function(w){
            message('Caught an warning!')
            print(w)
        },
        finally = {
            message('All done, quitting.')
      		  }
    		)    
		}


if(nSeq>1){


		myFirstAlignment=Alignment(c(mySequences[1],mySequences[2]))
		myFirstAlignment_reverse=Alignment(c(reverseComplement(mySequences[1]),mySequences[2]))
		

 		FA=msaConsensusSequence(myFirstAlignment)
 		FAR=msaConsensusSequence(myFirstAlignment_reverse)

		Consensus=str_count(FA, pattern = '\\\\?')
 		Consensus_reverse=str_count(FAR, pattern ='\\\\?')
	
		if(Consensus<Consensus_reverse){myFirstAlignment=msa(mySequences)}else{myFirstAlignment=msa(c(reverseComplement(mySequences[1]),mySequences[-1]))}


  fasta_file =\"$outAlign\"
 
  msaPrettyPrint(myFirstAlignment, output='tex', showConsensus = 'none', askForOverwrite=FALSE, verbose=FALSE,
                   alFile = fasta_file )

}

";

	
        print R_SCRIPT $R_script;
        close R_SCRIPT;
        system("chmod a+x $R_file"); 
        #print "Start $R_file\n";
		system("Rscript $R_file --slave");
		#print " END $R_file\n";
 		system("rm $R_file");
 		system("rm $outAlign_forR");
 		system("rm $outAlign_forRR");
 			
	
	}
	



###########################PROGRAM RUNS HERE###############################

print "
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

\t\t CREATING DATABASES FROM GENES OF INTEREST

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
  
 afork(@genes,$cores,\&subgene);
  
print "
-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-

\t\t DONE CREATING DATABASES FROM GENES OF INTEREST
  
-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-\n\n";


if ($species eq "SI"){ 
 
 
   print "
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

\t\t PREDICTING PLASMID AND CHROMOSOMES FOR ALL GENOMES

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
 
 
 
 afork(@seq,$cores,\&sub0);
 
 print "
-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-

\t\t DONE PREDICTING PLASMID AND CHROMOSOMES FOR ALL GENOMES\n
  
-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-\n\n";
 
 
 }else{ 

  print "
-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-

\t\t NO PREDICTION OF PLASMIDS WILL BE PERFORMED\n
  
-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-\n\n";
 
 }




 if ($species eq "NO"){ 
 
 afork(@seq,$cores,\&sub1_chr)}else{ 

 afork(@seq,$cores,\&sub1);
 
 }





 if ($species eq "NO"){ 
   print "
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

\t\t GENES PREDICTION FOR CHROMOSOME 

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
 opendir(BIN, $diroriginalgenomes) or die "Can't open $diroriginalgenomes: $!";
 @seq = grep { -T "$diroriginalgenomes/$_" } readdir BIN;
  afork(@seq,$cores,\&sub2_chr);
    print "
-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-

\t\t DONE GENES PREDICTION FOR CHROMOSOME 
  
-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-\n\n";
  
  
 }else{
   print "
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

\t\t EXTRACTING GENES FROM GENOMES 

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
 
 
 afork(@seq,$cores,\&sub2);
 
 
 print "
-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-

\t\t DONE EXTRACTING GENES FROM GENOMES 
  
-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-\n\n";
 
 }
#normal(@seq,\&sub2);




 if ($species eq "NO"){ 

  opendir(BIN, $diroriginalgenomes) or die "Can't open $diroriginalgenomes: $!";
 @seq = grep { -T "$diroriginalgenomes/$_" } readdir BIN;
 
 afork(@seq,$cores,\&sub3_chr);
 

 
 }else{
 
   print "
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

\t\t GENES PREDICTION FOR CHROMOSOME AND PLASMIDS  

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
afork(@seq,$cores,\&sub3);

 print "
-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-

\t\t DONE GENES PREDICTION FOR CHROMOSOME AND PLASMIDS 
  
-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-\n\n";
 }



    print "
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

\t\t MODIFICATION OF PRODIGAL OUTPUTS

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
 if ($species eq "NO"){ 
  opendir(BIN, $diroriginalgenomes) or die "Can't open $diroriginalgenomes: $!";
 @seq = grep { -T "$diroriginalgenomes/$_" } readdir BIN;
 afork(@seq,$cores,\&sub4_chr)}else{afork(@seq,$cores,\&sub4);}
  

print "
-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-

\t\t DONE  MODIFICATION OF PRODIGAL OUTPUTS
  
-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-\n\n";

 print "
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

\t\t BLASTING GENES 


>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
 
 if ($species eq "NO"){ opendir(BIN, $diroriginalgenomes) or die "Can't open $diroriginalgenomes: $!";
 @seq = grep { -T "$diroriginalgenomes/$_" } readdir BIN;
 afork(@seq,$cores,\&sub5_chr);}else{afork(@seq,$cores,\&sub5);}
 
 if ($species eq "NO"){ opendir(BIN, $diroriginalgenomes) or die "Can't open $diroriginalgenomes: $!";
 @seq = grep { -T "$diroriginalgenomes/$_" } readdir BIN;
 afork(@seq,$cores,\&sub5_1chr);}else{afork(@seq,$cores,\&sub5_1);}
 

print "
-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-

\t\t DONE BLASTING GENES  
  
-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-\n\n";

   print "
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

\t\t GENERATING CURATED ALIGNMENTS AND DESCRIPTIVE TABLE

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";

#normal(@seq,\&sub6);
if ($species eq "NO"){ opendir(BIN, $diroriginalgenomes) or die "Can't open $diroriginalgenomes: $!";
 @seq = grep { -T "$diroriginalgenomes/$_" } readdir BIN;
 afork(@seq,$cores,\&sub6_chr);}else{afork(@seq,$cores,\&sub6)}

print "
-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-

\t\t DONE GENERATING CURATED ALIGNMENTS AND DESCRIPTIVE TABLE  
  
-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-\n\n";


###########################CLEAN UP AND COMPRESS###############################
#normal(@seq);
	
	
	  print "
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

\t\t CONCATENATION OF ALIGMENTS 

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";

 afork(@genes,$cores,\&sub7);

print "
-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-

\t\t DONE CONCATENATION OF ALIGNMENTS 
  
-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-\n\n";

	  print "
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

\t\t CONCATENATION TABLE

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";



my $prefix=basename($out);
my $cmd_cat="for f in $dirresults/*description; do cat \$f >> $dirresults/$prefix.Allgenes_descriptive_table.temp ; rm \$f; done";

system($cmd_cat);

system("cat $dirresults/out.alignments.description.txt $dirresults/$prefix.Allgenes_descriptive_table.temp >$dirresults/$prefix.Allgenes_descriptive_table.txt"); 
system("rm $dirresults/$prefix.Allgenes_descriptive_table.temp");
system("rm $dirresults/*description*");
system("rm $dirresults/*tmp");
	

#system("cat $dirresults/out.alignments.description.txt $dirresults/*description > $dirresults/$prefix.Allgenes_descriptive_table.txt");
#system("rm $dirresults/out.alignments.description.txt");
#system("rm $dirresults/*description*");

print "
-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-

\t\t DONE CONCATENATION OF TABLE
  
-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-\n\n";

system("find $diroutput/*/. -size 0 |  xargs rm");
system("find  $dirresultspep/. -size 0 |  xargs rm");





print "Process completed\n";

print "
-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-\n";

print "Find intermedia files in: $output_folder\n";
print "Find final files in: $results_folder\n";
print "
-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-\n";

 	# Do stuff

system($cmd_fai);
system("rm $diroutput/*/*fai");

 
 
 my $duration = time - $start;
print "Execution time: $duration s\n";

	

 
 
 



























