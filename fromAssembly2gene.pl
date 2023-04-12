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


## run like ./fromAssembly2gene.pl -g gene/*fasta -a genome/* 
##debug ./fromAssembly2gene.pl -a test_assembly/*fasta -g test_genes/AA/*fasta -o test1_AA -c 1 -p -t prot
##debug ./fromAssembly2gene.pl -a test_assembly/*fasta -g test_genes/NT/*fasta -o test1_NT -c 1 -p 
#######################################################################################################################################################################################################################################################################################################################################################################
#global variables
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
my $plasmid_pred;
my $genes_type;
my $gene_out_alignment_fasta=();
#######################################################################################################################################################################################################################################################################################################################################################################
#arg options

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

# check flags 

if ($help){
&usage();
 exit(1); 
 }

if ($version){
print STDERR "\nVersion 2.22. Program Last updated 12th April 2023 \n\n";
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
 print STDERR "\nERROR: No output prefix given. This is required\n\n";
 exit(1); 
 }

# defaults

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

if (defined $plasmid ){
	if($plasmid){ $plasmid_pred="SI"};
}else{ 
$plasmid_pred="NO";
}


# help 

sub usage{
 print STDERR <<EOF;
 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
																																	
 From Assemblies to genes (Program Last updated 12th April 2023)																									
 																																	
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

#######################################################################################################################################################################################################################################################################################################################################################################

# directories paths

 #input folders
 my $rel_path_genes = $genes[0];
 my $abs_path_genes= File::Spec->rel2abs( $rel_path_genes ) ;
 my $dir_genes=dirname($abs_path_genes);
 my $rel_path_genomes = $seq[0];
 my $abs_path_genomes= File::Spec->rel2abs( $rel_path_genomes ) ;
 my $dir_genomes=dirname($abs_path_genomes);

 #output folders
 my $dir= getcwd();
 my $output_folder=$out."_output";
 my $results_folder=$out."_results";
 my $diroutput= $dir."/".$output_folder;
 my $dirresults= $dir."/".$results_folder;
 my $dirresultspep= $dirresults."/Peptides";
 my $dirModifiedGenes= $diroutput."/ModifiedGenes";
 my $dirModifiedGenomes= $diroutput."/ModifiedGenomes";
 my $dirscript= $dirresults."/scripts";

# check the folders dont exists
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

# create ouput folders

system("mkdir $diroutput");
system("mkdir $dirresults");
system("mkdir $dirresultspep");
system("mkdir $dirscript");
system("mkdir $dirModifiedGenes");
system("mkdir $dirModifiedGenomes");

print "
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";

print "folder for intermedia files: $output_folder\n";
print "folder for results files: $results_folder\n";
print "folder for modified inputs files: ModifiedGenes\n";
print "folder for modified inputs files: ModifiedGenomes\n";
print "
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n\n";

print "
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
print "\n\nPLEASE NOTE: Copies of the original files are stored in *output/Modified.. folders and some format modification are going to be performed on them.\n\n";
print "
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";


#######################################################################################################################################################################################################################################################################################################################################################################
#functions

sub check_exists_command { 
    my $check = `sh -c 'command -v $_[0]'`; 
    return $check;
}

sub sub1_genemodification{
	my $genes_name = shift;
	my $genes=$dirModifiedGenes."/".$genes_name;
	#print "$genes\n\n";
	unless (-e $genes | -s $genes){ die "ERROR: $genes file is empty or doesn't exist";}
	#check what type of file is input (different OS)
	my $typefile=`file $genes`;
	my $cmd_gene;
	my $linegene=$genes."tmp";
	if ( $typefile=~m/: ASCII*+CRLF/g){
		my $cmd_gene="awk '/^>/{gsub(\"\\r\",\"\\n\");printf(\"\%s\",\$0);next; }{gsub(\"\\r\",\"\");printf(\"\%s\",\$0);}  END {printf(\"\\n\");} ' < $genes | awk '/^>/{printf(\"\\n\%s\\n\",\$0);next; }{printf(\"\%s\",\$0);}  END {printf(\"\\n\");}' >$linegene";
		print "------CRLF $cmd_gene-------\n\n";
		system($cmd_gene);
		}elsif($typefile=~m/: ASCII*+CR/g){
		my $cmd_gene="awk '/^>/{gsub(\"\\r\",\"\\n\");printf(\"\%s\",\$0);next;}{gsub(\"\$\",\"\");printf(\"\%s\",\$0);}  END {printf(\"\\n\");} ' < $genes | awk '/^>/{printf(\"\\n\%s\\n\",\$0);next; }{printf(\"\%s\",\$0);}  END {printf(\"\\n\");}'  >$linegene";
		print "------CR $cmd_gene-------\n\n";
		system($cmd_gene);
		}else{ 
			my $cmd_gene="awk '/^>/{printf(\"\\n\%s\\n\",\$0);next; }{printf(\"\%s\",\$0);}  END {printf(\"\\n\");} ' < $genes >$linegene";
			print "------ELSE $cmd_gene-------\n\n";
			system($cmd_gene);
			}
	#check if gene is in nucleotide or in aminoacid
	# add last codon
	if( $genes_type eq "nucl"){
		my $cmd_grep=`grep -o ^M $linegene`;
		#print "\n---$cmd_grep---\n\n";
		if($cmd_grep ne ""){ die  print "ERROR: Gene $genes looks like an aminoacid sequence. The input genes are considered nucleotide sequences by default. please set -t 'protein' if your sequences are aminoacids. Do not mix sequences types" }
		my $lastcodon=`grep -o '..[^ ]\$' $linegene|tail -1`;
		chomp($lastcodon);
		$lastcodon =~ s/(.*?)\s?[\n?|\r?]/$1/;
		my $dec;
		my @array=("TAG","TAA","TGA","tag","taa","tga");
		my %params = map { $_ => 1 } @array;
		if(exists($params{$lastcodon}) | $genes_type eq 'prot') { system("echo '\n' >> $linegene")}elsif ( $genes_type eq 'nucl'){system("echo 'TAG\n' >> $linegene")} 
	}
	system("mv $linegene $genes");
	my $gene_tag=$genes_name;
   	$gene_tag=~s/.fasta$//g;
	$gene_out_alignment_fasta=$dirresults."/".$gene_tag."_alignment.fasta";
	system("cat $genes >$gene_out_alignment_fasta.tmp");
	#print "\n\n after processing $genes\n\n";
	unless (-e $genes | -s $genes){ die "ERROR: something went wrong while reformating $genes please check file";}
   	}

#was within subgene
sub sub2_gene_makeblastdb{
	my $genes_name = shift;
	my $genes=$dirModifiedGenes."/".$genes_name;
	#print "$genes\n\n";
   	my $gene_tag=$dirModifiedGenes."/".$genes_name;
   	$gene_tag=~s/.fasta$//g;
	print "\t\tRunning makeblastdb on $genes\n";
 	my $cmd_makeblastdb= "makeblastdb -in ".$genes." -dbtype ".$genes_type." -out ".$gene_tag;
  	print "You are asking to run\n $cmd_makeblastdb\n";
  	system($cmd_makeblastdb);
	if($genes_type eq "prot"){ 
		my $db_blast=$gene_tag.".pin";
		unless (-e $db_blast | -s $db_blast){ die "ERROR: something went wrong while running makeblastdb $genes please check file";}
	}
	if($genes_type eq "nucl"){ 
		my $db_blast=$gene_tag.".nin";
		unless (-e $db_blast | -s $db_blast){ die "ERROR: something went wrong while running makeblastdb $genes please check file";}
	}
}

sub sub3_chrplas_folders{
	my $seq = shift;
	$seq=$dirModifiedGenomes."/".$seq;
	my $genome=basename($seq);	
	my $dirprochr=$diroutput."/".$genome.".chromosome";
	my $fasta_chr=$dirprochr."/".$genome.".chromosome.fasta";
	my $diralignments=$dirprochr."/alignments";
	my $cmd_prod_chr="mkdir $dirprochr";
	my $cmd_blast_chr="mkdir $diralignments";
	system("$cmd_prod_chr");
	system("$cmd_blast_chr");
	if($plasmid_pred eq "NO"){
		my $cmd_cp="cp $seq $fasta_chr";
		system($cmd_cp);
		}
	if($plasmid_pred eq "SI"){
		my $dirproplas=$diroutput."/".$genome.".plasmid";
		my $diralignments=$dirproplas."/alignments";
		my $cmd_prod_plasmid="mkdir $dirproplas";
		my $cmd_blast_plasmid="mkdir $diralignments";
		system("$cmd_prod_plasmid");
		system("$cmd_blast_plasmid");
	}
 }
 

#was sub 0
sub sub4_plasclass_method{
	my $seq = shift;
	$seq=$dirModifiedGenomes."/".$seq;
	my $genome=basename($seq);
	my $dirproplas=$diroutput."/".$genome.".plasmid";	
	my $dirprochr=$diroutput."/".$genome.".chromosome";
	my $plaspred=$diroutput."/".$genome."_plas.prob.out";
	system("classify_fasta.py -f $seq -o $plaspred");
 }
 
sub sub5_plasmid_extraction{
	my $seq = shift;
	$seq=$dirModifiedGenomes."/".$seq;
	my $genome=basename($seq);
	my $dirproplas=$diroutput."/".$genome.".plasmid";	
	my $dirprochr=$diroutput."/".$genome.".chromosome";
	my $plaspred=$diroutput."/".$genome."_plas.prob.out";
    system("Rscript $results_folder/scripts/plasmid_chromosome_Detection.r $plaspred $dirprochr $dirproplas --slave");
 }


sub var_inputs{
	my $seq = shift;
	$seq=$dirModifiedGenomes."/".$seq;
    my $genome=basename($seq);
	my $dirproplas=$diroutput."/".$genome.".plasmid";	
	my $dirprochr=$diroutput."/".$genome.".chromosome";
	my $input_plasmid=$dirproplas."/".$genome.".plasmid_contigslist.txt";
	my $input_chr=$dirprochr."/".$genome.".chromosome_contigslist.txt";
	my $output_plasmid=$input_plasmid;
	my $output_chr=$input_chr;
	my $output_alignments_chr=$dirprochr."/alignments";
	my $output_alignments_plas=$dirproplas."/alignments";
	$output_plasmid=~s/_contigslist.txt/.fasta/g;
	$output_chr=~s/_contigslist.txt/.fasta/g;
	my $name_output_plasmid=$output_plasmid;
 	$name_output_plasmid=~s/.fasta$//;
	my $name_output_chr=$output_chr;
 	$name_output_chr=~s/.fasta$//;
	return ($dirproplas, $dirprochr,$input_plasmid,$input_chr,$output_plasmid,$output_chr,$name_output_plasmid,$name_output_chr,$output_alignments_chr,$output_alignments_plas);
}
 
#only needed when plasmid option selected extracts the corresponding contigs
 sub sub6_samtools{
	my $genome = shift;
	my $seq=$dirModifiedGenomes."/".$genome;
    my ($dirproplas, $dirprochr,$input_plasmid,$input_chr,$output_plasmid,$output_chr,$name_output_plasmid,$name_output_chr) = var_inputs($genome);
	if( -s $input_chr){
			my $temp_chr=$output_chr.".tmp";
			my $temp_chr_1=$output_chr.".tmp1";
			my $cmd_cut="cut -d \" \" -f1 $input_chr >$temp_chr";
			my $cmd_awk="xargs samtools faidx $seq< $temp_chr > $temp_chr_1";
			my $cmd_awk1="awk '/^>/{gsub(\" \",\"_\"); \$0=\$0\":$genome"."_chr\"}1' $temp_chr_1 >$output_chr";
			#print "\n\n\n$cmd_cut\n\n\n";
			system($cmd_cut);
			system($cmd_awk);
			system($cmd_awk1);
			system( "rm $temp_chr");
			system( "rm $temp_chr_1");
	}
	if($plasmid_pred eq "SI"){
		if(-s $input_plasmid){
			my $temp_plasmid=$output_plasmid.".tmp";
			my $temp_plasmid_1=$output_plasmid.".tmp1";
			my $cmd_cut="cut -d \" \" -f1 $input_plasmid >$temp_plasmid";
			my $cmd_awk="xargs samtools faidx $seq< $temp_plasmid > $temp_plasmid_1";
			my $cmd_awk1="awk '/^>/{gsub(\" \",\"_\"); \$0=\$0\":$genome"."_plasmid\"}1' $temp_plasmid_1 >$output_plasmid";
			#print "\n\n\n$cmd_cut\n\n\n";
			system($cmd_cut);
			system($cmd_awk);
			system($cmd_awk1);
			system("rm $temp_plasmid");
			system("rm $temp_plasmid_1");
		}
	}
		
		
}
 
sub sub7_prodigal{
	my $seq = shift;
	my $genome=basename($seq);
    my ($dirproplas, $dirprochr,$input_plasmid,$input_chr,$output_plasmid,$output_chr,$name_output_plasmid,$name_output_chr,$output_alignments_chr,$output_alignments_plas) = var_inputs($genome);
	if(-s $output_chr ){
		my $genes_gff_chr=$dirprochr."/".$genome.".chr.genes.gff";
		my $genes_gff_chrtmp=$dirprochr."/".$genome.".chr.genes.gfftmp";
		my $genes_faa_chr=$dirprochr."/".$genome.".chr.genes.faa";
		my $genes_faa_chrtmp=$dirprochr."/".$genome.".chr.genes.faatmp";
		my $genes_fasta_chr=$dirprochr."/".$genome.".chr.genes.fasta";
		my $genes_fasta_chrtmp=$dirprochr."/".$genome.".chr.genes.fastatmp";

		system("prodigal -i $output_chr -o $genes_gff_chrtmp -a $genes_faa_chr -d $genes_fasta_chr -f gff -q");
		system("\(cat $genes_gff_chrtmp; echo '##FASTA'; cat  $output_chr\) >$genes_gff_chr");              
		my $cmd_awk1="awk '/^>/{\$0=\$0\":$genome\"}1'  $genes_gff_chr > $genes_gff_chrtmp";
		my $cmd_awk2="awk '/^>/{\$0=\$0\":$genome\"}1'  $genes_faa_chr > $genes_faa_chrtmp";
		my $cmd_awk3="awk '/^>/{\$0=\$0\":$genome\"}1'  $genes_fasta_chr > $genes_fasta_chrtmp";
		system($cmd_awk1);
		system($cmd_awk2);
		system($cmd_awk3);
		system("mv $genes_gff_chrtmp $genes_gff_chr");
		system("mv $genes_faa_chrtmp $genes_faa_chr");
		system("mv $genes_fasta_chrtmp $genes_fasta_chr");
	}
	if($plasmid_pred eq "SI"){
		if(-s $output_plasmid){
			my $genes_gff_plas= $dirproplas."/".$genome.".plasmid.genes.gff";
			my $genes_gfftmp_plas= $dirproplas."/".$genome.".plasmid.genes.gfftmp";
			my $genes_faa_plas= $dirproplas."/".$genome.".plasmid.genes.faa";
			my $genes_faatmp_plas= $dirproplas."/".$genome.".plasmid.genes.faatmp";
			my $genes_fasta_plas= $dirproplas."/".$genome.".plasmid.genes.fasta";
			my $genes_fastatmp_plas= $dirproplas."/".$genome.".plasmid.genes.fastatmp";
	
			system("prodigal -i $output_plasmid -o $genes_gfftmp_plas -a $genes_faa_plas -d $genes_fasta_plas -f gff -q");
			system("\(cat $genes_gfftmp_plas; echo '##FASTA'; cat $output_plasmid\) >$genes_gff_plas");
			my $cmd_awk1="awk '/^>/{\$0=\$0\":$genome\"}1'  $genes_gff_plas > $genes_gfftmp_plas ";
			my $cmd_awk2="awk '/^>/{\$0=\$0\":$genome\"}1'  $genes_faa_plas > $genes_faatmp_plas ";
			my $cmd_awk3="awk '/^>/{\$0=\$0\":$genome\"}1'  $genes_fasta_plas > $genes_fastatmp_plas ";
			system($cmd_awk1);
			system($cmd_awk2);
			system($cmd_awk3);
			system("mv $genes_gfftmp_plas $genes_gff_plas ");
			system("mv $genes_faatmp_plas $genes_faa_plas ");
			system("mv $genes_fastatmp_plas $genes_fasta_plas");
		}
	}
}

sub modify_file{
	my $file=shift;
	my $tmp=$file.".tmp";
	my $cmd1="sed 's/ # /###/g' $file >$tmp";
	my $cmd2="awk '/^>/ {printf(\"\\n\%s\\n\",\$0);next; } { printf(\"\%s\",\$0);}  END {printf(\"\\n\");}' < $tmp > $file";	
	system($cmd1);
	system($cmd2);
	system("rm $tmp");
}

sub sub8_modifed_header_chrplas{
	my $seq = shift;
	my $genome=basename($seq);
    my ($dirproplas, $dirprochr,$input_plasmid,$input_chr,$output_plasmid,$output_chr,$name_output_plasmid,$name_output_chr,$output_alignments_chr,$output_alignments_plas) = var_inputs($genome);
	my $genes_gff_chr=$dirprochr."/".$genome.".chr.genes.gff";
	my $genes_faa_chr=$dirprochr."/".$genome.".chr.genes.faa";
	my $genes_fasta_chr=$dirprochr."/".$genome.".chr.genes.fasta";
	if(-s $genes_fasta_chr){
		modify_file($genes_faa_chr);
		modify_file($genes_fasta_chr);
	}
	if($plasmid_pred eq "SI"){
		my $genes_gff_plas= $dirproplas."/".$genome.".plasmid.genes.gff";
		my $genes_faa_plas= $dirproplas."/".$genome.".plasmid.genes.faa";	
		my $genes_fasta_plas= $dirproplas."/".$genome.".plasmid.genes.fasta";
		if(-s $genes_fasta_plas){
			modify_file($genes_faa_plas);
			modify_file($genes_fasta_plas);
		}
	}
}

sub sub9_blast{
	my $seq = shift;
	my $genome=basename($seq);
    my ($dirproplas, $dirprochr,$input_plasmid,$input_chr,$output_plasmid,$output_chr,$name_output_plasmid,$name_output_chr,$output_alignments_chr,$output_alignments_plas) = var_inputs($genome);
	
	my $genes_gff_plas= $dirproplas."/".$genome.".plasmid.genes.gff";
	my $genes_faa_plas= $dirproplas."/".$genome.".plasmid.genes.faa";
	my $genes_fasta_plas= $dirproplas."/".$genome.".plasmid.genes.fasta";
	my $temp=$genes_fasta_plas.".tmp";
	
	my $genes_gff_chr=$dirprochr."/".$genome.".chr.genes.gff";
	my $genes_faa_chr=$dirprochr."/".$genome.".chr.genes.faa";
	my $genes_fasta_chr=$dirprochr."/".$genome.".chr.genes.fasta";
	my $temp1=$genes_fasta_chr.".tmp";

	foreach my $gene(@genes){
		my $genename=basename($gene);
		my $gene_id=$dirModifiedGenes."/".$genename;
		$gene_id=~s/.fasta$//g;
		my $gene_tag=$genename;
		$gene_tag=~s/.fasta$//g;
			
		if ($genes_type eq "nucl"){	
			if(-s $genes_fasta_chr){
				my $outputblast_chr= $output_alignments_chr."/".$gene_tag.".Blast.txt";
				my $cmd_blast_chr="blastn -query $genes_fasta_chr -db $gene_id -dust no -outfmt \"6 qseqid sseqid pident qlen slen length qstart qend sstart send mismatch gapope evalue bitscore\" -perc_identity $percentage -out $outputblast_chr";
				system ($cmd_blast_chr);
			}
			if( $plasmid_pred eq "SI"){
				if(-s $genes_fasta_plas){
					my $outputblast_plas= $output_alignments_plas."/".$gene_tag.".Blast.txt";
					my $cmd_blast_plas="blastn -query $genes_fasta_plas -db $gene_id -dust no -outfmt \"6 qseqid sseqid pident qlen slen length qstart qend sstart send mismatch gapope evalue bitscore\" -perc_identity $percentage -out $outputblast_plas";
					system ($cmd_blast_plas);
					}
				}	
		}elsif ($genes_type eq "prot"){
			if(-s $genes_fasta_chr){
				my $outputblast_chr= $output_alignments_chr."/".$gene_tag.".Blast.txt";
				my $tempblast_chr=$outputblast_chr.".tmp";
				my $cmd_blast="blastp -query $genes_faa_chr -db $gene_id -outfmt \"6 qseqid sseqid pident qlen slen length qstart qend sstart send evalue bitscore qcovs qcovhsp scovhsp\"  -out $outputblast_chr";
				system ($cmd_blast);
				my $cmd_highper="awk '\$3>=".$percentage." && \$13>=".$percentage." {print \$0}' < $outputblast_chr >$tempblast_chr";
				system ($cmd_highper);
				system ("mv $tempblast_chr $outputblast_chr")
			}
			if( $plasmid_pred eq "SI"){
				if(-s $genes_fasta_plas){					
					my $outputblast_plas= $output_alignments_plas."/".$gene_tag.".Blast.txt";
					my $tempblast_plas=$outputblast_plas.".tmp";
					my $cmd_blast1="blastp -query $genes_faa_plas -db $gene_id -outfmt \"6 qseqid sseqid pident qlen slen length qstart qend sstart send evalue bitscore qcovs qcovhsp scovhsp\"  -out $outputblast_plas";
					system ($cmd_blast1);
					my $cmd_highper="awk '\$3>=".$percentage." && \$13>=".$percentage." {print \$0}' < $outputblast_plas >$tempblast_plas";
					system ($cmd_highper);
					system ("mv $tempblast_plas $outputblast_plas");
					}
				}
			}
	}
}

sub sub10_curate_alignment{
	my $seq = shift;
	my $genome=basename($seq);
    my ($dirproplas, $dirprochr,$input_plasmid,$input_chr,$output_plasmid,$output_chr,$name_output_plasmid,$name_output_chr,$output_alignments_chr,$output_alignments_plas) = var_inputs($genome);

	my $genes_gff_plas= $dirproplas."/".$genome.".plasmid.genes.gff";
	my $genes_faa_plas= $dirproplas."/".$genome.".plasmid.genes.faa";
	my $genes_fasta_plas= $dirproplas."/".$genome.".plasmid.genes.fasta";
	my $temp=$genes_fasta_plas.".tmp";

	my $genes_gff_chr=$dirprochr."/".$genome.".chr.genes.gff";
	my $genes_faa_chr=$dirprochr."/".$genome.".chr.genes.faa";
	my $genes_fasta_chr=$dirprochr."/".$genome.".chr.genes.fasta";
	my $temp1=$genes_fasta_chr.".tmp";
	
	foreach my $gene(@genes){
		my $genename=basename($gene);
		my $gene_id=$dirModifiedGenes."/".$genename;
		$gene_id=~s/.fasta$//g;
		my $gene_tag=$genename;
		$gene_tag=~s/.fasta$//g;

		my $gene_ref=$gene_id.".fasta";
		
		my $outputblast_chr= $output_alignments_chr."/".$gene_tag.".Blast.txt";
		my $tempblast_chr=$outputblast_chr.".tmp";
		
		if(-s $outputblast_chr){
			my $temp_chr=$output_alignments_chr."/".$gene_tag.".tmp";
			system("cut -f 1 $outputblast_chr >$temp_chr");
			my $output_gene_chr= $output_alignments_chr."/".$gene_tag.".fasta";
			if( $genes_type eq 'nucl'){
				system("xargs samtools faidx $genes_fasta_chr < $temp_chr >$output_gene_chr");
			}elsif($genes_type eq 'prot'){
				system("xargs samtools faidx $genes_faa_chr < $temp_chr >$output_gene_chr");
			}
				my $cmd_awk="awk '/^>/ {printf(\"\\n\%s\\n\",\$0);next; } { printf(\"\%s\",\$0);}  END {printf(\"\\n\");}' < $output_gene_chr > $output_gene_chr.tmp";	
				system($cmd_awk);
				system("mv $output_gene_chr.tmp $output_gene_chr");
				system("rm $temp_chr");
				my $output_ref_chr=$output_gene_chr.".plusRef.fasta";
				my $cmd="\(cat $gene_ref; cat $output_gene_chr\) > $output_ref_chr";
				system($cmd);
				my $out_alignment_description=$dirresults."/".$genome.".".$gene_tag.".chr.out.alignments.description";
				open OUT, ">$out_alignment_description" or die "Cannot open $out_alignment_description for writing\n";
				my $gene_genome_chr_out_alignment_AA_fasta=$dirresultspep."/".$gene_tag."_".$genome.".chr.AAsequence.fasta";
				open OUT1, ">$gene_genome_chr_out_alignment_AA_fasta" or die "Cannot open $gene_genome_chr_out_alignment_AA_fasta for writing\n";
				my $gene_genome_chr_out_alignment_nt_fasta=$dirresults."/".$gene_tag."_".$genome.".chr_alignment.fasta.tmp";
				open OUT2, ">$gene_genome_chr_out_alignment_nt_fasta" or die "Cannot open  $gene_genome_chr_out_alignment_nt_fasta for writing\n";
				my $rcmd=" Rscript $results_folder/scripts/curatingAlignments.r  $output_ref_chr $gene_tag $gene_genome_chr_out_alignment_nt_fasta $gene_genome_chr_out_alignment_AA_fasta $out_alignment_description --slave";
				print "-------$rcmd-------\n\n";
				system($rcmd);
		}else { 
 				my $out_alignment_description=$dirresults."/".$genome.".".$gene_tag.".chr.out.alignments.description";
				open OUT, ">$out_alignment_description" or die "Cannot open $out_alignment_description for writing\n";
				my $firstline_Ref=`grep ">" $gene_ref`;
				chomp($firstline_Ref);
				$firstline_Ref=~s/[\n\r]//g;
				$firstline_Ref=~s/>//g;
				print OUT "$genename\t\"$firstline_Ref\"\tNA\tNA\tNA\t\"$genome\"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tchromosome\tabsent\t0\t0\t0\n";
				close OUT;
	 	}
		my $outputblast_plas= $output_alignments_plas."/".$gene_tag.".Blast.txt";
		my $tempblast_plas=$outputblast_plas.".tmp";
		if($plasmid_pred eq "SI"){
			if(-s $outputblast_plas){
				
				my $temp_plasmid=$output_alignments_plas."/".$gene_tag.".tmp";
				system("cut -f 1 $outputblast_plas >$temp_plasmid");
				my $output_gene_plas= $output_alignments_plas."/".$gene_tag.".fasta";
				if( $genes_type eq 'nucl'){
					system("xargs samtools faidx $genes_fasta_plas< $temp_plasmid >$output_gene_plas");
				}elsif($genes_type eq 'prot'){
						system("xargs samtools faidx $genes_faa_plas< $temp_plasmid >$output_gene_plas");
				}
				my $cmd_awk="awk '/^>/ {printf(\"\\n\%s\\n\",\$0);next; } { printf(\"\%s\",\$0);}  END {printf(\"\\n\");}' < $output_gene_plas > $output_gene_plas.tmp";	
				system($cmd_awk);
				system("mv $output_gene_plas.tmp $output_gene_plas");
				system("rm $temp_plasmid");
				my $output_ref_plasmid=$output_gene_plas.".plusRef.fasta";
				my $cmd="\(cat $gene_ref; cat $output_gene_plas\) > $output_ref_plasmid";
				system($cmd);
				my $out_alignment_description=$dirresults."/".$genome.".".$gene_tag.".plasmid.out.alignments.description";
				open OUT, ">$out_alignment_description" or die "Cannot open $out_alignment_description for writing\n";
				my $gene_genome_plas_out_alignment_AA_fasta=$dirresultspep."/".$gene_tag."_".$genome.".plasm.AAsequence.fasta";
				open OUT1, ">$gene_genome_plas_out_alignment_AA_fasta" or die "Cannot open $gene_genome_plas_out_alignment_AA_fasta for writing\n";
				my $gene_genome_plas_out_alignment_nt_fasta=$dirresults."/".$gene_tag."_".$genome.".plasm_alignment.fasta.tmp";
				open OUT2, ">$gene_genome_plas_out_alignment_nt_fasta" or die "Cannot open  $gene_genome_plas_out_alignment_nt_fasta for writing\n";
				my $rcmd="Rscript $results_folder/scripts/curatingAlignments.r $output_ref_plasmid $gene_tag $gene_genome_plas_out_alignment_nt_fasta $gene_genome_plas_out_alignment_AA_fasta $out_alignment_description --slave";
				print "--------- $rcmd -------------\n\n";
				system($rcmd);
			}else { 
				my $out_alignment_description=$dirresults."/".$genome.".".$gene_tag.".plasmid.out.alignments.description";
				open OUT, ">$out_alignment_description" or die "Cannot open $out_alignment_description for writing\n";
				my $firstline_Ref=`grep ">" $gene_ref`;
				chomp($firstline_Ref);
				$firstline_Ref=~s/[\n\r]//g;
				$firstline_Ref=~s/>//g;
				print OUT "$genename\t\"$firstline_Ref\"\tNA\tNA\tNA\t\"$genome\"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tplasmid\tabsent\t0\t0\t0\n";
				close OUT;
				}	
		}
	}
}


####################################################################################################################################################################################################################################################################################################################################################################
#parallel

sub afork (\@$&) {
	my ($data, $max, $code) = @_;
	my $c = 0;
	my %hash;
	foreach my $data (@$data) {
		wait unless ++ $c <= $max;
		die "Fork failed: $!\n" unless defined (my $pid = fork);
		exit $code -> ($data) unless $pid;
		print($data);
		$hash{$pid}=$data;
	}
	1 until -1 == wait;
	#my $child_pid = wait();
	#my $status = $?;
	#if ( $status == -1  ) { die "wait failed File: ".$hash{$child_pid}."\n"; }
	#elsif ( $status & 127) { warn "Child $child_pid killed by signal File: ".$hash{$child_pid}."\n".( $? & 127)."\n"; }
	#elsif ( $status & 128) { warn "Child $child_pid. There was a core dump! File: ".$hash{$child_pid}."\n ".(  $? & 128)."\n"; }
	#elsif ( $status & 0x7F ) { warn "Child $child_pid killed by signal File: ".$hash{$child_pid}."\n".( $? & 0x7F )."\n"; }
	#elsif ( $status >> 8   ) { warn "Child $child_pid exited with error File: ".$hash{$child_pid}."\n".( $? >> 8 )."\n"; }
	#else                     { print "Child $child_pid exited successfully\n"; }

}

#######################################################################################################################################################################################################################################################################################################################################################################
#CREATING R SCRIPTS

# this bit can be use to use mlplasmid or plasclass method
 if($plasmid_pred eq "SI"){
 $R_file="$results_folder/scripts/plasmid_chromosome_Detection.r";
	
		open R_SCRIPT,">$R_file" or die "Cannot write $R_file script\n";

        	 $R_script= "
    rm(list=ls()); 
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
Sys.setenv('R_MAX_VSIZE'=32000000000)
suppressMessages(suppressWarnings(library(msa)))
suppressMessages(suppressWarnings(library(reshape2)))
suppressMessages(suppressWarnings(library(seqinr)))
suppressMessages(suppressWarnings(library( Biostrings)))
suppressMessages(suppressWarnings(library(stringr)))
options(warn=-1)
options(stringsAsFactors = FALSE)

args = commandArgs(trailingOnly=TRUE)
file=args[1]
genename=args[2]
out1=args[3]
out2=args[4]
out3=args[5]
setwd(dirname(file))
#print(getwd())



##Functions

findPotentialStartsAndStops <- function(DNA_string){ 
	positions =c()
   	types =c()
    if(length(positions)==0){
    	codons <- c('atg', 'taa', 'tag', 'tga','ATG', 'TAA', 'TAG', 'TGA')
     	for (i in  1:length(codons))
     	{
        	codon <- codons[i]
        	occurrences <- matchPattern(codon, DNA_string)
      		codonpositions <- occurrences\@ranges\@start  
        	numoccurrences <- length(codonpositions) 	
   			positions <- append(positions,codonpositions, after=length(positions))
     		types <- append(types,rep(codon, numoccurrences), after=length(types))
 		}
    indices <- order(positions)
    positions <- positions[indices]
    types <- types[indices]
    mylist <- list(positions,types)
 	return(mylist)
	}
}     

findORFsinSeq <- function(sequence){
	mylist <- findPotentialStartsAndStops(sequence)
    positions <- mylist[[1]]
    types <- mylist[[2]]
    orfstarts <- numeric()
    orfstops <- numeric()
    orflengths <- numeric()
    numpositions <- length(positions)
 
    if (numpositions >= 2){
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
  

seq2Fasta <- function(seq,name,filename){
  						sink(filename,append=TRUE)
    					cat(paste0('>', name),file=filename,append=TRUE)
    					cat('\\n',file=filename,append=TRUE)
    					the.sequence <- toString(seq)
    					the.sequence=gsub('N','-',the.sequence)
    					cat(the.sequence,file=filename,append=TRUE)
    					cat('\\n',file=filename,append=TRUE)  
  						sink(NULL)
						}

seqAA2Fasta <- function(seq,name,filename){
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



### voy aqui
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
Sys.setenv('R_MAX_VSIZE'=32000000000)
suppressMessages(suppressWarnings(library(msa)))
suppressMessages(suppressWarnings(library(reshape2)))
suppressMessages(suppressWarnings(library(seqinr)))
suppressMessages(suppressWarnings(library( Biostrings)))
suppressMessages(suppressWarnings(library(stringr)))
options(warn=-1)
options(stringsAsFactors = FALSE)
args = commandArgs(trailingOnly=TRUE)
file=args[1]
genename=args[2]
out1=args[3]
out2=args[4]
out3=args[5]
setwd(dirname(file))
#print(getwd())
				
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



sub sub11_cleanup_compress{
my $gene= shift;
	my $genename=basename($gene);
	$genename=~s/.fasta$//g;
	my $outAlign=$dirresults."/".$genename."_alignment.fasta";
	my $outAlign_forR=$dirresults."/".$genename."_alignment.fasta_forR";
	my $outAlign_forRR=$dirresults."/".$genename."_alignment.fasta_forR_forR";
	my $cmd_align="cat ".$dirresults."/".$genename."*_alignment.fasta.tmp >$outAlign";
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
suppressMessages(suppressWarnings(library(stringr)))
options(warn=-1)
				
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
	

######################################################>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PROGRAM RUNS HERE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<##########################################################

#######################################################################################################################################################################################################################################################################################################################################################################
#check dependencies

my @programs = ("cat", "mkdir", "sed","chmod","Rscript","rm","find","xargs","awk","grep","mv","makeblastdb","cp","cut","samtools","prodigal");
foreach my $program(@programs){print "checking for $program\n"; check_exists_command $program or die "$program doesn't exists"}

#######################################################################################################################################################################################################################################################################################################################################################################
#copy input files 
foreach my $files(@seq){ if (-e $files) { copy($files,$dirModifiedGenomes);} else { die "$files does not exist! Check yout command line\n";}};
foreach my $files(@genes){ if (-e $files){ copy($files,$dirModifiedGenes);}  else { die "$files does not exist! Check yout command line\n";}};

opendir(my $dir_open_seq, $dirModifiedGenomes) or die "Cannot open directory $dirModifiedGenomes: $!";
my @seq_array = grep { -T "$dirModifiedGenomes/$_" } readdir $dir_open_seq;
#@seq = readdir $dir_open_seq;
closedir $dir_open_seq;

opendir(my $dir_open_genes, $dirModifiedGenes) or die "Cannot open directory $dirModifiedGenes: $!";
my @genes_array = grep { -T "$dirModifiedGenes/$_" } readdir $dir_open_genes;
#@genes = readdir $dir_open_genes;
closedir $dir_open_genes;



#######################################################################################################################################################################################################################################################################################################################################################################


#######################################################################################################################################################################################################################################################################################################################################################################
#create tmp output file
my $out_alignment_description=$dirresults."/out.alignments.description.txt";
open OUT, ">$out_alignment_description" or die "Cannot open $out_alignment_description for writing\n";
print OUT "Gene\tRef\tRef_length_nt\tRef_length_AA\tRef_stopCodon\tGenomeFile\tQuery\tStart\tEnd\tStrand\tQuery_length_nt\tQuery_length_AA\tnt_gaps\tstopCodon\tSNPs\tGAPs\tINSERTION\tlocation\tclassification\tPercentageSimilarityAA\tPercentageSimilarityDNA\tN.copies\n";
close OUT;
#######################################################################################################################################################################################################################################################################################################################################################################

print "
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

\t\t READING AND MODIFYING GENES

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
  
afork(@genes_array,$cores,\&sub1_genemodification);
#=debug

print "
-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-

\t\t DONE READING AND MODIFYING GENES
  
-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-\n\n";
####################################################################################################################################################################################################################
print "
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

\t\t CREATING DATABASES FROM GENES OF INTEREST

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
  
afork(@genes_array,$cores,\&sub2_gene_makeblastdb);

print "
-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-

\t\t DONE CREATING DATABASES FROM GENES OF INTEREST
  
-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-\n\n";
####################################################################################################################################################################################################################
print "
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

\t\t CREATING IMPUT FOLDERS WITH MODIFIED INPUTS

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
 
afork(@seq_array,$cores,\&sub3_chrplas_folders);

print "
-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-

\t\t DONE CREATING IMPUT FOLDERS WITH MODIFIED INPUTS
  
-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-\n\n";

####################################################################################################################################################################################################################
 
   print "
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

\t\t PREDICTING PLASMID AND CHROMOSOMES FOR ALL GENOMES

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";

 if($plasmid_pred eq "SI"){
	 afork(@seq_array,$cores,\&sub4_plasclass_method);
	 afork(@seq_array,$cores,\&sub5_plasmid_extraction);
	 afork(@seq_array,$cores,\&sub6_samtools);
	  print "
-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-

\t\t DONE PREDICTING PLASMID AND CHROMOSOMES FOR ALL GENOMES\n
  
-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-\n\n";
 
 
 	}else{ 
  print "
-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-

\t\t NO PREDICTION OF PLASMIDS WAS PERFORMED\n
set -p option if you want to predict plasmids in your assembly
  
-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-\n\n";
 
 }
####################################################################################################################################################################################################################
   print "
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

\t\t PREDICTING GENES FROM GENOMES 

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";

 afork(@seq_array,$cores,\&sub7_prodigal);
  
print "
-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-

\t\t DONE PREDICTING GENES FROM GENOMES 
  
-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-\n\n";
 
####################################################################################################################################################################################################################
   print "
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

\t\t FORMATING

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";

 afork(@seq_array,$cores,\&sub8_modifed_header_chrplas);
  
print "
-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-

\t\t DONE  FORMATING
  
-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-\n\n";



####################################################################################################################################################################################################################
   print "
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

\t\t BLASTING

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
 
 afork(@seq_array,$cores,\&sub9_blast);
  
print "
-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-

\t\t DONE BLASTING
  
-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-\n\n";

####################################################################################################################################################################################################################
   print "
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

\t\t PREDICTING GENES FROM GENOMES 

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";

 afork(@seq_array,$cores,\&sub10_curate_alignment);
  
print "
-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-

\t\t DONE PREDICTING GENES FROM GENOMES 
  
-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-\n\n";
#####################################>>>>>>>>>>>>>>>>>> CLEAN UP AND COMPRESS <<<<<<<<<<<<<<<<<<<<<<################################################
#normal(@seq);
	
	
	  print "
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

\t\t CONCATENATION OF ALIGMENTS 

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";

 afork(@genes_array,$cores,\&sub11_cleanup_compress);

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

#system($cmd_fai);
system("rm $diroutput/*/*fai");

 
 
 my $duration = time - $start;
print "Execution time: $duration s\n";

	

 #=
 
 



























