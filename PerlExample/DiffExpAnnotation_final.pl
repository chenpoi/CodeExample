#!/usr/bin/perl
use strict;
use warnings;
use feature 'say';
use List::MoreUtils qw(uniq);
use BLAST;
use GO;
use DiffExp;
use Report;

my $path = "PATH/TO/DATA";
main();

sub main {
	
	# Hash reference to store gene-to-GO term mappings.
	my $gene_to_GO_ref;

	# Hash reference to store GO term-to-GO description mappings.
	my $GO_to_description_ref;

	# Hash reference to store BLAST mappings.
	my $transcript_to_protein_ref;

	$transcript_to_protein_ref = read_blast();
	$GO_to_description_ref = read_GO_desc();
	$gene_to_GO_ref = read_gene_to_GO();
	
	#print a report
	print_report($GO_to_description_ref, $transcript_to_protein_ref, $gene_to_GO_ref);
	
}

#this sub will store transcript ID and protein swissprot ID in a hash
sub read_blast {
	my $blast_file      = $path.'/blastp.outfmt6';
	open( my $blast_fh, '<', $blast_file ) or die $!;
		
	my %blastHash;
	while ( <$blast_fh> ) {
		chomp;
		my $blastTerm = BLAST->new(record => $_);
		$blastHash{$blastTerm->get_transcript()} = $blastTerm;
	}
	return \%blastHash;
}
#this sub will store protein swissprot ID and corresponding GO ID in a hash
sub read_gene_to_GO {
	my $gene_to_GO_file = $path.'gene_association_subset.txt';
	open( my $GENE_TO_GO_fn, '<', $gene_to_GO_file ) or die $!;
	my %gene_to_GO;
	#read gene_to_go file line by line
	#notice that there are some proteins belonging to more than one GO terms
	#I only choose the first one as its corresponding GO ID
	while ( my $line = <$GENE_TO_GO_fn>) {
		chomp $line;
		my @entity = split /\t/, $line;
		# this if statement will check if that protein has been recorded before
		# if yes, the related GO term will add to its reference array
		# if not, new reference will created
		if (defined $gene_to_GO{$entity[1]}) {
			push $gene_to_GO{$entity[1]}, $entity[4];
		} else {
			$gene_to_GO{$entity[1]} = [$entity[4]];
		}
	}
	return \%gene_to_GO;
}
#this sub will store GO ID and name in a hash
sub read_GO_desc {
	my $GO_desc_file    = $path.'go-basic.obo';
	open( my $GO_DESC_fh, '<', $GO_desc_file ) or die $!;
	my %GOHash;
	#sepqrate character is set as "[Term]" to split the file into different go terms
	local $/ = '[Term]';
	#read GO terms one by one
	#regular expression is applied to locate and extract GO ID and name
	while ( my $line = <$GO_DESC_fh> ) {
		chomp $line;
		my $GO_obj = GO->new(record => $line);
		if ($GO_obj->get_id()) {
			$GOHash{$GO_obj->get_id()} = $GO_obj;
		}
	}
	#store them in a hash
	return \%GOHash;
}
#this sub will firstly process differential expression file and then combine the results from four files into a new output file
sub print_report {
	# pass arguments
	my @args = @_;
	my %GO_to_description = %{$args[0]};
	my %transcript_to_protein = %{$args[1]};
	my %gene_to_GO = %{$args[2]};

	my $diff_exp_file   = $path.'diffExpr.P1e-3_C2.matrix';
	open( my $DIFF_EXP_fh, '<', $diff_exp_file ) or die $!;
	
	#initialze some variables
	my $beginning = 1;
	my %stats;
	my @transcripIDList;
	#transcript ID and analysis results are stored in two arrays in order
	while (my $line = <$DIFF_EXP_fh>) {
		chomp $line;
		if ($beginning) {
			$beginning = 0;
		} else {
			my $diff_obj = DiffExp->new(record => $line);
			$stats{$diff_obj->get_transcript()} = $diff_obj;
			push @transcripIDList, $diff_obj->get_transcript();
		}
	}
	#because there are some missing values, I use defined here to check them
	foreach my $i (@transcripIDList){	
		if (defined $transcript_to_protein{$i}) {
				my $proteinID = $transcript_to_protein{$i}->get_sp();
				my $GOIDList = $gene_to_GO{$proteinID};
				if (defined $GOIDList) {
					my @GOTermList = @GO_to_description{uniq @$GOIDList};
					my $report_obj = Report->new({
						diff_expressions	=> $stats{$i},
						protein_id			=> $proteinID,
						protein_desc		=> get_protein_info_from_blast_DB($proteinID),	
						GO_terms			=> \@GOTermList,	
					});
					$report_obj->print_all();	
				}
		}	
	}	
}
sub get_protein_info_from_blast_DB {
	my ( $protein_id ) = @_;
	my $db = '/blastDB/swissprot';
	my $exec = 'blastdbcmd -db '
	  . $db
	  . ' -entry '
	  . $protein_id
	  . ' -outfmt "%t" -target_only | ';

	open( SYSCALL, $exec ) or die "Can't open the SYSCALL ", $!;

	# Set default description as "NA". If protein is found in DB, overwrite 
	# description.
	my $protein_description = 'NA';
	while (<SYSCALL>) {
		chomp;
		if ( $_ =~ /RecName:\s+(.*)/i ) {
			$protein_description = $1;
		}
	}
	close SYSCALL;
	return $protein_description;
}


