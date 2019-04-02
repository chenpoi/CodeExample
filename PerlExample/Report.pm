package Report;

use Moose;
use MooseX::FollowPBP;
use feature qw(say);

open( REPORT, '>', 'final_report.txt' ) or die $!;

sub BUILD {
	my ($self, $args) = @_;
	$self->{diff_expressions} = $args->{diff_expressions};
	$self->{protein_id} = $args->{protein_id};
	$self->{protein_desc} = $args->{protein_desc};
	$self->{GO_terms} = $args->{GO_terms};	
}

sub print_all {
	my ($self) = @_;
	
	my $term_count = 0;
	
	foreach my $GO ( @{$self->get_GO_terms()} ) {
		$term_count++;	
		if (defined $GO) {
			if ( $term_count == 1 ) {
   			say REPORT join("\t",$self->get_diff_expressions()->get_transcript(), 
   								$self->get_diff_expressions()->get_sp_ds(),
   								$self->get_diff_expressions()->get_sp_hs(),
   								$self->get_diff_expressions()->get_sp_log(),
   								$self->get_diff_expressions()->get_sp_plat(),
   								$self->get_protein_id(),
   								$GO->get_id(),
   								$GO->get_name(),
   								$self->get_protein_desc(),
   								);
   			} else {
   				say REPORT "\t\t\t\t\t\t",$GO->get_id(),"\t",$GO->get_name();
   			}
		}
	}
}


has 'diff_expressions' => (
	is 			=> 'ro',
	isa			=> 'DiffExp',
	clearer		=> 'clear_diff_expressions',
	predicate	=> 'predicate_diff_expressions',
);

has 'protein_id' => (
	is 			=> 'ro',
	isa			=> 'Str',
	clearer		=> 'clear_protein_id',
	predicate	=> 'predicate_protein_id',
);

has 'protein_desc' => (
	is 			=> 'ro',
	isa			=> 'Str',
	clearer		=> 'clear_protein_desc',
	predicate	=> 'predicate_protein_desc',
);

has 'GO_terms' => (
	is 			=> 'ro',
	isa			=> 'ArrayRef',
	clearer		=> 'clear_GO_terms',
	predicate	=> 'predicate_GO_terms',
);


1;
