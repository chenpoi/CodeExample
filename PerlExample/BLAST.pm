package BLAST;

use Moose;
use MooseX::FollowPBP;
use feature qw(say);

has 'transcript' => (
	is 			=> 'ro',
	isa			=> 'Str',
	clearer		=> 'clear_transcript',
	predicate	=> 'predicate_transcript',
);

has 'isoform' => (
	is			=> 'ro',
	isa			=> 'Str',
	clearer		=> 'clear_isoform',
	predicate	=> 'predicate_isoform',
);

has 'gi' => (
	is			=> 'ro',
	isa			=> 'Int',
	clearer		=> 'clear_gi',
	predicate	=> 'predicate_gi',
);

has 'sp' => (
	is			=> 'ro',
	isa			=> 'Str',
	clearer		=> 'clear_sp',
	predicate	=> 'predicate_sp',
);

has 'prot' => (
	is			=> 'ro',
	isa			=> 'Str',
	clearer		=> 'clear_prot',
	predicate	=> 'predicate_prot',
);

has 'pident' => (
	is			=> 'ro',
	isa			=> 'Num',
	clearer		=> 'clear_pident',
	predicate	=> 'predicate_pident',
);

has 'len' => (
	is 			=> 'ro',
	isa			=> 'Int',
	clearer		=> 'clear_len',
	predicate	=> 'predicate_len',
);

has 'mismatch' => (
	is			=> 'ro',
	isa			=> 'Int',
	clearer		=> 'clear_mismatch',
	predicate	=> 'predicate_mismatch',
);

has 'gapopen' => (
	is			=> 'ro',
	isa			=> 'Int',
	clearer		=> 'clear_gapopen',
	predicate	=> 'predicate_gapopen',
);


sub BUILD {
	# pass argument
	my ($self, $args) = @_;
	my $record = $args->{record};
	# creat regular expression
	my $parsing_regex =
			qr/
			^(?<transcripID>.*?)\|
			(?<isoformID>.*)\t
			gi\|(?<giNumber>.*?)\|
			sp\|(?<swissProtID>.*?)\.\d\|
			(?<entryID>.*?)\t
			(?<identity>\d+)\t
			(?<Length>\d+)\t
			(?<mismatch>\d+)\t
			(?<gapOpen>\d+)\t
			/x;
	
	# get the value and pass them to accessor function
	if ($record =~ /$parsing_regex/) {
		$self->{transcript} = $+{transcripID};
		$self->{isoform} = $+{isoformID};
		$self->{gi} = $+{giNumber};
		$self->{sp} = $+{swissProtID};
		$self->{prot} = $+{entryID};
		$self->{pident} = $+{identity};
		$self->{len} = $+{Length};
		$self->{mismatch} = $+{mismatch};
		$self->{gapopen} = $+{gapOpen};
	}
}

sub print_all {
	my $self = shift;
	say join("\t",$self->get_transcript(), $self->get_isoform(),$self->get_gi(),
		$self->get_sp(),$self->get_prot(),$self->get_pident(),
		$self->get_len(),$self->get_mismatch(),$self->get_gapopen());
}



1;
