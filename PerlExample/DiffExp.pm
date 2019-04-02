package DiffExp;

use Moose;
use MooseX::FollowPBP;
use feature qw(say);

sub BUILD {
	my ($self,$args) = @_;
	my $record = $args->{record};
	my @words = split /\s+/, $record;
	
	$self->{transcript} = $words[0];
	$self->{sp_ds} = $words[1];
	$self->{sp_hs} = $words[2];
	$self->{sp_log} = $words[3];
	$self->{sp_plat} = $words[4];
	

}

has 'transcript' => (
	is 			=> 'ro',
	isa			=> 'Str',
	clearer		=> 'clear_transcript',
	predicate	=> 'predicate_transcript',
);


has 'sp_ds' => (
	is 			=> 'ro',
	isa			=> 'Num',
	clearer		=> 'clear_sp_ds',
	predicate	=> 'predicate_sp_ds',
);


has 'sp_hs' => (
	is 			=> 'ro',
	isa			=> 'Num',
	clearer		=> 'clear_sp_hs',
	predicate	=> 'predicate_sp_hs',
);



has 'sp_log' => (
	is 			=> 'ro',
	isa			=> 'Num',
	clearer		=> 'clear_sp_log',
	predicate	=> 'predicate_sp_log',
);


has 'sp_plat' => (
	is 			=> 'ro',
	isa			=> 'Num',
	clearer		=> 'clear_sp_plat',
	predicate	=> 'predicate_sp_plat',
);

1;
