package GO;

use Moose;
use MooseX::FollowPBP;
use feature qw(say);


has 'id' => (
	is 	=> 'ro',
	isa	=> 'Str',
	clearer		=> 'clear_id',
	predicate	=> 'predicate_id',
);

has 'name' => (
	is	=> 'ro',
	isa	=> 'Str',
	clearer		=> 'clear_name',
	predicate	=> 'predicate_name',
);

has 'namespace' => (
	is	=> 'ro',
	isa	=> 'Str',
	clearer		=> 'clear_namespace',
	predicate	=> 'predicate_namespace',
);

has 'definition' => (
	is	=> 'ro',
	isa	=> 'Str',
	clearer		=> 'clear_definition',
	predicate	=> 'predicate_definition',
);

has 'ISA' => (
	is	=> 'ro',
	isa	=> 'ArrayRef[Str]',
	default => sub {[]},
	lazy => 1,
	clearer		=> 'clear_ISA',
	predicate	=> 'predicate_ISA',
);

has 'alt_id' => (
	is	=> 'ro',
	isa	=> 'ArrayRef[Str]',
	default => sub {[]},
	lazy => 1,
	clearer		=> 'clear_alt_id',
	predicate	=> 'predicate_alt_id',
);


sub BUILD {
	# pass argument firstly
	my ($self, $args) = @_;
	my $record = $args->{record};
	# following three regulat expressions
	my $parsing_regex =
			qr/									
			^id:\s(?<GOID>GO:\d{7})\n			
			^name:\s(?<GOName>.*?)\n			
			^namespace:\s(?<nameSpace>[a-z]+_[a-z]+).*?\n	
			^def:\s"(?<def>.*?)"\s\[.*?\n			
			/msx;								
	my $findIsa = qr/
			^is_a:\s+(?<isa>.*?)\s+!
		/msx;
	my $findAltId = qr/
			^alt_id:\s+(?<alt_id>.*?)\s+
		/msx;
	if ( $record =~ /$parsing_regex/ ) {
			# accessor functions are used to set attrbute values
			$self->{id} = $+{GOID};
			$self->{name} = $+{GOName};
			$self->{namespace} = $+{nameSpace};
			$self->{definition} = $+{def};
			my @alt_ids = ();
			while ( $record =~ /$findAltId/g ) {
				push( @alt_ids, $+{alt_id} );
			}
			if ( @alt_ids ) {
				$self->{alt_id} = \@alt_ids;
			}
			my @isas = ();
			while ( $record =~ /$findIsa/g ) {
				push( @isas, $+{isa} );
			}
			$self->{ISA} = \@isas;
	}
}

sub printGO {
	my $self = shift;
	say join("\t",$self->get_id(),$self->get_name(),$self->get_namespace(),$self->get_definition());
	if (defined $self->get_ISA()) {
		foreach (@{$self->get_ISA()}) {
			say 'is_a:',"\t",$_;
		}
	}
	if (defined $self->get_alt_id()) {
		foreach (@{$self->get_alt_id()}) {
			say 'alt_id:',"\t",$_;
		}
	}
	say "";
}
1;
