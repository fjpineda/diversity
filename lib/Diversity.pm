package Diversity;

use strict;
use warnings;

use File::Basename;
use Data::Dumper;
use Bio::SeqIO;
# use Benchmark;

# ----------------------------------------
# symbols & matrices
# ----------------------------------------

my $_null;			# null symbol
my $_gap;			# gap symbol
my @_residues=();		# residues
my @_alphabet=();		# alphabet
my @_nonnull_symbols=();# all symbols except null

# ----------------------------------------
# diversity and variance calculations
# ----------------------------------------

my %_m_mat=();	# mismatch matrix          $_m_mat{$symbol1}{$symbol2}
my %_c_mat=();	# pairwise coverage matrix $_c_mat{$symbol1}{$symbol2}

my $_K;           # number of reads in the alignment file
my $_W;           # width (number of columns) of the alignment file
my $_M;     	# expected number of mismatches
my $_Z;		# expected pairwise width
my $_D;     	# diversity
my $_varD;		# variance of the diversity
my $_sigmaD;	# standard deviation of D

my @_p=();		# symbol probabilities $_p[$i]{$symbol}
my @_m=();
my @_z=();
my @_var=();   	# variances of the symbol probabilities
my @_cov=();   	# covariance of the symbol probabilities

# ----------------------------------
# new -- initialize symbols and matrices
# ----------------------------------

my %_options=();

sub new 
{ 
	my $self = shift;
	%_options = @_;
	
	my %defaults = (
		INDELS => 0
	);
	foreach ('INDELS') {
		if(! defined($_options{$_})) {$_options{$_}= $defaults{$_}};
	}

	$_null     = 'N';
	$_gap      = '-';
	@_residues = ('A','T','C','G');
	@_alphabet = sort (@_residues, $_gap, $_null);

	if($_options{INDELS} ==0) {
		_do_subs();
	}
	else {
		_do_indels();
	}
	
	return bless{};
}


# ----------------------------------
# initialize indicators
# ----------------------------------
sub _do_subs {

	warn("do_subs\n");
	%_m_mat = (   # mismatch matrix
		'A'   =>{'A'=>0, 'T'=>1, 'C'=>1, 'G'=>1, $_gap=>0, $_null=>0},
		'T'   =>{'A'=>1, 'T'=>0, 'C'=>1, 'G'=>1, $_gap=>0, $_null=>0},
		'C'   =>{'A'=>1, 'T'=>1, 'C'=>0, 'G'=>1, $_gap=>0, $_null=>0},
		'G'   =>{'A'=>1, 'T'=>1, 'C'=>1, 'G'=>0, $_gap=>0, $_null=>0},
		$_gap =>{'A'=>0, 'T'=>0, 'C'=>0, 'G'=>0, $_gap=>0, $_null=>0},
		$_null=>{'A'=>0, 'T'=>0, 'C'=>0, 'G'=>0, $_gap=>0, $_null=>0}
	);
	
	%_c_mat = (   # mismatch matrix
		'A'   =>{'A'=>1, 'T'=>1, 'C'=>1, 'G'=>1, $_gap=>0, $_null=>0},
		'T'   =>{'A'=>1, 'T'=>1, 'C'=>1, 'G'=>1, $_gap=>0, $_null=>0},
		'C'   =>{'A'=>1, 'T'=>1, 'C'=>1, 'G'=>1, $_gap=>0, $_null=>0},
		'G'   =>{'A'=>1, 'T'=>1, 'C'=>1, 'G'=>1, $_gap=>0, $_null=>0},
		$_gap =>{'A'=>0, 'T'=>0, 'C'=>0, 'G'=>0, $_gap=>0, $_null=>0},
		$_null=>{'A'=>0, 'T'=>0, 'C'=>0, 'G'=>0, $_gap=>0, $_null=>0}
	);

}

sub _do_indels {
	warn("do_indels\n");

	%_m_mat = (   # mismatch matrix
		'A'   =>{'A'=>0, 'T'=>1, 'C'=>1, 'G'=>1, $_gap=>1, $_null=>0},
		'T'   =>{'A'=>1, 'T'=>0, 'C'=>1, 'G'=>1, $_gap=>1, $_null=>0},
		'C'   =>{'A'=>1, 'T'=>1, 'C'=>0, 'G'=>1, $_gap=>1, $_null=>0},
		'G'   =>{'A'=>1, 'T'=>1, 'C'=>1, 'G'=>0, $_gap=>1, $_null=>0},
		$_gap =>{'A'=>1, 'T'=>1, 'C'=>1, 'G'=>1, $_gap=>0, $_null=>0},
		$_null=>{'A'=>0, 'T'=>0, 'C'=>0, 'G'=>0, $_gap=>0, $_null=>0}
	);
	
	%_c_mat = (   # mismatch matrix
		'A'   =>{'A'=>1, 'T'=>1, 'C'=>1, 'G'=>1, $_gap=>1, $_null=>0},
		'T'   =>{'A'=>1, 'T'=>1, 'C'=>1, 'G'=>1, $_gap=>1, $_null=>0},
		'C'   =>{'A'=>1, 'T'=>1, 'C'=>1, 'G'=>1, $_gap=>1, $_null=>0},
		'G'   =>{'A'=>1, 'T'=>1, 'C'=>1, 'G'=>1, $_gap=>1, $_null=>0},
		$_gap =>{'A'=>1, 'T'=>1, 'C'=>1, 'G'=>1, $_gap=>0, $_null=>0},
		$_null=>{'A'=>0, 'T'=>0, 'C'=>0, 'G'=>0, $_gap=>0, $_null=>0}
	);

}


# ----------------------------------
# initialize new diversity calculation
# ----------------------------------
my @_read_buffer; # holds preprocessed read strings

sub initialize {
	my ($self, $infilename) = @_;

	my $seqio_obj = Bio::SeqIO->new(	-file =>   $infilename,
										-format=>  'fasta',
										-alphabet=>'dna'
										);
	# reset variables
	$_K     = 0;
	$_W     = 0;
	$_M     = 0;
	$_Z     = 0;
	$_D     = 0;
	$_varD  = 0;
	$_sigmaD= 0;
	
	@_p   =();
	@_m   =();
	@_z   =();
	@_var =();
	@_cov =();

	# load the read buffer with standardized reads
	@_read_buffer = ();	
	while(my $seqio_obj = $seqio_obj->next_seq) {
		my $seq_string = _standardize_the_read($seqio_obj);
		$_read_buffer[$_K] = $seq_string;
		$_K++;
	}

	_estimate_probabilities();
	
	return(1);
}


sub _standardize_the_read {
	# standardize one read 
	# *** warning this is for DNA alphabet only. It needs to be generalized ***
	# 1. replace N's with null symbol
	# 2. replace exterior gaps with the null symbol because exterior gaps are 
	#    past the limits of the read, while interior gaps are actual gaps. 

	my $seq_object = shift;
	
	my ($read_id) = ($seq_object->display_id =~/_(\w+)_/);
	my $seq_string = $seq_object->seq;
	
	$_W = length($seq_string); # width of the alignment (should be the same for all reads in the file)

	# 1) trim leading gaps (replace with null symbol)
	my ($leading_gaps) = ($seq_string =~ m/^(-+)/);
	if(defined($leading_gaps)) {
		$leading_gaps =~ s/-/$_null/g;
		$seq_string =~ s/^-+/$leading_gaps/;
	}

	# 2) trim trailing gaps (replace with null symbol)
	my ($trailing_gaps) = ($seq_string =~ m/(-+)$/);
	if(defined($trailing_gaps)) {
		$trailing_gaps =~ s/-/$_null/g;
		$seq_string =~ s/-+$/$trailing_gaps/;
	}
		
	# 3) replace everything that is not in @_nonnull_symbols with the null symbol
	$seq_string =~ s/[^ATCG-]/$_null/g; 
	
	return $seq_string;
}

sub _estimate_probabilities {

	my @frequency =();

	# zero out the accumulators
	for(my $i=0; $i<$_W; $i++) {
		foreach my $symbol ((@_alphabet)) {
			$frequency[$i]{$symbol} = 0;
			$_p[$i]{$symbol}  = 0;
		}
	}
	
	# accumulate symbol counts	
	for(my $k=0; $k<$_K; $k++) { 
		my @symbol_sequence = split //,$_read_buffer[$k];
		for(my $i=0; $i<$_W; $i++) {
			my $symbol = $symbol_sequence[$i];			
			$frequency[$i]{$symbol}++;
		}
	}
		
	# calc probabilities, variances and covariances
	for(my $i=0; $i<$_W; $i++) {
		foreach my $beta (@_alphabet) {
			my $p = $frequency[$i]{$beta}/$_K;
			$_p[$i]{$beta} = $p;
			$_var[$i]{$beta} = $p*(1-$p)/$_K;
			foreach my $alpha (@_alphabet) {
				if($alpha lt $beta) {
					$_cov[$i]{$alpha}{$beta} = -$_p[$i]{$alpha}*$_p[$i]{$beta}/$_K;
					$_cov[$i]{$beta}{$alpha} = $_cov[$i]{$alpha}{$beta};
				}
			}
		}
	}
	return 1;
}


sub _calculate_diversity {

	# calc mismatches
	for(my $i=0; $i<$_W; $i++) {
		foreach my $alpha (@_alphabet) {
			$_m[$i]{$alpha} = 0;
			foreach my $beta (@_alphabet) {
				$_m[$i]{$alpha}  += $_m_mat{$alpha}{$beta}*$_p[$i]{$beta};
			}
		}
	}
	
	$_M=0;
	for(my $i=0; $i<$_W; $i++) {
		foreach my $alpha (@_alphabet) {
			$_M += $_m[$i]{$alpha}*$_p[$i]{$alpha}
		}
	}
	
	# calc expected pairwise coverage
	for(my $i=0; $i<$_W; $i++) {
		foreach my $alpha (@_alphabet) {
			$_z[$i]{$alpha} = 0;
			foreach my $beta (@_alphabet) {
				$_z[$i]{$alpha}  += $_c_mat{$alpha}{$beta}*$_p[$i]{$beta};
			}
		}
	}
	
	$_Z=0;
	for(my $i=0; $i<$_W; $i++) {
		foreach my $alpha (@_alphabet) {
			$_Z += $_z[$i]{$alpha}*$_p[$i]{$alpha}
		}
	}

	# diversity
	$_D = $_M/$_Z;

	# ------------------
	# calc variance of the Diversity via error propagation
	# ------------------
	
	$_varD = 0;	
	for(my $i=0; $i<$_W; $i++) {
		my %Dp  = ();
		my %Dpc = ();
		foreach my $alpha (@_alphabet) {
			$Dp{$alpha} = 2*($_Z*$_m[$i]{$alpha} - $_M*$_z[$i]{$alpha})/($_Z*$_Z);

			my $sum = $Dp{$alpha}*$_var[$i]{$alpha};
			foreach my $beta (@_alphabet) {
				if($beta lt $alpha) {
					$sum += 2*$Dp{$beta}*$_cov[$i]{$alpha}{$beta};
				}
			}
			$_varD += $Dp{$alpha}*$sum;
		}
	}
	$_sigmaD = sqrt($_varD);
	
	return 1;
}

# ----------------------------------
# brute-force diversity (robust average pairwise difference)
# ----------------------------------

sub apd {

	my $m_sum = 0;
	my $c_sum = 0;
	for(my $i=0; $i<$_K; $i++) { 
		# warn("$i\n");
		my @seq_i = split //,$_read_buffer[$i];
		for(my $j=$i; $j<$_K; $j++) { 
			my @seq_j = split //,$_read_buffer[$j];
			for(my $w=0; $w<$_W; $w++) {
				$m_sum += $_m_mat{$seq_i[$w]}{$seq_j[$w]};
				$c_sum += $_c_mat{$seq_i[$w]}{$seq_j[$w]};
				#printf("\t%3s\t%3s\t%3d\n",
				#	$seq_i[$w], 
				#	$seq_j[$w],
				#	$_m_mat{$seq_i[$w]}{$seq_j[$w]}
				#	);
			}
		}
	}
	
	my $apd = $m_sum/$c_sum;
	# warn("m_sum = $m_sum\n");
	# warn("c_sum = $c_sum\n");
	# warn("m_sum/c_sum = $apd\n");
	return($apd);
}

# ----------------------------------
# getters
# ----------------------------------

sub epd {
	my $self=shift;
	
	_calculate_diversity();
	
	return ($_D, $_sigmaD);
}

sub diversity {
	my $self=shift;
	return $_D;
}

sub variance {
	my $self=shift;
	return $_varD;
}

sub sigma {
	my $self=shift;
	return $_sigmaD;
}

sub Z {
	my $self=shift;
	return $_Z;
}

sub n_reads {
	my $self=shift;
	return $_K
}

1;

__END__



# -----------------------------------------------------
# consensus sequences
# -----------------------------------------------------

sub consensus_alignment {
	my $self;
	return _calculate_aligned_consensus();
}

sub consensus {
	my $self;
	my $seq = _calculate_aligned_consensus();
	$seq =~ s/-//g;
	return $seq;
}


sub _argmax {
	my $probs = shift;

	my $argmax = ' '; 
	my $max   = -1;
	foreach my $symbol (@_alphabet) {
		if($probs->{$symbol} > $max) {
			$max = $probs->{$symbol};
			$argmax = $symbol;
		}
	}
	return $argmax;
}

sub _calculate_aligned_consensus {
	my $seq_string='';
	for(my $i=0; $i < $_width; $i++) {
		my $symbol;
#		if($_K[$i]==0) {
#			$symbol = 'N';
#		}
#		else {
			$symbol = _argmax($_p_hat[$i]);
#		}
		$seq_string .= $symbol;
	}
	return $seq_string;
}


