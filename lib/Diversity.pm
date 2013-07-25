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
my $_gap_null;		# gap and null symbols
my @_residues;		# residues array
my $_residues;		# residues string
my @_alphabet;		# alphabet
my $_residues_and_gap;   # all symbols except null

# ----------------------------------------
# diversity and variance calculations
# ----------------------------------------
my $_diversity_type; # indels, substitutions only, synonymous, etc.
my $_alphabet_type;  # 'dna' or 'aa' (amino acids)

my %_m_mat;	# mismatch matrix          $_m_mat{$symbol1}{$symbol2}
my %_c_mat;	# pairwise coverage matrix $_c_mat{$symbol1}{$symbol2}

my $_K;           # number of reads in the alignment file
my $_W;           # width (number of columns) of the alignment file
my $_M;     	# expected number of mismatches
my $_Z;		# expected pairwise width
my $_D;     	# diversity
my $_varD;		# variance of the diversity
my $_sigmaD;	# standard deviation of D

my @_p;		# symbol probabilities $_p[$i]{$symbol}
my @_m;
my @_z;
my @_var;   	# variances of the symbol probabilities
my @_cov;   	# covariance of the symbol probabilities
my @_valid_positions;  	# array positions in the alginment that can be included in
				# the diversity calculation
my $_gap_threshold;	# if proportion of gap symbols exceeds _gap_threshold, the position is no included in diversity 
my $_null_threshold;	# if proportion of null symbols exceeds _gap_threshold, the position is no included in diversity


# ----------------------------------
# new -- initialize symbols and matrices
# ----------------------------------

my %_options=();

sub new 
{ 
	my $self = shift;

	$_gap_threshold = 1;
	$_null_threshold = 1;

	
	return bless{};
}


# ----------------------------------
# initialize indicators
# ----------------------------------
sub _do_dna_subs {

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

sub _do_dna_indels {
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

sub _do_dna_synonymous {

	%_m_mat = (   # mismatch matrix
		'A'   =>{'A'=>0,'T'=>1,'C'=>1,'G'=>1,'a'=>0,'t'=>0,'c'=>0,'g'=>0,$_gap=>0,$_null=>0},
		'T'   =>{'A'=>1,'T'=>0,'C'=>1,'G'=>1,'a'=>0,'t'=>0,'c'=>0,'g'=>0,$_gap=>0,$_null=>0},
		'C'   =>{'A'=>1,'T'=>1,'C'=>0,'G'=>1,'a'=>0,'t'=>0,'c'=>0,'g'=>0,$_gap=>0,$_null=>0},
		'G'   =>{'A'=>1,'T'=>1,'C'=>1,'G'=>0,'a'=>0,'t'=>0,'c'=>0,'g'=>0,$_gap=>0,$_null=>0},
		'a'   =>{'A'=>0,'T'=>0,'C'=>0,'G'=>0,'a'=>0,'t'=>0,'c'=>0,'g'=>0,$_gap=>0,$_null=>0},
		't'   =>{'A'=>0,'T'=>0,'C'=>0,'G'=>0,'a'=>0,'t'=>0,'c'=>0,'g'=>0,$_gap=>0,$_null=>0},
		'c'   =>{'A'=>0,'T'=>0,'C'=>0,'G'=>0,'a'=>0,'t'=>0,'c'=>0,'g'=>0,$_gap=>0,$_null=>0},
		'g'   =>{'A'=>0,'T'=>0,'C'=>0,'G'=>0,'a'=>0,'t'=>0,'c'=>0,'g'=>0,$_gap=>0,$_null=>0},
		$_gap =>{'A'=>0,'T'=>0,'C'=>0,'G'=>0,'a'=>0,'t'=>0,'c'=>0,'g'=>0,$_gap=>0,$_null=>0},
		$_null=>{'A'=>0,'T'=>0,'C'=>0,'G'=>0,'a'=>0,'t'=>0,'c'=>0,'g'=>0,$_gap=>0,$_null=>0}
	);
	
	%_c_mat = (   # pairwise coverage matrix
		'A'   =>{'A'=>1,'T'=>1,'C'=>1,'G'=>1,'a'=>0,'t'=>0,'c'=>0,'g'=>0,$_gap=>0,$_null=>0},
		'T'   =>{'A'=>1,'T'=>1,'C'=>1,'G'=>1,'a'=>0,'t'=>0,'c'=>0,'g'=>0,$_gap=>0,$_null=>0},
		'C'   =>{'A'=>1,'T'=>1,'C'=>1,'G'=>1,'a'=>0,'t'=>0,'c'=>0,'g'=>0,$_gap=>0,$_null=>0},
		'G'   =>{'A'=>1,'T'=>1,'C'=>1,'G'=>1,'a'=>0,'t'=>0,'c'=>0,'g'=>0,$_gap=>0,$_null=>0},
		'a'   =>{'A'=>0,'T'=>0,'C'=>0,'G'=>0,'a'=>0,'t'=>0,'c'=>0,'g'=>0,$_gap=>0,$_null=>0},
		't'   =>{'A'=>0,'T'=>0,'C'=>0,'G'=>0,'a'=>0,'t'=>0,'c'=>0,'g'=>0,$_gap=>0,$_null=>0},
		'c'   =>{'A'=>0,'T'=>0,'C'=>0,'G'=>0,'a'=>0,'t'=>0,'c'=>0,'g'=>0,$_gap=>0,$_null=>0},
		'g'   =>{'A'=>0,'T'=>0,'C'=>0,'G'=>0,'a'=>0,'t'=>0,'c'=>0,'g'=>0,$_gap=>0,$_null=>0},
		$_gap =>{'A'=>0,'T'=>0,'C'=>0,'G'=>0,'a'=>0,'t'=>0,'c'=>0,'g'=>0,$_gap=>0,$_null=>0},
		$_null=>{'A'=>0,'T'=>0,'C'=>0,'G'=>0,'a'=>0,'t'=>0,'c'=>0,'g'=>0,$_gap=>0,$_null=>0}
	);

}


# ----------------------------------
# initialize new diversity calculation
# ----------------------------------
my @_read_buffer; # holds preprocessed read strings

sub initialize {
	my ($self, $infilename, $alphabet_type) = @_;

	my $seqio_obj;

	if(($alphabet_type eq 'dna')) {
		$_null     = 'N';
		$_gap      = '-';
		$_gap_null = $_gap.$_null;
		
		@_residues = ('A','T','C','G');

		$_residues = join '',@_residues;
		$_residues_and_gap = $_residues.$_gap;
		@_alphabet = sort (@_residues, $_gap, $_null);
		$seqio_obj = Bio::SeqIO->new(	-file =>   $infilename,
										-format=>  'fasta',
										-alphabet=>'dna'
							);
	}

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
	# and calculate the width ($_W) and number of rows ($_K) in the alignment file
	
	@_read_buffer = ();	
	while(my $seqio_obj = $seqio_obj->next_seq) {
		$_read_buffer[$_K] =  _standardize_the_read($seqio_obj);
		$_K++;
	}

	_estimate_probabilities();

	
	return(1);
}


sub _standardize_the_read {

	my $seq_object = shift;
	
	my ($read_id) = ($seq_object->display_id =~/_(\w+)_/);
	my $seq_string = $seq_object->seq;
	
	$_W = length($seq_string); # width of the alignment (should be the same for all reads in the file)

	# 1) trim leading nonresidue symbols with null symbol (because could be masked or padded with gaps)
	my ($leader) = ($seq_string =~ m/^([$_gap_null]+)/);
	if(defined($leader)) {
		my $len     = length($leader);
		$seq_string = ($_null x $len) . substr($seq_string, $len); 
	}

	# 2) trim trailing nonresidue symbols with null symbol (because could be masked or padded with gaps)
	my ($trailer) = ($seq_string =~ m/([$_gap_null]+)$/);
	if(defined($trailer)) {
		my $len     = length($trailer);
		$seq_string = substr($seq_string, 0, -$len) . ($_null x $len);
	}
		
	# 3) finally replace everything that is not a residue or a gap with the null symbol
	# $seq_string = uc $seq_string;
	# warn("commented out uc()");
	# $seq_string =~ s/[^$_residues_and_gap]/$_null/ig;

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

	# build array of alignment postions that are below the null threshold and below the gap threshold
	# these are the positions in the alignment that will be used to calculate the diversity
	@_valid_positions=();
	for(my $i=0; $i< $_W; $i++) {
		if(  ($_p[$i]{$_gap} <= $_gap_threshold) & ($_p[$i]{$_null} <= $_null_threshold)) {
			push @_valid_positions, $i;
		}
	}
	
	
	
	# calc mismatches
	foreach my $i (@_valid_positions) {
		foreach my $alpha (@_alphabet) {
			$_m[$i]{$alpha} = 0;
			foreach my $beta (@_alphabet) {
				$_m[$i]{$alpha}  += $_m_mat{$alpha}{$beta}*$_p[$i]{$beta};
			}
		}
	}
	
	$_M=0;
	foreach my $i (@_valid_positions) {
		foreach my $alpha (@_alphabet) {
			$_M += $_m[$i]{$alpha}*$_p[$i]{$alpha}
		}
	}
	
	# calc expected pairwise coverage
	foreach my $i (@_valid_positions) {
		foreach my $alpha (@_alphabet) {
			$_z[$i]{$alpha} = 0;
			foreach my $beta (@_alphabet) {
				$_z[$i]{$alpha}  += $_c_mat{$alpha}{$beta}*$_p[$i]{$beta};
			}
		}
	}
	
	$_Z=0;
	foreach my $i (@_valid_positions) {
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
	# for(my $i=0; $i<$_W; $i++) {
	foreach my $i (@_valid_positions) {
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
			foreach my $w (@_valid_positions) {
				$m_sum += $_m_mat{$seq_i[$w]}{$seq_j[$w]};
				$c_sum += $_c_mat{$seq_i[$w]}{$seq_j[$w]};
			}
		}
	}
	
	my $apd = $m_sum/$c_sum;
	return($apd);
}

# ----------------------------------
# setters and getters
# ----------------------------------

sub set_gap_threshold {
	my $self = shift;
	$_gap_threshold = shift;
}

sub set_null_threshold {
	my $self = shift;
	$_null_threshold = shift;
}

sub epd {
	my $self=shift;
	my %arguments = @_;
	if(!defined($arguments{'type'})) {
		_do_dna_subs();
	}
	elsif($arguments{'type'} eq 'indels') {
		_do_dna_indels();
	}
	elsif($arguments{'type'} eq 'no_indels') {
		_do_dna_subs();
	}
	elsif($arguments{'type'} eq 'synonymous') {
		_do_dna_synonymous();
	}

	
	
	my $_diversity_mode = shift;
	
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


