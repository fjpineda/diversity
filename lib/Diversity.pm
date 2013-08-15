package Diversity;

use strict;
use warnings;

use File::Basename;
use Data::Dumper;
use Bio::SeqIO;
# use Benchmark;

# ----------------------------------------
# debug
# ----------------------------------------

sub _debug_freq {
	my $freq_ref = shift;
	# print Dumper($freq_ref); die();

	warn("\$freq_ref = $freq_ref\n");
	
	my %symbols=();
	foreach my $row (@$freq_ref) {
		# print Dumper($row);
		foreach my $symbol (keys %$row) {
			$symbols{$symbol} += $row->{$symbol};
		}
	}
	print Dumper(\%symbols);
}

# ----------------------------------------
# symbols & matrices
# ----------------------------------------

my $_null;			# null symbol
my $_gap;			# gap symbol
my $_gap_null;		# gap and null symbols
my @_residues;		# residues array
my $_residues;		# residues string
my @_alphabet;		# alphabet
my @_observed_symbols;	# all symbols found in the input file
my $_residues_and_gap;  # all symbols except null

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

my @_freq;        # symbol frequencies in input file
my @_p;		# symbol probabilities $_p[$i]{$symbol}
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

	# warn("do_subs\n");

	@_residues = ('A','T','C','G');

	$_residues = join '',@_residues;
	$_residues_and_gap = $_residues.$_gap;
	@_alphabet = sort (@_residues, $_gap, $_null);
		

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
	# warn("do_indels\n");

	@_residues = ('A','T','C','G');

	$_residues = join '',@_residues;
	$_residues_and_gap = $_residues.$_gap;
	@_alphabet = sort (@_residues, $_gap, $_null);

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
	
	# warn("do_synonymous\n");

	@_residues = ('A','T','C','G','a','t','c','g');

	$_residues = join '',@_residues;
	$_residues_and_gap = $_residues.$_gap;
	@_alphabet = sort (@_residues, $_gap, $_null);

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
	@_var =();
	@_cov =();
	
	# load the read buffer with standardized reads
	# and calculate the width ($_W) and number of rows ($_K) in the alignment file
	
	@_read_buffer = ();	
	while(my $seqio_obj = $seqio_obj->next_seq) {
		$_read_buffer[$_K] =  _standardize_the_read($seqio_obj);
		$_K++;
	}

	_accumulate_symbol_frequencies();

	
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
		
	return $seq_string;
}

sub _accumulate_symbol_frequencies {
	# warn("_accumulate_symbol_frequencies()\n");

	@_freq =();
	my %counts=();

	# initialize the accumulators
	for(my $i=0; $i<$_W; $i++) {
		$_freq[$i]=undef;
	}
	
	# accumulate symbol counts	
	for(my $k=0; $k<$_K; $k++) { 
		my @symbol_sequence = split //,$_read_buffer[$k];
		for(my $i=0; $i<$_W; $i++) {
			my $symbol = $symbol_sequence[$i];			
			$_freq[$i]{$symbol}++;
			$counts{$symbol}++;
		}
	}
	@_observed_symbols = keys %counts;
	
	return 1;
}


sub _calculate_diversity {

	# copy the array of frequencies and fill-in with 0 the frequency of missing
	# symbols in positions that don't have the full complement of symbols
	my @frequency =();	
	for(my $i=0; $i<$_W; $i++) {
		foreach my $symbol (@_observed_symbols) {
			if(defined($_freq[$i]{$symbol})) {
				$frequency[$i]{$symbol} = $_freq[$i]{$symbol};
			}
			else{
				$frequency[$i]{$symbol} = 0; 
			}
		}
	}

	# fill-in with 0's the frequency of any alphabet symbols that did not appear
	# in the input file
	for(my $i=0; $i<$_W; $i++) {
		foreach my $alpha (@_residues) {
			if(!defined($frequency[$i]{$alpha})) {
				$frequency[$i]{$alpha} = 0;
			}
		}	
	}	
	
	# if the input file has lower and upper case symbols collapse them all to 
	# upper case, except if we are calculating synonymous diversity in which case
	# the lower case symbols mark the amino acid residues with synonymous mutations
	if($_diversity_type ne 'syn') {
		for(my $i=0; $i<$_W; $i++) {
			foreach my $alpha (@_residues) {
				my $lower_case_alpha = lc($alpha);
				if(defined($frequency[$i]{$lower_case_alpha})) {
					$frequency[$i]{$alpha} += $frequency[$i]{$lower_case_alpha};
					delete($frequency[$i]{$lower_case_alpha});
				}

			}
		}
	}

	
	for(my $i=0; $i<$_W; $i++) {
		if(!defined($_freq[$i]{$_gap})) {
				$frequency[$i]{$_gap}=0;	
		}
		if(!defined($_freq[$i]{$_null})) {
				$frequency[$i]{$_null}=0;	
		}
	}
	
	# _debug_freq(\@frequency); 
	
	# calculate probabilities, variances and covariances
	@_p=();
	@_var=();
	@_cov=();
	for(my $i=0; $i<$_W; $i++) {
		foreach my $beta (@_alphabet) {
			if(!defined($frequency[$i]{$beta})){
				warn("undefined element for i = $i and beta = $beta");
			}
		
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

	# build array of alignment postions that are below the null threshold and below the gap threshold
	# these are the positions in the alignment that will be used to calculate the diversity
	@_valid_positions=();
	for(my $i=0; $i< $_W; $i++) {
		if(  ($_p[$i]{$_gap} <= $_gap_threshold) & ($_p[$i]{$_null} <= $_null_threshold)) {
			push @_valid_positions, $i;
		}
	}
	


	# calc mismatches
	my @Sm   =();
	foreach my $i (@_valid_positions) {
		foreach my $alpha (@_alphabet) {
			$Sm[$i]{$alpha} = 0;
			foreach my $beta (@_alphabet) {
				$Sm[$i]{$alpha}  += $_m_mat{$alpha}{$beta}*$_p[$i]{$beta};
			}
		}
	}
	
	$_M=0;
	foreach my $i (@_valid_positions) {
		foreach my $alpha (@_alphabet) {
			$_M += ($Sm[$i]{$alpha}-$_m_mat{$alpha}{$alpha}/$_K)*$_p[$i]{$alpha}
		}
	}
	$_M *= $_K/($_K-1);
	
	# calc expected pairwise coverage
	my @Sc=();
	foreach my $i (@_valid_positions) {
		foreach my $alpha (@_alphabet) {
			$Sc[$i]{$alpha} = 0;
			foreach my $beta (@_alphabet) {
				$Sc[$i]{$alpha}  += $_c_mat{$alpha}{$beta}*$_p[$i]{$beta};
			}
		}
	}
	
	$_Z=0;
	foreach my $i (@_valid_positions) {
		foreach my $alpha (@_alphabet) {
			$_Z += ($Sc[$i]{$alpha}-$_c_mat{$alpha}{$alpha}/$_K)*$_p[$i]{$alpha}
		}
	}
	$_Z *= $_K/($_K-1);

	# diversity
	$_D = $_M/$_Z;

	# ------------------
	# calc variance of the Diversity via error propagation
	# ------------------
	
	$_varD = 0;	
	foreach my $i (@_valid_positions) {
		my %Dp = ();
		foreach my $alpha (@_alphabet) {
			my $Mp = ($_K/($_K-1) ) *(2*$Sm[$i]{$alpha}-$_m_mat{$alpha}{$alpha}/$_K);
			my $Zp = ($_K/($_K-1) ) *(2*$Sc[$i]{$alpha}-$_c_mat{$alpha}{$alpha}/$_K);
			$Dp{$alpha} = ($_Z*$Mp - $_M*$Zp)/($_Z*$_Z);
		}
		foreach my $alpha (@_alphabet) {
			my $sum = $Dp{$alpha}*$_var[$i]{$alpha};
			foreach my $beta (@_alphabet) {
				if($beta lt $alpha) {
					$sum += 2*$Dp{$beta}*$_cov[$i]{$alpha}{$beta};
				}
			}
			$_varD += $Dp{$alpha} *$sum;
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
		for(my $j=$i+1; $j<$_K; $j++) { 
			my @seq_j = split //,$_read_buffer[$j];
			for(my $w=0; $w< $_W; $w++) {
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
	$_diversity_type = $arguments{'type'};




	if(!defined($_diversity_type)) {
		$_diversity_type = 'no_indels';
		_do_dna_subs();
	}
	elsif($_diversity_type eq 'with_indels') {
		_do_dna_indels();
	}
	elsif($_diversity_type eq 'no_indels') {
		_do_dna_subs();
	}
	elsif($_diversity_type eq 'syn') {
		_do_dna_synonymous();
	}
	else {
		die("unknown diversity type: '$_diversity_type'");
	}
	
	_calculate_diversity();
	
	return ($_D, $_sigmaD, $_Z);
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


