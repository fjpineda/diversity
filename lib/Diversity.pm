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

my $_null;			    # null symbol
my $_gap;			    # gap symbol
my $_end;               # symbol for leading or trailing gaps or nulls due to read length
my $_gap_null;		    # gap and null symbols
my @_residues;		    # residues array
my $_residues;		    # residues string
my @_alphabet;		    # alphabet array
my $_alphabet;			# alphabet string
my @_observed_symbols;	# all symbols found in the input file
my $_residues_and_gap;  # all symbols except null
my @_mask;	            # binary mask defining which positions to analyze
my %_alpha_mask;		# binary mask hash defining which symbols to compare for fast_apd

# ----------------------------------------
# diversity and variance calculations
# ----------------------------------------
my $_diversity_type; # including indels or substitutions only
my $_alphabet_type;  # 'dna', 'rna', or 'protein'
my %_m_mat;	# mismatch matrix                        $_m_mat{$symbol1}{$symbol2}
my %_c_mat;	# pairwise coverage matrix               $_c_mat{$symbol1}{$symbol2}

my $_K;         # number of reads in the alignment file
my $_W;         # width (number of columns) of the alignment file
my $_M;     	# expected number of mismatches
my $_Z;		    # expected pairwise width
my $_P;		    # number of pairs of symbols
my $_D;     	# diversity
my $_varD;		# variance of the diversity
my $_sigmaD;	# standard deviation of D
my @_freq;      # symbol frequencies in input file
my @_p;		    # symbol probabilities $_p[$i]{$symbol}
my @_var;   	# variances of the symbol probabilities
my @_cov;   	# covariance of the symbol probabilities
my @_valid_positions;  	# array positions in the alignment that can be included 
                        #   in the diversity calculation
my $_gap_threshold;	    # if proportion of gap symbols exceeds _gap_threshold, 
                        #   the position is not included in diversity calculations
my $_null_threshold;	# if proportion of null symbols exceeds _gap_threshold, 
                        #   the position is not included in diversity calculations
my $_end_threshold;	    # if proportion of end symbols exceeds _end_threshold, 
                        #   the position is not included in diversity calculations

# ----------------------------------
# new -- initialize symbols and matrices
# ----------------------------------

my %_options=();

sub new 
{ 
	my $self = shift;

	$_gap_threshold  = 1;
	$_null_threshold = 1;
	$_end_threshold  = 1;
	
	return bless{};
}


# ----------------------------------
# initialize indicators
# ----------------------------------
sub _do_subs {
	# warn("do_subs\n");		

	my $_m_mat_ref = _build_matrix('subs','mismatch'); # mismatch matrix
	%_m_mat = %$_m_mat_ref;
	my $_c_mat_ref = _build_matrix('subs','coverage'); # coverage matrix
	%_c_mat = %$_c_mat_ref;
	
	foreach my $res (@_residues) {$_alpha_mask{$res} = 1}
	@_alpha_mask{$_gap, $_null, $_end} = (0) x 3;
}

sub _do_indels {
	# warn("do_na_indels\n");

	my $_m_mat_ref = _build_matrix('indels','mismatch'); # mismatch matrix
	%_m_mat = %$_m_mat_ref;
	my $_c_mat_ref = _build_matrix('indels','coverage'); # coverage matrix
	%_c_mat = %$_c_mat_ref;
	
	foreach my $res (@_residues) {$_alpha_mask{$res} = 1}
	@_alpha_mask{$_gap, $_null, $_end} = (1,0,0);
}

sub _build_matrix {
	my ($div_type, $matrix_type) = @_;
	my %matrix;
	my $false_conditions;
	
	foreach my $row (@_alphabet) {
		my %cols;
    	foreach my $col (@_alphabet) {
    	
    		#define conditions under which each matrix is not true (zero)
    		if ($div_type eq 'subs' && $matrix_type eq 'mismatch') {
    			$false_conditions = $row eq $_null || $col eq $_null || $row eq $_gap 
    				|| $col eq $_gap || $row eq $col || $row eq $_end || $col eq $_end;
    		}
    		elsif ($div_type eq 'subs' && $matrix_type eq 'coverage') {
    			$false_conditions = $row eq $_null || $col eq $_null || $row eq $_gap 
    				|| $col eq $_gap || $row eq $_end || $col eq $_end;
    		}
    		elsif($div_type eq 'indels' && $matrix_type eq 'mismatch') {
    			$false_conditions = $row eq $_null || $col eq $_null || $row eq $col
    			 	|| $row eq $_end || $col eq $_end;
    		}
    		elsif($div_type eq 'indels' && $matrix_type eq 'coverage') {
    			$false_conditions = $row eq $_null || $col eq $_null || $row eq $_end 
    				|| $col eq $_end ||($row eq $_gap && $col eq $_gap);
    		}
    				
    		my $value = 1;
    		if ($false_conditions) { $value = 0 }
    		$cols{$col} = $value;
    	}
    	$matrix{$row} = \%cols;
	}
	
	#foreach my $row (keys %matrix) {
	#	print $row, "\t"; 
	#	foreach my $col (keys %{$matrix{$row}}) { 
	#		print $col, " ", $matrix{$row}->{$col}, "\t";
	#	}
	#	print "\n";
	#}
	#print "\n";

	return (\%matrix);
}


# ----------------------------------
# initialize new diversity calculation
# ----------------------------------
my @_read_buffer; # holds preprocessed read strings

sub initialize {
	(my $self, my $infilename, $_alphabet_type, my $mask) = @_;

	my $seqio_obj;

	if($_alphabet_type eq 'dna') {
		$_null     = 'N';
		@_residues = qw(A C G T);
	}
	elsif($_alphabet_type eq 'rna') {
		$_null     = 'N';
		@_residues = qw(A C G U);
	}
	elsif($_alphabet_type eq 'protein') {
		$_null     = 'X';
		@_residues = qw(A C D E F G H I K L M N P Q R S T V W Y *);
	}
	else {
		$_alphabet_type = 'dna';
		$_null     = 'N';
		@_residues = qw(A C G T);
	}	
	$_gap      = '-~.';
	$_end      = '#';
	$_gap_null = $_gap.$_null;
	$_residues = join '',@_residues;

	$seqio_obj = Bio::SeqIO->new(	-file =>   $infilename,
									-format=>  'fasta',
									-alphabet=>$_alphabet_type
						);

	# reset variables
	$_K     = 0;
	$_W     = 0;
	$_M     = 0;
	$_Z     = 0;
	$_P     = 0;	
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
				
	# detect which of the possible gap symbols is present
	foreach my $symbol (@_observed_symbols) {
		if($symbol =~ /[$_gap]/) { $_gap = $symbol }
	}
	
	$_gap_null = $_gap.$_null;
	$_residues_and_gap = $_residues.$_gap;
	@_alphabet = sort (@_residues, $_gap, $_null, $_end);
	$_alphabet = join("", @_alphabet);

	_qc_symbols();
	
	# populate the analysis mask array: 0 = skip; 1 = analyze; (default = analyze all positions)
	if ($mask) { @_mask = split(//,$mask) }
	else { @_mask = (1) x $_W }

	return(1);
}

sub _standardize_the_read {

	my $seq_object = shift;
	
	my ($read_id) = ($seq_object->display_id =~/_(\w+)_/);
	my $seq_string = $seq_object->seq;
	
	$_W = length($seq_string); # width of the alignment (should be the same for all reads in the file)

	# 1) trim leading nonresidue symbols with end symbol (because could be masked or padded with gaps)
	my ($leader) = ($seq_string =~ m/^([$_gap_null]+)/);
	if(defined($leader)) {
		my $len     = length($leader);
		$seq_string = ($_end x $len) . substr($seq_string, $len); 
	}

	# 2) trim trailing nonresidue symbols with end symbol (because could be masked or padded with gaps)
	my ($trailer) = ($seq_string =~ m/([$_gap_null]+)$/);
	if(defined($trailer)) {
		my $len     = length($trailer);
		$seq_string = substr($seq_string, 0, -$len) . ($_end x $len);
	}
	
	# Convert any lowercase alphabet symbols to uppercase
	$seq_string = uc $seq_string;	
	
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

sub _qc_symbols {
	
	# remove end symbol for printing
	my  (@no_end_alpha, @no_end_obs);
	for (@_observed_symbols) {unless ($_ eq $_end) {push(@no_end_obs, $_)}}
	for (@_alphabet) {unless ($_ eq $_end) {push(@no_end_alpha, $_)}}	

	# warn if alphabet type is 'protein' but input looks like nucleic acid
	if ($_alphabet_type eq 'protein' && @_observed_symbols <= length('ACGT$_gap$_null$_end')) {
		my $msg = "\nWarning: You have specified the alphabet type as 'protein', \n";
		$msg .=  "  but the input data may be dna or rna.\n"; 
		$msg .= "An incorrect alphabet will cause erroneous results.\n";
		$msg .= "Expected alphabet: @no_end_alpha\n";
		$msg .= "Found alphabet: @no_end_obs\n\n";
		warn($msg);
	} 
	
	# die if unexpected symbols are detected
	my $observed_symbols = join("", @_observed_symbols);
	if ($observed_symbols =~ /[^$_alphabet]/) {
		my $msg = "\nUnexpected symbols detected in the input file!\n";
		$msg .= "Expected alphabet: @no_end_alpha\n";
		$msg .= "Found alphabet: @no_end_obs\n\n";
		die($msg);
	}
	
	# fill-in with 0 the frequency of missing symbols in positions that don't 
    # have the full complement of symbols
	for(my $i=0; $i<$_W; $i++) {
		foreach my $symbol (@_observed_symbols) {
			unless(defined($_freq[$i]{$symbol})) {
				$_freq[$i]{$symbol} = 0; 
			}
		}
	}

	# fill-in with 0's the frequency of any alphabet symbols that did not appear
	# in the input file
	for(my $i=0; $i<$_W; $i++) {
		#foreach my $alpha (@_residues) {
		foreach my $alpha (@_alphabet) {
			if(!defined($_freq[$i]{$alpha})) {
				$_freq[$i]{$alpha} = 0;
			}
		}	
	}	

	return(1);
}

sub calculate_diversity {
	my ($self, $_diversity_algorithm, $_diversity_type) = @_;
		
	# build array of alignment positions that are below the null, gap, and end thresholds 
	#   and are not masked
	# these are the positions in the alignment that will be used to calculate the diversity
	@_valid_positions=();
	for(my $i=0; $i< $_W; $i++) {
		if(  ($_freq[$i]{$_gap}   <= $_gap_threshold * $_K) 
			& ($_freq[$i]{$_null} <= $_null_threshold * $_K)
			& ($_freq[$i]{$_end}  <= $_end_threshold * $_K)			
			& $_mask[$i]) {
			push @_valid_positions, $i;
		}
	}
	
	# if there are no valid positions to analyze return undef	
	unless (@_valid_positions) {
		return undef;
	}

	# Build matrices for the designated type of calculation
	if(!defined($_diversity_type)) {
		$_diversity_type = 'no_indels';
		_do_subs();
	}
	elsif($_diversity_type eq 'with_indels') {
		_do_indels;
	}
	elsif($_diversity_type eq 'no_indels') {
		_do_subs();
	}
	else {
		die("unknown diversity type: '$_diversity_type'");
	}

	# Calculate pairwise difference with the appropriate algorithm
	if(!defined($_diversity_algorithm)) {
		$_diversity_algorithm = 'epd';
		_epd();
	}
	elsif($_diversity_algorithm eq 'epd') {
		_epd();
	}
	elsif($_diversity_algorithm eq 'fast_apd') {
		_fast_apd();
	}
	elsif($_diversity_algorithm eq 'fast_apd2') {
		_fast_apd2();
	}
	elsif($_diversity_algorithm eq 'apd') {
		_apd();
	}
	elsif($_diversity_algorithm eq 'conventional_apd') {
		_conventional_apd();
	}
	else {
		die("unknown diversity algorithm: '$_diversity_algorithm'");
	}
	
	return($_D, $_sigmaD);
}
	
# ----------------------------------
# Expected pairwise difference calculation
# ----------------------------------	
	
sub _epd {
	
	# calculate probabilities, variances and covariances
	@_p=();
	@_var=();
	@_cov=();
	for(my $i=0; $i<$_W; $i++) {
		foreach my $beta (@_alphabet) {
			if(!defined($_freq[$i]{$beta})){
				warn("undefined element for i = $i and beta = $beta");
			}
		
			my $p = $_freq[$i]{$beta}/$_K;
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
# Combinatorial average pairwise difference calculation
# ----------------------------------	
	
sub _fast_apd {
	my ($sum_freqs, $sum_matches, $_P) = (0) x 3;
	my %matches;	
	
	# sum the number of matches and pairs
	foreach my $i (@_valid_positions) {
		foreach my $alpha (@_alphabet) {
			$sum_freqs += $_freq[$i]{$alpha}*$_alpha_mask{$alpha};
			$matches{$alpha} += _choose_2($_freq[$i]{$alpha}*$_alpha_mask{$alpha});
		}
		$_P += _choose_2($sum_freqs);
		$sum_freqs = 0;
	}
	
	# adjust for gap to gap comparisons
	$_P -= $matches{$_gap};
	foreach my $res (@_residues) {$sum_matches += $matches{$res}}
	
	# mismatches
	$_M = $_P - $sum_matches;

	# average pairwise diversity
	$_D = $_M/$_P;
	
	return 1;
}

# ----------------------------------
# brute-force diversity (robust average pairwise difference)
# ----------------------------------

sub _apd {
	my $m_sum = 0;
	my $c_sum = 0;
	for(my $i=0; $i<$_K; $i++) { 
		my @seq_i = split //,$_read_buffer[$i];
		for(my $j=$i+1; $j<$_K; $j++) { 
			my @seq_j = split //,$_read_buffer[$j];
			for(my $w=0; $w< $_W; $w++) {
				$m_sum += $_m_mat{$seq_i[$w]}{$seq_j[$w]};
				$c_sum += $_c_mat{$seq_i[$w]}{$seq_j[$w]};
			}
		}
	}
	
	$_D = $m_sum/$c_sum;
	return 1;
}

# ----------------------------------
# conventional brute-force average pairwise difference
# ----------------------------------

sub _conventional_apd {

	my ($m_sum, $c_sum, $sum, $read_pairs) = (0) x 4;
	for(my $i=0; $i<$_K; $i++) { 
		my @seq_i = split //,$_read_buffer[$i];
		for(my $j=$i+1; $j<$_K; $j++) { 
			my @seq_j = split //,$_read_buffer[$j];
			for(my $w=0; $w< $_W; $w++) {
				$m_sum += $_m_mat{$seq_i[$w]}{$seq_j[$w]};
				$c_sum += $_c_mat{$seq_i[$w]}{$seq_j[$w]};
			}
			if ($c_sum) {
				$sum += $m_sum/$c_sum;
				$read_pairs++;
				($m_sum, $c_sum) = (0) x 2;
			}
		}
	}
	
	$_D = $sum/$read_pairs;
	return 1;
}

# ----------------------------------
# choose function
# ----------------------------------

sub _choose_2 {
	my $n = shift;
	return($n * ($n - 1) / 2);
}


# ----------------------------------
# setters and getters
# ----------------------------------

sub gap_threshold {
	my $self = shift;
	my $set_thresh = shift;
	if ($set_thresh) {$_gap_threshold = $set_thresh}
	return($_gap_threshold);
}

sub null_threshold {
	my $self = shift;
	my $set_thresh = shift;
	if ($set_thresh) {$_null_threshold = $set_thresh}
	return($_null_threshold);
}

sub end_threshold {
	my $self = shift;
	my $set_thresh = shift;
	if ($set_thresh) {$_end_threshold = $set_thresh}
	return($_end_threshold);
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
	return $_K;
}

sub valid_positions {
	my $self=shift;
	my $size = @_valid_positions;
	return $size;
}

sub width {
	my $self=shift;
	return $_W;
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


