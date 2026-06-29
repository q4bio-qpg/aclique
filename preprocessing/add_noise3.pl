#!/usr/bin/perl -w
use strict;

# Adds random noise to columns, simulating changes that are non-inherited.
# Eg homoplasy - convergent evolution rather than inherited.

# Usage add_noise fraction seed count file.fa
# Produces file.fa.<n> for count consecutive <n>s (0 upwards)

my $fract = shift(@ARGV);
srand(shift(@ARGV));
my $count = shift(@ARGV);
my $outfn = $ARGV[0];

# Orig =                    4-5min
# 1e-1, correlation 0.9  => 2m42s  RF=1   SCORE=+1   switches, +0  siblings s=1
# 1e-1, correlation 0.95 => 2m31s  RF=1   SCORE=+0   switches, +0  siblings s=1
# 1e-1, correlation 0.99 => 0m13s  RF=53  SCORE=+20  switches, +27 siblings s=1
						     
# 1e-1, correlation 1.00 => 0m2.4s RF=1   SCORE=+1   switches, +0  siblings s=1
# 1e-1, correlation 1.00 => 1m33s  RF=103 SCORE=+29  switches, +73 siblings s=2
# 1e-1, correlation 1.00 => 0m2.4s RF=1   SCORE=+1   switches, +0  siblings s=3
# 1e-1, correlation 1.00 => 0m2.0s RF=1   SCORE=+1   switches, +0  siblings s=4
# 1e-1, correlation 1.00 => 0m1.9s RF=100 SCORE=+39  switches, +76 siblings s=5

my $correlation = 0.99;#0.95; # eg 0.7 => 70% correlated to the last variant

my @bases = qw/A C G T/;

# Load seqs
my @name;
my @seq;
while (<>) {
    chomp($_);
    if (/^>/) {
	push(@name, $_);
    } else {
	push(@seq, $_);
    }
}

# Pick columns to add, sorted in descending order to make adding easy
my $len = length($seq[0]);
my $nedits = int($len * $fract);
my @pos;   # position of bases
my @base1; # base 1 at pos
my @base2; # base 2 at pos
my @freq;  # probability of base1 (vs base2)

for (my $i = 0; $i < $nedits; $i++) {
    $pos[$i] = int(rand($len));
    $base1[$i] = $bases[int(rand(4))];
    do {
	$base2[$i] = $bases[int(rand(4))];
    } while $base2[$i] eq $base1[$i];
    $freq[$i] = rand();
}
@pos = sort {$b <=> $a} @pos;

for (my $n = 0; $n < $count; $n++) {
    open(FH, ">", $outfn . "." . $n) || die;

    my @last_allele = ();
    for (my $i = 0; $i < $nedits; $i++) {
	my @allele;
	if ($i > 0 && rand() < $correlation) {
	    @allele = @last_allele;
	} else {
	    for (my $snum = 0; $snum < scalar(@name); $snum++) {
		$allele[$snum] = rand() < $freq[$i] ? 1 : 0;
	    }
	}
	for (my $snum = 0; $snum < scalar(@name); $snum++) {
	    my $b = $allele[$snum] ? $base1[$i] : $base2[$i];
	    substr($seq[$snum], $pos[$i], 0) = $b;
	}

	@last_allele = @allele;
    }

    for (my $snum = 0; $snum < scalar(@name); $snum++) {
	my $name = $name[$snum];
	my $seq  = $seq[$snum];
	print FH "$name\n$seq\n";
    }

    close(FH);
}
