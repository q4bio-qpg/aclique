#!/usr/bin/perl -w

# Load graph
my $p=<>;
my @v;
my @e;
my $nvertex = 0;
my $nedge = 0;
my @vmap; # map weighted to first unweighted node
while (<>) {
    if (/^v\s+(\d+)\s+(\d+)/) {
	$v[$1]=$2;
	$vmap[$1]=$nvertex+1;
	$nvertex += $2;
    }
    if (/^e\s+(\d+)\s+(\d+)/) {
	if ($2 > $1) {
	    push(@{$e[$1]},$2);
	    $nedge++;
	}
    }
}

my @e2;
for (my $v=1; $v < scalar(@v); $v++) {
    my $u = $vmap[$v];
    for (my $u1 = $vmap[$v]; $u1<$vmap[$v]+$v[$v]; $u1++) {
	for (my $u2 = $u1+1; $u2<$vmap[$v]+$v[$v]; $u2++) {
		push(@e2, "e $u1 $u2");
	}
	foreach (@{$e[$v]}) {
	    for (my $u2 = $vmap[$_]; $u2<$vmap[$_]+$v[$_]; $u2++) {
		push(@e2, "e $u1 $u2");
	    }
	}
    }
}

print "p edge $nvertex ",scalar(@e2),"\n";
print join("\n",@e2),"\n";

