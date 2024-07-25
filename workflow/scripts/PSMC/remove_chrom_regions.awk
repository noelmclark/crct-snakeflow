# simple awk script that removes the y chr (or any chrom passed to chrom=) regions for all given
# scatter groups.  Then it prints the ones that pass in bed format If it picks out none, then it throws an error.

# use -v on invocation to set the value of scat to the chrom id


BEGIN {
	FS="\t"
	OFS="\t"
}

# check that the first row has columns in the correct order
NR == 1 {
	if($1 != "id") {
		kill = 1;
		print "First column of scatter_intervals file not named id" > "/dev/stderr"
	}
    if($2 != "scatter_idx") {
		kill = 1;
		print "Second column of scatter_intervals file not named scatter_idx" > "/dev/stderr"
	}
	if($3 != "chrom") {
		kill = 1;
		print "Third column of scatter_intervals file not named chrom" > "/dev/stderr"
	}
	if($4 != "start") {
		kill = 1;
		print "Fourth column of scatter_intervals file not named start" > "/dev/stderr"
	}
	if($5 != "end") {
		kill = 1;
		print "Fifth column of scatter_intervals file not named end" > "/dev/stderr"
	}
    if($6 != "scatter_length") {
		kill = 1;
		print "Sixth column of scatter_intervals file not named scatter_length" > "/dev/stderr"
	}
	if(kill > 0) {
		exit(1);
	}
	next
}

# all others just print out the three columns needed for the regions file that don't match the chrom listed
$3 != chrom {
	print $3, $4, $5;
	n++
}

END {
	if(n==0) {
		print "Didn't find any elements from chrom", chrom ". Bailing out!" > "/dev/stderr"
		exit(1);
	}
}