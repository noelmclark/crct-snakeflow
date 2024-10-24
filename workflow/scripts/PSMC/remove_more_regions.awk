# simple awk script that removes the y chr (or any chrom passed to chrom=) regions for all given
# scatter groups.  Then it prints the ones that pass in bed format If it picks out none, then it throws an error.

# use -v on invocation to set the value of scat to the chrom id


BEGIN {
	FS="\t"
	OFS="\t"
}


# all others just print out the three columns needed for the regions file that don't match the chrom listed
$1 != chrom {
	print $1, $2, $3;
	n++
}

END {
	if(n==0) {
		print "Didn't find any elements from chrom", chrom ". Bailing out!" > "/dev/stderr"
		exit(1);
	}
}