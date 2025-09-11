# simple awk script that pulls the first two columns of the given tsv file and saves it to a new file


BEGIN {
	FS="\t"
	OFS="\t"
}

# check that the first row has columns in the correct order
NR == 1 {
	if($1 != "chrom") {
		kill = 1;
		print "First column of file not named chrom" > "/dev/stderr"
	}
    if($2 != "pos") {
		kill = 1;
		print "Second column of file not named pos" > "/dev/stderr"
	}
	if(kill > 0) {
		exit(1);
	}
	next
}

# print out the first two columns needed for the targets file 
{ print $1, $2; n++ }

END {
	if(n==0) {
		print "Didn't find any columns. Bailing out!" > "/dev/stderr"
		exit(1);
	}
}