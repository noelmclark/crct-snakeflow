## I copied this from Eric's post bcf snakeflow https://github.com/noelmclark/mega-post-bcf-exploratory-snakeflows/blob/ebf68d44d45b20ac4f5a2e75eb2a4eb60915b5ee/workflow/scripts/beagle3header.awk
# simple awk script.  Pass it a file with sample names
# one per line and it makes the Beagle3 header


BEGIN {printf("marker\tallele1\tallele2")} 

{for(i=1;i<=3;i++) printf("\t%s", $1)} 

END {printf("\n")}