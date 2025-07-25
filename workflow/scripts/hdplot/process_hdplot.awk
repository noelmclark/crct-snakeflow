#!/usr/bin/awk -f

BEGIN {
  OFS = "\t"
  print "CHROM", "POS", "ID", "depth_a", "depth_b", "ratio", "num_hets", "num_samples", "num_called", "H_all", "H", "std", "D"
}

# Read the genotype file (first file)
FNR == NR {
  for (i = 4; i <= NF; i++) {
    gt[FNR, i] = $i
  }
  chrom[FNR] = $1
  pos[FNR] = $2
  id[FNR] = $3
  next
}

# Process the depth file (second file)
{
  n = split($0, ad, "\t")
  num_samples = n
  a_sum = 0
  b_sum = 0
  hets = 0
  called = 0

  for (i = 1; i <= n; i++) {
    split(ad[i], ab, ",")
    a = ab[1] + 0
    b = ab[2] + 0

    g = gt[FNR, i + 3]
    if (g == "0/1" || g == "1/0") {
      hets++
      a_sum += a
      b_sum += b
    }

    if (g != "./." && g != ".") {
      called++
    }
  }

  total_reads = a_sum + b_sum
  ratio = (total_reads > 0) ? a_sum / total_reads : "NA"
  std = (total_reads > 0) ? sqrt(total_reads * 0.25) : "NA"
  D = (std > 0) ? -(total_reads / 2 - a_sum) / std : "NA"
  H_all = hets / num_samples
  H = (called > 0) ? hets / called : "NA"

  print chrom[FNR], pos[FNR], id[FNR], a_sum, b_sum, ratio, hets, num_samples, called, H_all, H, std, D
}
