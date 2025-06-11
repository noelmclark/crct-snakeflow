#!/usr/bin/awk -f

BEGIN {
  OFS = "\t"
  file_index = 0
  current_line = 0
  print "CHROM", "POS", "ID", "depth_a", "depth_b", "ratio", "num_hets", "num_samples", "num_called", "H_all", "H", "std", "D"
}

FNR == NR {
  # Reading genotype file (first file)
  file_index = 1
  g_lines[FNR] = $0
  next
}

{
  # Reading depth file (second file)
  current_line++

  if (!(current_line in g_lines)) {
    print "ERROR: Genotype file has fewer lines than depth file!" > "/dev/stderr"
    exit 1
  }

  # Process paired lines
  split(g_lines[current_line], g_fields, "\t")
  split($0, d_fields, "\t")

  chrom = g_fields[1]
  pos = g_fields[2]
  id = g_fields[3]

  hets = 0
  called = 0
  a_sum = 0
  b_sum = 0
  n_gt = length(g_fields)
  n_dp = length(d_fields)

  expected_samples = n_gt - 3
  if (n_dp != expected_samples) {
    print "ERROR: Sample count mismatch at line " current_line ": " \
          expected_samples " genotypes vs " n_dp " depth values" > "/dev/stderr"
    exit 1
  }

  for (i = 1; i <= expected_samples; i++) {
    gt = g_fields[i + 3]
    split(d_fields[i], ab, ",")
    a = ab[1] + 0
    b = ab[2] + 0

    if (gt == "0/1" || gt == "1/0") {
      hets++
      a_sum += a
      b_sum += b
    }

    if (gt != "./." && gt != ".") {
      called++
    }
  }

  total_reads = a_sum + b_sum
  ratio = (total_reads > 0) ? a_sum / total_reads : "NA"
  std = (total_reads > 0) ? sqrt(total_reads * 0.25) : "NA"
  D = (std > 0) ? -(total_reads / 2 - a_sum) / std : "NA"
  H_all = hets / expected_samples
  H = (called > 0) ? hets / called : "NA"

  print chrom, pos, id, a_sum, b_sum, ratio, hets, expected_samples, called, H_all, H, std, D
}

END {
  if (length(g_lines) != current_line) {
    print "ERROR: Depth file has fewer lines than genotype file!" > "/dev/stderr"
    exit 1
  }
}
