#!/usr/bin/awk -f

BEGIN {
  OFS = "\t"
  qc_file = "qc_skipped_rows.tsv"
  print "CHROM", "POS", "ID", "depth_a", "depth_b", "ratio", "num_hets", "num_samples", "num_called", "H_all", "H", "std", "D"
  print "CHROM", "POS", "ID", "REASON" > qc_file
}

FNR == NR {
  # Reading genotype file (first file)
  g_lines[FNR] = $0
  next
}

{
  current_line++
  
  if (!(current_line in g_lines)) {
    print "ERROR: Genotype file has fewer lines than depth file!" > "/dev/stderr"
    exit 1
  }

  # Process paired lines
  split(g_lines[current_line], g_fields, "\t")
  split($0, d_fields, "\t")

  chrom = g_fields[1]
  pos   = g_fields[2]
  id    = g_fields[3]

  hets = 0
  called = 0
  a_sum = 0
  b_sum = 0
  n_gt = length(g_fields)
  n_dp = length(d_fields)

  expected_samples = n_gt - 3
  if (n_dp != expected_samples) {
    print "ERROR: Sample count mismatch at line " current_line ": " \
          expected_samples " genotypes versus " n_dp " depth values" > "/dev/stderr"
    print chrom, pos, id, "Sample count mismatch" >> qc_file
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
  num_samples = expected_samples

  if (total_reads > 0) {
    ratio = a_sum / total_reads
    std = sqrt(total_reads * 0.25)

    if (std > 0) {
      D = -(total_reads / 2 - a_sum) / std
    } else {
      D = "NA"
      print "WARNING: std=0 at line " current_line " (" chrom ":" pos ") — cannot compute D" > "/dev/stderr"
      print chrom, pos, id, "std=0" >> qc_file
    }
  } else {
    ratio = "NA"
    std = "NA"
    D = "NA"
    print "WARNING: total_reads=0 at line " current_line " (" chrom ":" pos ") — skipping calculations" > "/dev/stderr"
    print chrom, pos, id, "total_reads=0" >> qc_file
  }

  H_all = hets / num_samples
  H = (called > 0) ? hets / called : "NA"

  print chrom, pos, id, a_sum, b_sum, ratio, hets, num_samples, called, H_all, H, std, D
}

END {
  if (length(g_lines) != current_line) {
    print "ERROR: Depth file has fewer lines than genotype file!" > "/dev/stderr"
    exit 1
  }
}
