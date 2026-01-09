#!/bin/bash
set -euo pipefail

cd tests

# test primary input
diff --color test_1_output.fastq \
    <(../flexiplex -k barcodes.txt -b GGGG -x TTT test_1_input.fastq 2>/dev/null)

# test segfault when no UMI is provided
../flexiplex -x ATCGGCGTACGACT -b \"????????\" -x ATCCACGTGCTTGAGACTGTGG -k test_23_barcodes.txt -f 2 -e 1 test_2_input.fastq 2>/dev/null >/dev/null

# test prefixed UMI - sample read
diff --color test_3_output.fastq \
    <(../flexiplex -u \"??????????\" -b \"????????\" -x GTGGCCGATGTTTCGCATCGGCGTACGACT -k test_23_barcodes.txt -f 4 -e 1 test_3_input.fastq 2>/dev/null)

rm flexiplex_reads_barcodes.txt