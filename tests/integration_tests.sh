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

# test -a: an assigned read is reported exactly as it is without -a, with the
# uncorrected barcode and UMI added to the comment as CR and UR
diff --color test_1_output.fastq \
    <(../flexiplex -k barcodes.txt -b GGGG -x TTT -a true test_1_input.fastq 2>/dev/null \
      | sed 's/\tCR:Z:[^\t]*//; s/\tUR:Z:[^\t]*//')

# test -a: a read whose barcode is not in the known list is still reported, once,
# with '-' for the barcode and the barcode as observed in the read in CR
diff --color \
    <(printf '@-_TAAGACGTAT#SRR13948564.233_-1of1\tCB:Z:-\tCR:Z:AATCCGTC\tUB:Z:TAAGACGTAT\tUR:Z:TAAGACGTAT\n') \
    <(../flexiplex -u '??????????' -b '????????' -x GTGGCCGATGTTTCGCATCGGCGTACGACT -k CCCCCCCC -f 4 -e 0 -a true test_3_input.fastq 2>/dev/null \
      | awk 'NR%4==1')

# test -a: that same read is not reported at all without -a
unassigned=$(../flexiplex -u '??????????' -b '????????' -x GTGGCCGATGTTTCGCATCGGCGTACGACT -k CCCCCCCC -f 4 -e 0 test_3_input.fastq 2>/dev/null | wc -l)
if [ "$unassigned" -ne 0 ]; then
    echo "a read without a barcode match was reported without -a" >&2
    exit 1
fi

rm flexiplex_reads_barcodes.txt