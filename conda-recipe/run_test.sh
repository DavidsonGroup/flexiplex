set -eu
msg() {
  echo "=> $@"
}

shell() {
  msg Testing bash curl etc
  (
  diff <(curl -sL https://raw.githubusercontent.com/mritchielab/FLAMES/devel/inst/extdata/bc_allow.tsv.gz | zcat) \
    <(curl -sL https://raw.githubusercontent.com/mritchielab/FLAMES/devel/inst/extdata/bc_allow.tsv.gz | zcat) > /dev/null
)
}

testflexiplex() {
  msg Testing flexiplex demultiplex output 
  diff <(flexiplex -x CTACACGACGCTCTTCCGATCT -b '????????????????' -u '????????????' -x TTTTTTTTT -e 2 -f 8 -k <(curl -sL https://raw.githubusercontent.com/mritchielab/FLAMES/devel/inst/extdata/bc_allow.tsv.gz | zcat) <(curl -sL https://raw.githubusercontent.com/mritchielab/FLAMES/devel/inst/extdata/fastq/musc_rps24.fastq.gz | zcat)) <(curl -sL https://raw.githubusercontent.com/mritchielab/FLAMES/devel/tests/testthat/demultiplexed.fq) > /dev/null
  msg done
}

testflexiplexfail() {
  msg Testing fail
  diff <(flexiplex -x CTACACGACGCTCTTCCGATCT -b '????????????????' -u '????????????' -x TTTTTTTTT -e 2 -f 1 -k <(curl -sL https://raw.githubusercontent.com/mritchielab/FLAMES/devel/inst/extdata/bc_allow.tsv.gz | zcat) <(curl -sL https://raw.githubusercontent.com/mritchielab/FLAMES/devel/inst/extdata/fastq/musc_rps24.fastq.gz | zcat) 2>/dev/null) <(curl -sL https://raw.githubusercontent.com/mritchielab/FLAMES/devel/tests/testthat/demultiplexed.fq) > /dev/null
  msg unexpected line reached
}


msg Running in $PWD
if [[ -n "${1:-}" ]]; then
  # Run single specified code-quality tool.
  $1
else
  # Run all code-quality tools.
  shell
  testflexiplex
  testflexiplexfail
fi
