# lobSTR
allelotype \
  --command classify \
  --bam test.sorted.bam \
  --noise_model PATH_TO_LOBSTR/models/illumina_v3.pcrfree \
  --out test \
  --strinfo hg19_v3.0.2/lobstr_v3.0.2_hg19_strinfo.tab \
  --index-prefix hg19_v3.0.2/lobstr_v3.0.2_hg19_ref/lobSTR_