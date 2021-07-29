# HPP SV Calling

This repository holds Docker build scripts for SV calling
used by the [Human Pangenome Project](https://humanpangenome.org)
and benchmark of different SV callers.

# Comparison

We compare our callsets with the GIAB v0.6 Tier 1 SV benchmark set
for HG002.

## Run pbsv

Docker image: `wwliao/hpp_pbsv:1.0.0--5ac93f8f6bd021a33feeb3fbca50466f27847b15`

```sh
$ pbsv discover -s HG002 \
                --tandem-repeats GRCh38_no_alt.trf.bed \
                hifi_alignments/HG002.GRCh38_no_alt.bam \
                calls/pbsv/outputs/HG002.GRCh38_no_alt.svsig.gz \
  && pbsv call --ccs \
               --preserve-non-acgt \
               -t DEL,INS,INV,DUP,BND \
               -m 40 \
               -j 12 \
               GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
               calls/pbsv/outputs/HG002.GRCh38_no_alt.svsig.gz \
               calls/pbsv/HG002.GRCh38_no_alt.pbsv.vcf \
  && bcftools sort -m4G \
                   -Oz \
                   -o calls/pbsv/HG002.GRCh38_no_alt.pbsv.vcf.gz \
                   calls/pbsv/HG002.GRCh38_no_alt.pbsv.vcf \
  && bcftools index -t calls/pbsv/HG002.GRCh38_no_alt.pbsv.vcf.gz \
  && rm calls/pbsv/HG002.GRCh38_no_alt.pbsv.vcf
```

## Run Sniffles and Iris

Docker image: `wwliao/hpp_sniffles:1.0.0--de39905b020c976bbebdd99d32c5df9b2c812ab4`

```sh
$ sniffles -s 4 \
           -l 40 \
           -n -1 \
           --cluster \
           --ccs_reads \
           -t 12 \
           -m hifi_alignments/HG002.GRCh38_no_alt.bam \
           -v calls/sniffles/outputs/HG002.GRCh38_no_alt.sniffles.raw.vcf \
  && echo HG002 > HG002.txt \
  && bcftools reheader -s HG002.txt \
                       -o calls/sniffles/outputs/HG002.GRCh38_no_alt.sniffles.unrefined.vcf \
                       calls/sniffles/outputs/HG002.GRCh38_no_alt.sniffles.raw.vcf \
  && rm HG002.txt calls/sniffles/outputs/HG002.GRCh38_no_alt.sniffles.raw.vcf \
  && iris genome_in=GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
          vcf_in=calls/sniffles/outputs/HG002.GRCh38_no_alt.sniffles.unrefined.vcf \
          reads_in=hifi_alignments/HG002.GRCh38_no_alt.bam \
          vcf_out=calls/sniffles/HG002.GRCh38_no_alt.sniffles.vcf \
          threads=12 \
          out_dir=calls/sniffles/tmp/HG002 \
          --also_deletions \
          --hifi \
	  --rerunracon \
          --keep_long_variants \
  && bcftools sort -m4G \
                   -Oz \
                   -o calls/sniffles/HG002.GRCh38_no_alt.sniffles.vcf.gz \
                   calls/sniffles/HG002.GRCh38_no_alt.sniffles.vcf \
  && bcftools index -t calls/sniffles/HG002.GRCh38_no_alt.sniffles.vcf.gz \
  && rm -rf calls/sniffles/tmp/HG002 \
  && rm calls/sniffles/HG002.GRCh38_no_alt.sniffles.vcf
```

## Run SVIM

Docker image: `wwliao/hpp_svim:1.0.0--93e7668e2ecd1fc315126e69da862fd365339509`

```sh
$ svim alignment --interspersed_duplications_as_insertions \
                 --read_names \
                 --zmws \
                 --cluster_max_distance 0.5 \
                 --minimum_depth 4 \
                 --min_sv_size 40 \
                 --sample HG002 \
                 calls/svim/outputs/HG002 \
                 hifi_alignments/HG002.GRCh38_no_alt.bam \
                 GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
  && bcftools sort -Oz \
                   -o calls/svim/HG002.GRCh38_no_alt.svim.vcf.gz \
                   calls/svim/outputs/HG002/variants.vcf \
  && bcftools index -t calls/svim/HG002.GRCh38_no_alt.svim.vcf.gz
```

## Run SVIM-asm

Docker image: `wwliao/hpp_svimasm:1.0.0--896384d810154ae24154f2143fe327e9e57b326a`

```sh
$ svim-asm diploid --interspersed_duplications_as_insertions \
                   --query_names \
                   --min_sv_size 40 \
                   --sample HG002 \
                   calls/svim-asm/outputs/HG002 \
                   asm_alignments/HG002.paternal.GRCh38_no_alt.bam \
                   asm_alignments/HG002.maternal.GRCh38_no_alt.bam \
                   GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
  && bcftools sort -Oz \
                   -o calls/svim-asm/HG002.GRCh38_no_alt.svim-asm.vcf.gz \
                   calls/svim-asm/outputs/HG002/variants.vcf \
  && bcftools index -t calls/svim-asm/HG002.GRCh38_no_alt.svim-asm.vcf.gz
```

## Run Truvari

Docker image: `wwliao/hpp_truvari:1.0.0--9352c3b9ffbe4680b1712ee465acdaec7f6f909f`

```sh
$ truvari bench -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
                -b HG002_SVs_Tier1_v0.6.GRCh38.vcf.gz \
                -c HG002.GRCh38_no_alt.${CALLER}.filtered.vcf.gz \
                --includebed HG002_SVs_Tier1_v0.6.GRCh38.sorted.bed \
                -o bench-${CALLER} \
                -r 1000 \
                -p 0.01 \
                --passonly \
                --no-ref a \
                --giabreport
```
