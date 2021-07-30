## Tandem repeat annotation for GRCh38 no alt

```sh
$ wget -O GRCh38_no_alt.trf.bed \
       https://raw.githubusercontent.com/PacificBiosciences/pbsv/master/annotations/human_GRCh38_no_alt_analysis_set.trf.bed
```

## Tandem repeat annotation for CHM13 v1.1 + chrY + EBV

```sh
$ wget http://t2t.gi.ucsc.edu/chm13/hub/t2t-chm13-v1.1/trf/trf.bigBed
$ wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed
$ chmod +x bigBedToBed
$ ./bigBedToBed trf.bigBed CHM13.trf.bed
$ cat <(cut -f1-3 CHM13.trf.bed) \
      <(grep -w "^chrY" GRCh38_no_alt.trf.bed) \
      <(grep -w "^chrEBV" GRCh38_no_alt.trf.bed) \
      > CHM13Y_EBV.trf.bed
```
