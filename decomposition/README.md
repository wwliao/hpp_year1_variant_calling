# Graph variant decomposition

Work in progress. Please don't use the script.

1. Extract a given sample from the VG-deconstructed VCF:

    ```sh
    PREFIX=$(basename $VCF .vcf.gz)
    bcftools view -a -I -s $SAMPLE -Ou $VCF \
        | bcftools view -e 'GT="ref" | GT="0|." | GT=".|0" | GT=".|." | GT="."' -Ou \
        | bcftools norm -m - -Oz -o $SAMPLE.$PREFIX.vcf.gz \
	    && bcftools index -t $SAMPLE.$PREFIX.vcf.gz
    ```

