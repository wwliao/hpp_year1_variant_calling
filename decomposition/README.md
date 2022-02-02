# Graph variant decomposition (work in progress)

1. Extract a given sample from the VG-deconstructed VCF

    ```sh
    PREFIX=$(basename $VCF .vcf.gz)
    bcftools view -a -I -s $SAMPLE -Ou $VCF \
        | bcftools view -e 'GT="ref" | GT="0|." | GT=".|0" | GT=".|." | GT="."' -Ou \
        | bcftools norm -m - -Oz -o $SAMPLE.$PREFIX.vcf.gz \
	    && bcftools index -t $SAMPLE.$PREFIX.vcf.gz
    ```

2. Convert the format of HPRC pangenome graphs from GFA to HashGraph

    There are two HPRC pangenome graphs:
    
    - Minigraph-Cactus: [hprc-v1.0-mc-grch38.gfa.gz](https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.0-mc-grch38.gfa.gz)
    - PGGB: [hprc-v1.0-pggb.gfa.gz](https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/pggb/hprc-v1.0-pggb.gfa.gz)

    Using VG to convert GFA to HashGraph

    ```sh
    vg convert -g -a -t 36 <input GFA> > <output HashGraph>
    ```

3. Decompose graph variants
    
    ```sh
    python3 decompose_graph_variants.py -o $SAMPLE.$GRAPH.decomposed.vcf.gz \
                                           $GRAPH.hg \
                                           $SAMPLE.$GRAPH.vcf.gz \
        && bcftools sort -m 10G \
                         -T $SAMPLE_sort_tmp/ \
                         -Oz -o $SAMPLE.$GRAPH.decomposed.sorted.vcf.gz \
                         $SAMPLE.$GRAPH.decomposed.vcf.gz \
        && bcftools index -t $SAMPLE.$GRAPH.decomposed.sorted.vcf.gz
    ```

