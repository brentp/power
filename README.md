## mospow

`mospow` is a tool for quickly determining your power to detect mutation in WGS alignments

specifically, `mospow` calculates the fraction of the genome covered by at least X reads in a single individual, across a trio of samples, or across a tumor-normal pair

```
mospow [options] BED \
    -d --depth_cutoff           depth cutoff (default = 10)
    -b --bed_regions            BED file of regions to restrict power calculation 
    -r --region                 instead of -b, single region <chr:start-end>
```

## what

given an exome, targetted sequencing, or WGS we want to know, simply, how well were (the) regions captured?

but that's not so simple.

the questions to ask/answer are (with '10X' a place-holder for a generic coverage cutoff)

1. what fraction of bases where covered at 10X? (this is output by mosdepth)
2. what fraction of regions were at 10X?
3. across $n samples (e.g. a trio) what fraction of bases/regions were at 10X in all samples
4. across a cohort, which samples appear to be problematic--e.g. were missing too much data
   in critical regions?
5. given a set of important variants, (e.g. clinvar), how well are the exons (regions) that harbor
   those variants covered?

Especially for 5, we want a single number that we can use in a cohort to see which samples
are poorly covered--e.g. exome samples that might need WGS or another lane of sequencing.

## output (planned)

a user should be able to quickly see which samples have critical regions below some threshold.

### Tables

#### High-level coverage stats.

given some regions of interest (ROI) we will have text and html outout of:

sample | median-coverage | % bases at 1X | % bases at 10X | % ROI at 10X 
------ | --------------- | ------------- | -------------- | ------------
mom    |  31             |  97.4         |    94.3        |  98.4
dad    |  34             |  98.1         |    95.3        |  98.9
kid    |  24             |  94.1         |    92.1        |  98.9
...    |  ..             |  ....         |    ....        |  ....

#### Per-gene (region) coverage stats.

for each region, the 4th column is assumed to be the gene. we will aggregate coverage across
the gene and report:


gene   |  median coverage       |  % ROI bases at 10X 
       |  mom | dad | kid | ... | mom | dad | kid | ... 
------ | ---- | --- | --- | --- | --- | --- | --- | ... 
TERT   |  31  |  21 | 22  | ... | 95  | 88  | 22  | ... 
TNFA   |  34  |  25 | 17  | ... | 87  | 72  | 22  | ...
KCNQ2  |  24  |  25 | 22  | ... | 99  | 99  | 99  | ...
...    |  ..  | ... | ..  | ... | ... | ..  | ..  | ...


in the HTML report cells with low coverage will be colored red by the severity.

this will also output JSON for per-base coverage of each gene in the ROI list. This could be a 
lot of data so it might limit to 

#### Per-gene coverage plots

In the table above, when a cell is clicked(?), a plot will show the per-base coverage ofeach
that gene and that sample in a plot where introns are minmized. additional samples can be added $somehow.

These plots will show the location of clinvar variants in that gene; they'll be plotly.js plots.

igv.js [supports tabix](https://github.com/igvteam/igv.js/blob/master/test/testTabix.js) so if mosdepth
output .tbi instead of .csi, this would be suported easily...

## Invocation

given the above requirements, invocation would look like:

```
  mospow \
      --variants $clinvar.vcf \
      -d 10 \
      -b regions_of_interest.bed \ # 4th col is gene name
      $sample.mosdepth-output.per-base.bed.gz
```



## previous work

[chanjo](https://www.chanjo.co/) does a lot of this using [sambamba](https://github.com/lomereiter/sambamba)

chanjo covers some of the same use-cases as `mospow` but `chanjo` is region-based and doesn't enable some of those
per-base stuff considered in `mospow`
