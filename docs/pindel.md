# Pindel

Analysis steps involved:

![Steps Pindel](images/pindel.png)

Pindel is a tool for detecting larger structural variants, such as insertions and deletions, from short-read 
sequencing data.
Originally developed by [Ye et al. in 2009](https://watermark.silverchair.com/bioinformatics_25_21_2865.pdf?token=AQECAHi208BE49Ooan9kkhW_Ercy7Dm3ZL_9Cf3qfKAc485ysgAAA7UwggOxBgkqhkiG9w0BBwagggOiMIIDngIBADCCA5cGCSqGSIb3DQEHATAeBglghkgBZQMEAS4wEQQMhI2zgV-N3r0Ju31AAgEQgIIDaH1PGdGJcXfM8ZSks51MJtt6eIK-WTeIh0xvteS5tAv7w_qS5YJUXBDf3cxIvUM6x1JWVWYgd0dgMoap18UpGS7yUn_qR8ukkh8W60HNvBR3X5pjiPVH4iZMyAwFmI0MGu43zOw-sJ_i7LlGLFsq4M5KkU5lUC3Md41DAx2efDHbmEzXz9VPqScxfbNxhuuUJDHMkVS9lsPnWIt4-ZTvgdWARgFrG4IjnR0tGHJ-2ACOsYuTZ6m60dJSKEw_S6QNIP7UP-NFML0orkSWECpmbeuu9JLlo2TnJmqW74aqhb1Jngx2FRcJMcrvsNHEFd161syFqUS6Oeswj-A-0eRZdT588afyUeKz7oiRBIpwTXo7pSeSRN8HCstxUokISNFcvIdsvu87bymQNcYfL4YZhQjQkJowJatbT0JwwVmqprzlfbh3iTLJTZmo8D_bNQ4T_jFOhks2AtmbF9yS0BpEjtdcptIQr4914dCCRD00AiFGYcWSrJibfj7tmGpJjSfOYNSp5uQwp6uvs-CcD1SkOilDtZIL_hai0zw1i7206qipbbGTRbXTr9L0YlxL2mz8w2zCtOX1taHL-m5jOeVFOgZQwIo_uml7RJMejA6sAr1v_JDIoUuCErBaH3tUBn01r2sQimzP5zKn5rGkxmEQMpmWaPnx35XZ6jKb6Q5iyO3pIk4uF_RzZC0pH2WL04a5HkIRm2oW0qgHG__pTuwByJzprKYDmfvWFxn2ZhV__i90HMx07d9VvwtlfmdZceSm6jClBczqjiwXgD6txCELsith1ggmfPmnZewf80btn91tz0cY4YHOqmHo-uK9sIOZSJO6twBXJO3FkcJo0J6m-8rCSN4egov5PczPI9yZg-Mxx2X5vT6Vt4f3ER4bP-i5mcuOjkk9RZ_Za-dF8NN-0F7Glhczf5XurHsN_O11MK_ZRzBKgBQAho3j0SM0lVTAIO7kYsQD4zgKDA7Vz8tjvDj2BUWYL-WMtLEPam6A5AwrAQJatgFv4oCaFgbrP9m77r0Xt-8g12gt3qOLOsgvHx6yMCSQg9RgSDf6bEUuYxVEFMFvuccSofX5WvitGE2vbNvtJjnCantwQYV9VYxkV6RzqMXAz8o1mKu53Occ-3OApGtYacCuIWsNLO8O32pdE7zf9MUh7wRc), it is particularly effective for identifying complex indels that may 
not be captured by standard variant calling tools like GATK Mutect2 or VarDict. 
Pindel uses a split-read approach, where it
analyzes reads that are split across the boundaries of an indel to accurately determine its size and location.
For more details, see the documentation of [Poppy GMS]().

In the context of the poppy-uppsala pipeline, Pindel is used to call larger indels in specific regions of 
interest that involve the following genes:

- *SRP72*
- *NPM1*
- *FLT3*
- *RUNX1*
- *TP53*
- *NOTCH1*
- *WT1*
- *NF1*
- *UBTF*
- *PPM1D*
- *CALR*
- *ASXL1*
- *CEBPA*
- *GATA1*

Pindel is [known to run (very) slowly and to have a high memory usage](https://github.com/genome/pindel/blob/master/FAQ),
therefore the analysis is performed in a restricted set of regions in the Twist Myeloid panel to reduce the 
computational cost of the analysis.

After variant calling with Pindel, hard (1.) and soft (2.) filters are applied to the calls: 

1. Variant calls with VAF <= 1% or DP < 100, or AD < 5 are filtered out
2. Intron and splices variants as annotated by VEP, germline calls, artifacts as well as variants on the UBTF gene
are filtered out. That is, they are hidden by default upon opening of the Excel report.

Moreover, the AF field must be edited in all variants as it is set by Pindel as the "Allele count based on AD-field"
and not as a frequency.

Note that the artifacts are annotated w.r.t. the machines used at CGU (NextSeq and NovaSeqX), they may not be 
relevant if other sequencing machines or other setups are used.

## Implementation

Pindel is run with the default options in poppy-uppsala.
It is seen as a variant calling tool for CNVs and SVs in hydra-genetics, however the VCF output of Pindel can be 
processed with tools for indels such as VT for decomposition and normalization.

The output of Pindel requires additional processing to be usable in the context of the poppy-uppsala pipeline,
for instance:
- Reformat the output to regular VCF files,
- Compute the variant allele frequency from the allelic depth and 
  add the AF value to the INFO fields in the VCF file.

The Python scripts fixing these processing steps and the related Snakemake rules can be found in the 
`poppy/workflow/scripts/` and `poppy/workflow/rules/` directories of the Poppy GMS repository:
- `create_artifact_file_pindel.py` and `pindel_processing_*.py` scripts,
- `pindel_processing.smk` set of rules.

As the [source code of Pindel](https://github.com/genome/pindel/tree/master) is no longer maintained in its open-source repository, we would ideally like to
replace it with a more modern tool in the future.

### References

- Ye, K., Schulz, M. H., Long, Q., Apweiler, R., & Ning, Z. (2009). Pindel: a pattern growth approach to detect break points of large deletions and medium sized insertions from paired-end short reads. *Bioinformatics*, 25(21), 2865-2871. doi:10.1093/bioinformatics/btp394