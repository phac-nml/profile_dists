[![PyPI](https://img.shields.io/badge/Install%20with-PyPI-blue)](https://pypi.org/project/profile_dists/#description)
[![Bioconda](https://img.shields.io/badge/Install%20with-bioconda-green)](https://anaconda.org/bioconda/profile_dists)
[![Conda](https://img.shields.io/conda/dn/bioconda/profile_dists?color=green)](https://anaconda.org/bioconda/profile_dists)
[![License: Apache-2.0](https://img.shields.io/github/license/phac-nml/profile_dists)](https://www.apache.org/licenses/LICENSE-2.0)

# Profile Dists
![alt text](https://github.com/phac-nml/profile_dists/blob/main/logo.png?raw=true)

- [Introduction](#introduction)
  * [Citation](#citation)
  * [Contact](#contact)
- [Install](#install)
    + [Compatibility](#compatibility)
- [Getting Started](#getting-started)
  * [Usage](#usage)
  * [Configuration and Settings](#configuration-and-settings)
  * [Data Input](#data-input)
  * [Data Output](#data-output)
- [Troubleshooting and FAQs](#troubleshooting-and-faqs)
- [Other information](#other-information)
- [Legal and Compliance Information](#legal-and-compliance-information)
- [Updates and Release Notes](#updates-and-release-notes)

<small><i><a href='http://ecotrust-canada.github.io/markdown-toc/'>Table of contents generated with markdown-toc</a></i></small>

# Introduction

Surveillance and outbreak analysis of pathogens has been operationalized by multiple public health laboratories around the world using gene-by-gene approaches. Gene by gene comparisons can be considered the next iteration of multi-locus seuence typing (MLST) due to the compression of large amounts of sequence information into a set of integer allele identifiers according to a defined schema. There are a variety of software approaches available to call allele profiles from schemes such as chewBBACA. Additionally there are tools available which can calculate a distant matrix using the standardized allele profile format of ['name', 'locus_1','locus_2'...'locus_n'] in the form of GrapeTree. However, this method is not tolerant of non-integer data in loci columns and one of the strengths of chewBBACA is that it will provide additional information in columns in the case of partial matches etc. Furthermore, there is interest in utilizing hashing of allele sequences for deterministic labelling of alleles without the need for a centralized allele nomenclature service. Additionally, GrapeTree does not provide functionality for querying samples against a database of samples. Furthermore, the memory requirements of GrapeTree grow significantly with the number of samples being compared as the distance matrix is held in memory. As datasets grow, using text based formats such as CSV, TSV represent significant amounts of runtime in terms of reading, parsing and writing. New formats such as parquet support compression and are optimized for efficiency of storage and retrieval of data.

To address needs of users within our team we have designed an integrated solution for calculating distance matricies and querying of genetically similar samples within a defined threshold to support outbreak and surveillance activities. We provide the flexibility to have standard text based outputs as well as parquet. It is implemented in pure python and currently is only available in a single threaded version but later refinements may include the support for multiprocessing. To facilitate integration of the tool into larger workflows it will also be implemented as a nextflow workflow.
Citation
## Citation

Robertson, James, Wells, Matthew, Schonfeld, Justin, Reimer, Aleisha. Profile Dists: Convenient package for comparing genetic similarity of samples based on allelic profiles. 2023. https://github.com/phac-nml/profile_dists
## Contact

[James Robertson] : <james.robertson@phac-aspc.gc.ca>

# Install

Install the latest released version from conda:

    conda create -c bioconda -c conda-forge -n profile_dists profile_dists

Install using pip:

    pip install profile_dists

Install the latest master branch version directly from Github:

    pip install git+https://github.com/phac-nml/profile_dists.git
    
### Compatibility

All required packages are listed in the setup.py. The package does require Numba to JIT compile certain functions, if any issues during installation are encountered please make an issue on the repository.

# Getting Started

## Usage

### **Distance matrix calculation**

The default behaviour of profile dists is to construct a pairwise distance matrix with the columns representing the
reference sequences and rows representing the query sequences. For clustering analysis, you can construct a symmertic/square
distance matrix by supplying the same profile as reference and query.

```
profile_dists --query samples.profile --ref samples.profile --outdir results
```

### Fast Matching

Profile dists also supports querying of a set of sample profiles against a reference database and reporting that as either a
matrix or three column table of [ query_id, reference_id, dist ]

```
profile_dists --query queries.profile --ref reference.profile --outdir results --match_threshold 10
```

## Configuration and Settings

If you run ``profile_dists --help``, you should see the following usage statement describing each parameter `profile_dists` accepts:

        Profile Dists: Calculate genetic distances based on allele profiles v. 1.0.0

        options:
                -h, --help            show this help message and exit
                --query QUERY, -q QUERY
                                        Query allelic profiles (default: None)
                --ref REF, -r REF     Reference allelic profiles (default: None)
                --outdir OUTDIR, -o OUTDIR
                                        Result output files (default: None)
                --outfmt OUTFMT, -u OUTFMT
                                        Out format [matrix, pairwise] (default: matrix)
                --file_type FILE_TYPE, -e FILE_TYPE
                                        Out format [text, parquet] (default: text)
                --distm DISTM, -d DISTM
                                        Distance method raw hamming or scaled difference [hamming, scaled] (default: scaled)
                --missing_thresh MISSING_THRESH, -t MISSING_THRESH
                                        Maximum percentage of missing data allowed per locus (0 - 1) (default: 1.0)
                --sample_qual_thresh SAMPLE_QUAL_THRESH, -c SAMPLE_QUAL_THRESH
                                        Maximum percentage of missing data allowed per sample (0 - 1) (default: 1.0)
                --match_threshold MATCH_THRESHOLD, -a MATCH_THRESHOLD
                                        Either a integer or float depending on what distance method is used (only used with pairwise format (default: -1)
                --mapping_file MAPPING_FILE, -m MAPPING_FILE
                                        json formatted allele mapping (default: None)
                --batch_size BATCH_SIZE
                                        Manual selection of how many records should be included in a batch, default=auto (default: None)
                --max_mem MAX_MEM     Maximum amount of memory to use (default: None)
                --force, -f           Overwrite existing directory (default: False)
                -s, --skip            Skip QA/QC steps (default: False)
                --columns COLUMNS     Single column file with one column name per line or list of columns comma separate (default: None)
                -n, --count_missing   Count missing as differences (default: False)
                -p CPUS, --cpus CPUS  Count missing as differences (default: 1)
                -V, --version         show program's version number and exit 

## Data Input

Profile dists supports several different formats as input: 

**Native**

|  id  |  locus_1  |  locus_2  |  locus_3  |  locus_4  |  locus_5  |  locus_6  |  locus_7  | 
| ----------- | ----------- |----------- | ----------- | ----------- |----------- | ----------- | ----------- |
|  S1  |    1  |  1  |  1  |  1  |  1  |  1  |  1  | 
|  S2  |    1  |  1  |  2  |  2  |  ?  |  4  |  1  | 
|  S3  |    1  |  2  |  2  |  2  |  1  |  5  |  1  | 
|  S4  |    1  |  2  |  3  |  2  |  1  |  6  |  1  | 
|  S5  |    1  |  2  |  ?  |  2  |  1  |  8  |  1  | 
|  S6  |    2  |  3 |  3  |  -  |  ?  |  9  |  0  | 

- Direct support for missing data in the form of ?, 0, -, None or space


**chewBBACA**

|  id  |  locus_1  |  locus_2  |  locus_3  |  locus_4  |  locus_5  |  locus_6  |  locus_7  | 
| ----------- | ----------- |----------- | ----------- | ----------- |----------- | ----------- | ----------- |
|  S1  |    1  |  INF-2  |  1  |  1  |  1  |  1  |  1  | 
|  S2  |    1  |  1  |  2  |  2  |  NIPH  |  4  |  1  | 
|  S3  |    1  |  2  |  2  |  2  |  1  |  5  |  1  | 
|  S4  |    1  |  LNF  |  3  |  2  |  1  |  6  |  1  | 
|  S5  |    1  |  2  |  ASM  |  2  |  1  |  8  |  1  | 
|  S6  |    2  |  INF-8  |  3  |  PLOT3  |  PLOT5  |  9  |  NIPH  | 

- All non integer fields will be converted into missing data '0'

**Hashes**

|  id  |  locus_1  |  locus_2  |  locus_3  |  locus_4  |  locus_5  |  locus_6  |  locus_7  | 
| ----------- | ----------- |----------- | ----------- | ----------- |----------- | ----------- | ----------- |
|  S1  |    dc0a04880d1ad381ffd54ce9f6ad1e7a |  -  |  b9f94bf167f34b9fcf45d79cab0e750a  |  8a07b9cb0ab7560ad07b817ca34036bb  |  80c8156d77d724ac0bb16ec60993bc84  |  7a1d0a48f16fa25910cddfea38dab229  |  e1ee776b32c2f6131a7238ce50b75469  | 
|  S2  |    dc0a04880d1ad381ffd54ce9f6ad1e7a  |  0af06522a32865cd2db2cf5a854d195b  |  9fc502308c616ae34146d7f7b0081bd8  |  4577dec2c840472800a3b104c88bb0ef  |  -  |  bba24c25c28c08058d6f32ecfbf509e9  |  e1ee776b32c2f6131a7238ce50b75469  | 
|  S3  |    dc0a04880d1ad381ffd54ce9f6ad1e7a  |  04d45219ee5f6065caf426ba740215e5  |  9fc502308c616ae34146d7f7b0081bd8  |  4577dec2c840472800a3b104c88bb0ef  |  80c8156d77d724ac0bb16ec60993bc84  |  874225c0dec5219dd64584ba32938dbd  |  e1ee776b32c2f6131a7238ce50b75469  | 
|  S4  |    dc0a04880d1ad381ffd54ce9f6ad1e7a  |  -  |  e79562c280691c321612ecdf0dadad9e  |  4577dec2c840472800a3b104c88bb0ef  |  80c8156d77d724ac0bb16ec60993bc84  |  c8087ad8b01d9f88e8eb2c3775ef2e64  |  e1ee776b32c2f6131a7238ce50b75469  | 
|  S5  |    dc0a04880d1ad381ffd54ce9f6ad1e7a  |  04d45219ee5f6065caf426ba740215e5  |  -  |  4577dec2c840472800a3b104c88bb0ef  |  80c8156d77d724ac0bb16ec60993bc84  |  4d547ea59e90173e8385005e706aae96  | e1ee776b32c2f6131a7238ce50b75469  | 
|  S6  |    8214e9d02d1b11e6239d6a55d4acd993  |  -  |  e79562c280691c321612ecdf0dadad9e  |  -  |  -  |  e3088425be5e7de8d9a95da8e59a9ea8  |  -  | 

- Direct support for missing data in the form of ?, 0, - or space

## Data Output

### Output Profile

**Native**

|  id  |  locus_1  |  locus_2  |  locus_3  |  locus_4  |  locus_5  |  locus_6  |  locus_7  | 
| ----------- | ----------- |----------- | ----------- | ----------- |----------- | ----------- | ----------- |
|  S1  |    1  |  1  |  1  |  1  |  1  |  1  |  1  | 
|  S2  |    1  |  1  |  2  |  2  |  0  |  4  |  1  | 
|  S3  |    1  |  2  |  2  |  2  |  1  |  5  |  1  | 
|  S4  |    1  |  2  |  3  |  2  |  1  |  6  |  1  | 
|  S5  |    1  |  2  |  0  |  2  |  1  |  8  |  1  | 
|  S6  |    2  |  3 |  3  |  0  |  0  |  9  |  0  | 

- All columns are converted to contain only integers with missing data represented as a 0

### Output Directory Structure

```
{Output folder name}
├── allele_map.json - Mapping of allele hash string to integer id, can be use for reruning of analyses or masking specific alleles
├── query_profile.{text|parquet}  - Standardized allele profile for query sequences
├── ref_profile.{text|parquet}  - Standardized allele profile for reference sequences
├── results.{text|parquet} - Either symmetric distance matrix or three column file of [query_id, ref_if, distance]
└── run.json - Contains logging information for the run including parameters and quality information
```

# Troubleshooting and FAQs

No issues are currently reported.

# Other information

# Legal and Compliance Information

Copyright Government of Canada 2023

Written by: National Microbiology Laboratory, Public Health Agency of Canada

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this work except in compliance with the License. You may obtain a copy of the License at:

[http://www.apache.org/licenses/LICENSE-2.0](http://www.apache.org/licenses/LICENSE-2.0)

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.

# Updates and Release Notes

v1.0.0 

Please see the `CHANGELOG.md`.
