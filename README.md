<img src="https://github.com/sunbeam-labs/sunbeam/blob/stable/docs/images/sunbeam_logo.gif" width=120, height=120 align="left" />

# sbx_cenote_taker

<!-- badges: start -->
[![Tests](https://github.com/sunbeam-labs/sbx_cenote_taker/actions/workflows/tests.yml/badge.svg)](https://github.com/sunbeam-labs/sbx_cenote_taker/actions/workflows/tests.yml)
[![DockerHub](https://img.shields.io/docker/pulls/sunbeamlabs/sbx_cenote_taker)](https://hub.docker.com/repository/docker/sunbeamlabs/sbx_cenote_taker/)
<!-- badges: end -->

## Introduction

sbx_cenote_taker is a [sunbeam](https://github.com/sunbeam-labs/sunbeam) extension for identifying viruses in samples with [Cenote-Taker3](https://github.com/jedvachey/Cenote-Taker3). This pipeline uses [MEGAHIT](https://github.com/voutcn/megahit) for assembly of contigs and then processes assemblies with Cenote-Taker3.

N.B. This extension requires also having sbx_assembly installed.

### Installation

```
sunbeam extend https://github.com/sunbeam-labs/sbx_assembly.git
sunbeam extend https://github.com/sunbeam-labs/sbx_cenote_taker.git
```

### Cenote-Taker database

sbx_cenote_taker expects the Cenote-Taker3 reference database to be available locally. Download the database following the official instructions, for example:

```
conda activate cenote-taker
get_ct3_dbs -o /path/to/ct3_db --hmm T --hallmark_tax T --refseq_tax T --mmseqs_cdd T --domain_list T --hhCDD T --hhPFAM T --hhPDB T
```

Update the `cenote_taker_db` entry in your Sunbeam configuration to point at the resulting directory.

### Running

Run with sunbeam on the target `all_cenote_taker`:

```
sunbeam run --profile /path/to/project/ all_cenote_taker
```

### Options for config.yml

  - blast_db: path to blast db (default: "") (NOTE: this should be the database file not just the directory it's in)
  - blastx_threads: number of threads for running blastx (default: 4)
  - bowtie2_build_threads: number of threads for running bowtie2-build (default: 4)
  - cenote_taker_db: path to cenote-taker3 db (default: "") (NOTE: this should be a directory)
  - include_phages: Whether to include phages in the output (default: False)
