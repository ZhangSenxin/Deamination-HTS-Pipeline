# README

***Note:***

_FILE: a file path

_PATH: a folder path

<br/>

# Barcode prepare

**I. Cluster inbarcode.**

**Timing:** ~19min for 12000 sequences

The script `cluster_inbarcode.py` extracts a 9-nt sequence on both sides to perform clustering method allowing less than 2 nt mismatch. The output from this step will be used in the following barcode assignment step. The sequence file (specified with the `-i` or `--input_path` option) and the output path (specified with the `-o` or `--output_path` option) are required. Once the program finishes, it will print the number indicating the largest cluster size, and output a `.pkl` file containing the cluster information and the input data for the next step.

```shell
python3 cluster_inbarcode.py -i SEQ_FILE -o CLUSTER_PATH

python3 cluster_inbarcode.py --input_path SEQ_FILE --output_path CLUSTER_PATH
```

<br/>

**II. Assign barcode to each sequenceh.**

**Timing:** <1min for 12000 sequences

To determine a two-strand barcode list based on the maximum number printed at the previous step, run `barcode_assignment.py` to assign a barcode to each sequence. Sequences in the same cluster will be assigned different barcodes. The cluster `.pkl` file’s folder path (specified with the `-i` or `--input_path` option), the output path (specified with the `-o` or `--output_path` option), and the barcode file (specified with the `-m` or `--meta_file` option) are all required. 

```shell
python3 barcode_assignment.py -i CLUSTER_PATH -o SAVE_FILE -m BARCODE_FILE

python3 barcode_assignment.py --input_path CLUSTER_PATH --output_path SAVE_FILE --meta_file BARCODE_FILE
```

<br/>

**III. Perform simulation test to demultiplex each sequence according to the in-barcodes.**

**Timing:** ~1min for simulating 12 million data (-r 1000)

To simulate sequencing data, run the `fastq_simulation.py` script. This requires specifying the input sequence file with barcode as the output of previous step (using the `-i` or `--input_path` option) and the output path (using the `-o` or `--output_path` option). Optionally, you can use `-r` or `--repeat` to specify the simulated sequencing depth (default is `-r 10`).

```shell
python3 fastq_simulation.py -i SEQ_FILE -o FASTQ_PATH

python3 fastq_simulation.py --input_path SEQ_FILE --output_path FASTQ_PATH
```

<br/><br/>

# Deamination-HTS-Pipeline

**I. Multx reads according to inbarcode**

**Timing:** ~60 min for ~30 million reads

In this step, we split the meta into 24 sub-metas, each containing a maximum of 500 sequences. This is because the program `fastq-multx` limits the number of sequences in each run. Then we run `fastq-multx` to split each read into different files according to barcode and inbarcode. Note that reads in each file are considered from the same reference sequence.

```shell
fastq-multx -x -B META_FILE -m 1 -d 1 -b R1.gz R2.gz -o OUT_R1.gz OUT_R2.gz
```

<br/>

**II. Annotation and generating mutation profile files**

**Timing:** ~75 min for ~30 million reads

To perform sequence annotation and mutation computation, we use `read_fq.gz.py`. The annotation step requires a reference file (shown as /meta/meta_12k_ref.txt, input as `-m` or `--meta_file`) to calculate mutation information. We specify the data path (input as `-f` or `--file_path`) and save path (input as `-s` or `--save_path`) to generate mutation profile files. Optionally, we can use `-c` or `--chain` to specify the top strand (`-c 1`, default) or bottom strand (`-c 0`).

```shell
python3 read_fq.gz.py -m META_FILE -f SEQ_PATH -s SAVE_PATH -c 1

python3 read_fq.gz.py --meta_file META_FILE --file_path SEQ_PATH --save_path SAVE_PATH --chain 1
```

<br/>

**III. Basic mutation statistic**

Timing: ~10 min for ~30 million reads

To perform mutation statistics at base resolution (A, G, T, C, WRC, AGC, AGCT, WGCW, NWRC), run the program `mut_summary.py`. This program calculates the sum of mutations for all sequences and positions.

You need a dedicated region reference file (shown as /meta/IGHV_ref_kabat.txt) to specify region information for each read. Input this file as `-m` or `--meta_file`. Specify the data path (input as `-f` or `--file_path`) and save path (input as `-s` or `--save_path`) to generate a basic mutation statistic file.

```shell
python3 mut_summary.py -m META_FILE -f SEQ_PATH -s SAVE_PATH

python3 mut_summary.py --meta_file META_FILE --file_path SEQ_PATH --save_path SAVE_PATH
```

<br/>

**IV. Motif mutation statistic**

Timing: ~15 min for ~30 million reads

In contrast to the previous step, we use an additional program called `mut_summary_2.py` to gather mutation statistics. This program records mutation information from each target position of each reference sequence. The parameters required are identical to those of the previous step.

```shell
python3 mut_summary_2.py -m META_FILE -f SEQ_PATH -s SAVE_PATH

python3 mut_summary_2.py --meta_file META_FILE --file_path SEQ_PATH --save_path SAVE_PATH
```

<br/>

**V. Mapping to full-length reference**

Timing: <1 min for ~30 million reads

To perform further analysis, we use `mapping2full_length_sequence.py`. This tool combines each 50-nt segment from the same V gene to obtain complete mutation information. To recover the full sequence, you will need the V gene reference (shown as /meta/V_region_kabat.csv, input as `-m` or `--reference_path`) and the region reference file (input as `-r` or `--region_path`). By specifying the data path (input as `-f` or `--file_path`), save path (input as `-s` or `--save_path`), and another mismatch path (input as `-i` or `--mismatch_path`), all the full length information will be generated in a specialized data format.

```shell
python3 mapping2full_length_sequence.py -m META_FILE -r REGION_FILE -f SEQ_PATH -s SAVE_PATH -i MISMATCH_FILE

python3 mapping2full_length_sequence.py --reference_path META_FILE --region_path REGION_FILE --file_path SEQ_PATH --save_path SAVE_PATH --mismatch_path MISMATCH_FILE
```

<br/>

**VI. Prepare visualization**

Timing: <1 min for ~30 million reads

To prepare data for visualization, we need to run `data_produce.py`. This script summarizes the output from `mapping2full_length_sequence.py` into several matrix data. To run `data_produce.py`, we need the following meta files: V gene reference (input as `-m` or `--reference_path`), region reference file (input as `-r` or `--region_path`), and species reference (input as `-p` or `--species_path`). Once we specify the full-length output data path (input as `-f` or `--file_path`) and save path (input as `-s` or `--save_path`), many motif mutation data for visualization will be generated.

```shell
python3 data_produce.py -m META_FILE -r REGION_FILE -f SEQ_PATH -s SAVE_PATH -p SPECIES_FILE

python3 data_produce.py --reference_path META_FILE --region_path REGION_FILE --file_path SEQ_PATH --save_path SAVE_PATH --species_path SPECIES_FILE
```

<br/>

**VII. Plot figures**

Timing: <1 min for ~30 million reads

Once all necessary preparations have been made, run `figure_plot.py` in order to generate the box plot. Note that a species reference meta file (input as `-p` or `--species_path`) is required. To complete the visualization work, specify the output path of `data_produce.py` as the data path (input as `-f` or `--file_path`) and the save path (input as `-s` or `--save_path`).

```shell
python3 figure_plot.py -f SEQ_PATH -s SAVE_PATH -p SPECIES_FILE

python3 figure_plot.py --file_path SEQ_PATH --save_path SAVE_PATH --species_path SPECIES_FILE
```
