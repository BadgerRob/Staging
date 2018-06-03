# Basecalling, demultiplexing and adaptor trimming

## Base calling using albacore

All fast5 files [pass + skip + QC] should be copied into a 'pass' directory. Basecall reads from fast5 files.

```
>mkdir pass
>cp path/to/skip -t path/to/pass -r
```

`-t` -output path  
`-r` -recursive



Call albacore `read_fast5_basecaller.py`

```
>read_fast5_basecaller.py -i path/to/fast5_dir/pass -t n -s path/to/basecalled_fast5/output -f FLO-MIN106 -k SQK-LSK108 --barcoding -o fast5 -r
```

`-i` [input file]  
`-t` [threads]  
`-s` [output/path]  
`-f` [flow_cell_chemistry] -9.4.1 (FLO-MIN106)  
`-k` [sequencing_kit] -ligation sequencing kit  
`-o` [output_file_type] -fastq  
`-r` [recursive]`  
`--barcoding` -demultiplex flag  


Basecalled and demultiplexed fastq reads saved in `output/workspace/`


## Adapter trimming and demultiplexing using porechop

Porechop can be used to remove sequencing adapters and barcode from fastq reads. 
Chimeric reads can be discarded. Barcoding binning can be set to require one or both barcodes per read. 
Barcode discrepancies between albacore and porechop can be discarded.  

More info at: https://github.com/rrwick/Porechop#quick-usage-examples

This code sets stringent binning with high barcode threshold requiring end to end sequences.

Call porechop `porechop`

```
>porechop -i path/to/albacore/basecalled_fastq/barcodes/ -b path/to/porechop/output --discard_middle --barcode_threshold 85 --require_two_barcodes
```

`-i` [input_dir] -from ablacore  
`-b` [demultiplexing_flag] -followed by `output/dir`  
`--discard_middle` -remove middle adaptors  
`--barcode_threshold 85` 
-stringent barcoding threshold from default  
`--require_two_barcodes` -require two barcodes from each read

Reads trimmed and demultiplexed for onward analysis in `/porechop/output/`

