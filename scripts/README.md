Scripts and other post-analysis for TRUST4
=======

This folder contains scripts and comments for post processing TRUST4's results. The detailed information for each python scripts can be found by running "python3 the_script.py -h". The detailed information for perl scripts can be found by running "perl the_script.pl".

#### trust-stats.py

Compute the clonotype diversity in a sample

#### trust-cluster.py
Cluster similar CDR3s based on the trust_cdr3.out or trust_report.tsv file.

#### barcoderep-filter.py
Filter the lowly expressed clonotype if it is identical to another highly expressed clonotype in another cell. This strategy is inspired by 10x VDJ pipeline to remove contaminations from diffused mRNAs: https://support.10xgenomics.com/single-cell-vdj/software/pipelines/latest/algorithms/cell-calling.

#### barcoderep-expand.py
Put the secondary chains in the barcode_report.tsv file in its own barcode line. This function is useful if a barcode contains multiple cells. trust-airr.pl script can then be used to create AIRR entries for the secondary chains.

#### trust-barcoderep-to-10X.pl
Convert the barcode_report.tsv file to 10x CellRanger vdj format (the format used at least in CellRanger3).

#### airr-imgtgap.py
Add the gaps defined in the IMGT file to sequence_alignment and germline_alignment fields in the AIRR output.

#### epitope annotation
TRUST4's report file and the AIRR output format is compatible with eptiope prediction methods, such as [TCRMatch](https://github.com/IEDB/TCRMatch). You can use command like

`./tcrmatch -i trust_report.tsv -d CEDAR_data.tsv -t 8 -r > trust_with_epitope.txt`

or

`./tcrmatch -i trust_airr_report.tsv -d CEDAR_data.tsv -t 8 -a > trust_with_epitope.txt`

You can also upload the trust_report.tsv to TCRMatch's web portal on IEDB: [http://tools.iedb.org/tcrmatch/](http://tools.iedb.org/tcrmatch/).  
