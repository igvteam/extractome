# extractome

Creates an "extracted" genome assembly from a  fasta + list of regions in "bed" format

## Usage

```python

extractome regions.bed <options>

or

python extractome/extract.py regions.bed <options>

```

**Required inputs**

* Region file in "bed" format, described [here](https://genome.ucsc.edu/FAQ/FAQformat#format1).  Only the first 3 columns are required.
* Either a fasta file or IGV genome identifer (see Options below)

**Options**

* --fasta reference fasta file, required if --genome is not specified
* --genome igv.js genome id (e.g. hg38), required if --fasta is not specified
* --name base name for output files, default=Xome
* --output output directory name, default=output


## Output

The script creates 3 output files

* base_name.fa
* base_name.regions.bed  - the input regions file lifted over to extracted fasta
* base_name.chain  - a UCSC "chain" file. Can be used to liftover files to the extracted fasta with tools such as [CrossMap](http://crossmap.sourceforge.net/)


## Example

```
extractome  test/data/cpgIsland_mm10.bed --genome mm10 --name CPG_mm10 --output output 

```

