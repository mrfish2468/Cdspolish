# Introduction
CDSpolish is a genome polisher developed for generating a high-quality genome for Nanopore bacteria sequence. Nanopore systematic errors are corrected by retreiving coding sequences from ENA database and polished by an Random Forest. When paired with Racon, Medaka and homopolish, the genome quality can slightly increase and frameshift in plasmids can be removed. For metagenome, indels can be reduced up to 34.7%.

# Installation
CDSpolish is recommendated to install and run within a conda environment

	git clone https://github.com/ythuang0522/homopolish.git
	cd cdspolish
	conda env create -f environment.yml
	conda activate cdspolish

# Quick usage

Cdspolish should be run with a pre-trained model (R9.4.pkl/R10.3.pkl for Nanopore and pb.pkl for PacBio CLR) and one sketch (virus, bacteria, or fungi). For Nanopore sequencing, Cdspolish should be run after the Racon-Medaka-Homopolish pipeline as it only removes indel errors. Note that the taxonomy of the draft genome need to be identified to genus level in advance for Cdspolish to work. For instance, if your Homopolish-polished genome (yourgenome.fasta) is bacteria and sequenced by R9.4 flowcell and genus is known(yourgenus), please type
```
python3 cdspolish.py yourgenome.fasta R9.4.pkl yourgenus youroutput
```
