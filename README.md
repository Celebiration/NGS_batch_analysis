

# About

This program analizes NGS data parallelly. By specifying a simple mapping table file it automatically does adapter truncation, reads merging, demultiplexing, editing outcome statistics tabulating and sample comparison, or it tabulates sequences occurrences for randomized DNA library.

# Usage

Originally the code is implemented in a web-server. However, it can serve as a command line program with a little modification. Following these steps to analyze your data:

## step1: prepare a mapping file

download the mapping.csv template and fill the table with your NGS library information.

## step2: prepare your data

create a directory named "fq" which contains all your NGS sequencing data. Make sure it matches your mapping.csv file(i.e. file names specified in mapping.csv correspond to files in fq/). 

## step3: run the program

Run the program using command:

```bash
python process.py [input directory] [output directory]
```

