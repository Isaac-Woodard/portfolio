Profile Page Link:

### Overview

### Usage

Sort
<0/1(no/with Cache)> <degree> <gbk file> <sequence length> [<cache size>] [<debug level>]
DNAfilename.gbk.btree.data.sequencelength.degree

Search
<0/1(no/with Cache)\> <btree file\> <query file\> [<cache size\>] [<debug level\>]

Example
java GeneBankCreateBTree 0 2 DNAExample.gbk 5
java GeneBankSearch 0 DNAExample.gbk.btree.data.5.2 QueryExample.txt

### Expected Frequencies from DNAExample
2 gaagc (4 reported)
0 tcccc (won't be printed to console)
5 tagat (6 reported)
