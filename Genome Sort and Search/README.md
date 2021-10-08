Profile Page Link: https://sites.google.com/view/isaacwoodard/projects/genome-sort-search

# Overview
The final project of my CS 321 Data Structures course required sorting the human genome using a balanced-tree (B-tree) data structure. The goal was to then search for the frequency of particular DNA sequences in the genome. I worked on the project over the course of a month with my partner Aaron Parish.

### Process
The first step in  the process is creating the B-tree file with GeneBankCreateBTree.java. The steps to make the B-tree file can be summarized as:
 - Parse segment of Genebank file
 - Convert into binary
 - Sort through B-tree
 - Write to disk

This is done in segments due to the size of the human genome. Each DNA subsequence of a user specified length is stored in a binary B-tree file. A counter for each subsequence is also stored to keep track of how many times it was found.

After creating the B-tree file, a search can be made using GeneBankSearch.java. The DNA sequences to search for are specified using a query file. The sequences are read in, converted to binary, and then searched for in the B-tree file. If found, the frequency of the sequence is printed to console. 

The counting doesn't quite work, as can be seen in the example below.

### Files
**Main Classes:**
 - **GeneBankCreateBTree.java:** Sorts a genome file and stores as a binary B-tree file.
 - **GeneBankSearch.java:** Searches for a DNA sequence in a query file.

**B-Tree Classes:**
 - **IUDoubleLinkedList.java:** Manages a double-linked list of objects held in nodes.
 - **LinearNode2.java:** A node holding an object and references to other nodes in a double-linked list.
 - **IndexedUnsortedList.java:** Interface which defines the methods for IUDoubleLinkedList.java.

**List Classes (From CS 221):**
 - **BTree.java:** Manages the B-Tree file. Has inner classes BTreeNode.java and Cache.java.
 - **TreeObject.java:** Handles keys for B-tree nodes.

**Other Files**
 - **Genebank:** Genome from NCBI.
 - **B-tree:** Binary file containing sorted genome (created by GeneBankCreateBTree.java).
 - **Query:** Text document containing DNA sequence to search for.

# Usage
First GeneBankCreateBtree.java must be run to count DNA sequence frequencies and create the B-tree file; then GeneBankSearch.java can be run to search for dna sequences specified in a query file.

### GeneBankCreateBTree

Command line args: <0/1(no/with Cache)\> <degree\> <gbk file\> <sequence length\> [<cache size\>] [<debug level\>]
 - Degree should be at least 2.

Creates a binary file in the working directory with the naming convention: DNAfilename.gbk.btree.data.sequencelength.degree

### GeneBankSearch

Command line args: <0/1(no/with Cache)\> <btree file\> <query file\> [<cache size\>] [<debug level\>]

### Example

java GeneBankCreateBTree 0 2 DNAExample.gbk 5

java GeneBankSearch 0 DNAExample.gbk.btree.data.5.2 QueryExample.txt

Expected Frequencies from DNAExample:
 - gaagc 2x (4 reported)
 - tcccc 0x (won't be printed to console)
 - tagat 5x (6 reported)
