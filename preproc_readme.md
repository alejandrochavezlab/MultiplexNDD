How to get count data from a barcode sequencing experiment. 
Assumptions: Unix (mac or linux) based machine. RStudio. Usegalaxy.org. 

1. Download the FastQ files according to the indexes you used into a folder. Unzip.
2. Concatenate all the different files into 1 master file. (Unix)

	cat folderpath/*.fastq > folderpath/masterfile.fastq

3. To keep everything simple, we want to make sure that we are only using reads that are the
same length. If we ordered 75 bp reads we expect 75 bp. You may want to check that the reads
are truly 75 bp. 

	awk 'BEGIN {OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) == 75) {print header, seq, qheader, qseq}}' < masterfile.fastq > filteredmasterfile.fastq

4. Now you want to extract the useful information from the fastq file. There are usually
three useful pieces of information in a fastq file for this purpose. One is the illumina 
index (found in the first line of every 4 lines). Another is the internal index during
1st round amplification (found in the 2nd line of every 4 lines). The last is the barcode 
sequence (alse found in 2nd line). I do this in Unix.

I use this paste command which you see below. You want to get an idea of how the paste command makes your file look:

	head filteredmasterfile.fastq | paste - - - -
	
This command splits every 4 lines into a different column. the command awk which we use below sees "fields"
"fields" are separated by spaces or tabs or blank space. So the 4 lines turn into different fields
"Fields" can be indexed by a dollarsign and a number as you see below
Also below this you will see a substr command, which takes only a piece of a string of letters.
This substring command takes 3 arguments (the field, the place you want to start, and the number of places after where you start to extract)

	cat filteredmasterfile.fastq | paste - - - - | awk '{print substr($2, 7, 8), substr($3, 5,5), substr($3, 50,20)}' > filtered_col_extract.txt

cat /home/pkrupav/Horfeome/NGS_Analysis/03Sept2024_Analysis/RawOutputs_Zip/ORF17_SingleRound_Val_S34_L002_R1_001.fastq | paste - - - - | awk '{print substr($3, 5, 5), substr($3, 101, 4), substr($3, 41, 20)}' > ORF17_SingleRound_Extract.txt
	
So the above command prints 3 things for each read. One is the illumina index which is in the second field
and starts at position 7 and is 8 characters long [substr($2, 7, 8)]. Second it prints the internal
index which is in the 3rd field and starts at position 5 and is 5 characters long [substr($3, 5,5)].
Lastly it takes a barcode which in this case follows the same example.
Caution: if you samples are double indexed you will have to add in an extra step but i trust you can figure that out.

5. Great, now we have a file (filtered_col_extract.txt) that has 1 row for each read and 3 columns.
We want to be able to attribute each read to a specific well that it came from and we want to annotate that.
My wells are defined by some combination of the internal index and the illumina index. So I simply
combine the 2 letter combinations into 1.

awk '{print $1$2"\t"$3}' < /home/pkrupav/Horfeome/NGS_Analysis/03Sept2024_Analysis/AnalysisOutputs/ORF17_SingleRound_Extract.txt > /home/pkrupav/Horfeome/NGS_Analysis/03Sept2024_Analysis/AnalysisOutputs/ORF17_SingleRound_Extract_Merged.txt

6. Now we have a very simple two column file. And we want to replace the combination of the 
indexes with a name that is easier to read and carries information about the well. In the
case of a drug screen you might want to replace the index combo with a name like "Cyclohexamide"
that way you can bin all the reads attributed to the cyclohexamide condition. We turn to R for this
because it is much better at this sort of data merging. Open RStudio. Start a script.

seqs <- read.table("pathtofile/filtered_col_merged.txt", header = F, sep = "\t")
colnames(seqs) <- c("merged_index", "seq")
seqs <- data.frame(lapply(seqs, as.character), stringsAsFactors=FALSE)

The above commands reads in your text file, gives it column names, and then converts all the 
items in the data frame into characters (a data type we need them to be in).

7. Now you want to (in excel) define what combinations of indexes mean what well. like below
HSPH1_glu_1	CGTACTAGAATGC
HSPA13_glu_1	CGTACTAGGTTGC
HSPA9_glu_1	CGTACTAGGATCA
HSPA14_glu_1	CGTACTAGTGATC

you want to save this is a .csv file. (well_ids.csv)

8. You want to read this file into Rstudio.

well_id <- read.csv(file = "pathtofile/well_ids.csv")
colnames(well_id) <- c("ID", "merged_index")
well_id$merged_index <- as.character(well_id$merged_index)

9. Now you want to add this information to the seqs dataset. we use the merge command
reads <- merge(seqs, well_id)

10. This should have multiple columns but you really only care about two. the sequence and
the well ID. I will just select out these two things.
reads <- data.frame(reads$ID, reads$seq)

11. I want to write this out to a text file so I can convert into a fasta file which will
be used to align the sequenced DNA to the reference sequences (the list of barcodes and what they mean)

write.table(reads, file= "pathtofile/well_indexed_reads.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

12. I want to convert that file to a fasta file. Back to Unix.

awk '{ printf ">%s\n%s\n",$1,$2 }' well_indexed_reads.txt > well_indexed_reads.fasta

13. Now we have our fasta file. and we want to make our reference fasta file. This is more
straight forward and two steps because presumably there are not more than 1 million barcodes.
You want to make the list in excel. An example looks like this: 
426-GAL-mCHERRY_A20	TAGGTCATCCTCAGCGGTCT
426-GAL-mCHERRY_A297	ATGATGCGGTCACCTTCCGT
416-GAL-AMN1_A137	CAACTCCGAGGGATGTCAGT
416-GAL-AMN1_A41	CTGCAACTGATCGTCATTCT
416-GAL-AMN1_A42	GCAGCAGCGTGTACGATTCT

And this is a .csv file called reference.csv

14. Read this into Rstudio

ref <- read.csv(file = "pathtofile/reference.csv", header = T, stringsAsFactors = FALSE)
ref <- read.csv(file = "C:\Users\moldev\OneDrive\Desktop\Chavez_Lab\hOrfeome Screen\ReferenceDocs\7May2024_BarcodeReference.csv", header = T, stringsAsFactors = FALSE)

15. Then write it out as a fasta file.

install.packages("seqRFLP") 
library(seqRFLP)
reffasta <- dataframe2fas(ref, file = "pathtofile/reference.fasta")

16. Now we move on to galaxy which hosts a cloud version of an alignment software. Upload
the well_indexed_reads.fasta and the reference.fasta. Go to NGS alignment and bowtie2.
Run bowtie2 with defaults but building your own reference.

17. This should take a while to run. take a break and when it is done you will see that
you can download a .bam file. This is what you want. Download it to somewhere and open R.
In r you need to import some packages that allow you to interface with bams.