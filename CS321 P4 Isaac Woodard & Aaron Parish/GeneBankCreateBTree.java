import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.RandomAccessFile;
import java.util.Scanner;

/**
 * Converts a genebank file (.gbk) into a BTree. Sequences of DNA
 * are represented as long values.
 * 
 * @author Isaac Woodard
 *
 */
public class GeneBankCreateBTree {
	public static void main(String[] args) {//<0/1(no/with Cache)> <degree> <gbk file> <sequence length> [<cache size>] [<debug level>]
		try {
			//Parse and process command line arguments	//TODO: put command line parsing into a method
			int cacheQ = Integer.parseInt(args[0]);
			boolean cache = false;
			int cacheSize = 0;
			int degree = Integer.parseInt(args[1]);
			File geneBank = new File(args[2]);
			int subLength = Integer.parseInt(args[3]);
			int debugLevel = 0;
			if (subLength < 1 || subLength > 31) { //checks if the subsequence length is valid
				System.err.println("Error: invalid subsequence length");
				System.exit(0);
			}
			if (cacheQ == 1) { //checks for if the user wants a cache
				cacheSize = Integer.parseInt(args[4]);
				cache = true;
			}
			if (cacheQ == 0 && args.length > 4) { //other check for if the user wants to debug
				debugLevel = Integer.parseInt(args[4]);
			}
			if (args.length > 5) { //checks for if the user wants to debug
				debugLevel = Integer.parseInt(args[5]);
			}
		
			//Make the btree file on disk
			String fileName = args[2] + ".btree.data." + subLength + "." + degree; //format: xyz.gbk.btree.data.k.t.
			File file = new File(fileName);
			if (file.exists()){
			    file.delete();
			}
			RandomAccessFile rFile = new RandomAccessFile(fileName, "rw"); //TODO: does this actually make a file?
			
			//Parse genebank file and store sub sequences in BTree
			BTree bTree = new BTree(degree, rFile, cacheSize);
			parseGBK(geneBank, subLength, bTree);
			bTree.writeRoot();
			
			if (debugLevel == 1) {
				try{
					bTree.inOrder("dump", subLength);				
				} catch (IOException e) {
					e.printStackTrace();
					System.exit(0);
				}
			}
			rFile.close();
		} catch (NumberFormatException e) { //catches command line parsing problems
			printUsage();
		} catch (FileNotFoundException e1) {
			e1.printStackTrace();
		} catch (IOException e1) {
			e1.printStackTrace();
		}
	}
	
	/**
	 * Parses the genebank file and inserts subsequences
	 * into the BTree
	 * @param gbk
	 * @throws IOException 
	 */
	private static void parseGBK(File gbk, int subLength, BTree bTree) throws IOException {
		Scanner fileScanner = new Scanner(gbk); //TODO: the class BufferReader would be more efficient than a scanner
		StringBuilder mainSeq = new StringBuilder();
		StringBuilder subSeq = new StringBuilder();
		Boolean append = false;
		final int LINES = 20; //Number of lines to read from file at a time
		final int CARRY = 100; //Number of characters to carry over between buffer refills
		//int counter = 1; //DEBUG: counts number of sequences
		//int ncounter = 0; //DEBUG: counts number of n and N characters
		
		while (fileScanner.hasNextLine()) {
			//skip to next "ORIGIN"
			if (append == false) {
				while (fileScanner.hasNextLine() && !fileScanner.nextLine().contains("ORIGIN")); //go to start of dna sequence
				append = true;
			}
			//read up to LINES number of lines
			int i = 0;
			String tempS = new String();
			while (i < LINES && fileScanner.hasNextLine() && append == true) {
				tempS = fileScanner.nextLine();
				mainSeq.append(tempS);
				if (tempS.contains("//")) { //end of dna sequence
					append = false;
				}
				i++;
			}
			//get subsequences
			char tempC;
			for (i = 0; i < (append ? mainSeq.length()-CARRY : mainSeq.length()); i++) {
				int j = 0; //counter for the number of characters found for the subsequence
				int k = 0; //counter for the increment from position i for the inner loop
				tempC = mainSeq.charAt(i);
				if(!(isDNA(tempC))) { //if the current char isn't DNA, move on without going into the inner loop
					j = subLength;
				}
				while (j < subLength && i+k < mainSeq.length()) {
					tempC = mainSeq.charAt(i+k);
					if (isBreaker(tempC)) { //break at '/' or 'N'
						j = subLength;
						//ncounter++; //DEBUGGING
					}
					if(isDNA(tempC)) {//ignore numbers, spaces, newlines
						subSeq.append(tempC);
						j++;
					}
					k++;
				}
				//if subsequence is complete, insert to BTree
				if (subSeq.length() == subLength) {
					long binary = convertToBinary(subSeq.toString());
					bTree.insert(binary);
					//System.out.println("(" + counter + ") " + binary + ", " + subSeq.toString()); //DEBUGGING
					//counter++; //DEBUGGING
				}
				subSeq.delete(0, subSeq.length()); //the endpoint is exclusive
			}
			//refresh buffer
			if (mainSeq.length() > CARRY) //to avoid deleting a negative range (and throwing an exception)
				mainSeq.delete(0, mainSeq.length()-CARRY);
			if (append == false) //only delete everything if we reached the end of the sequence
				mainSeq.delete(0, mainSeq.length());
		}
		fileScanner.close();
		//System.out.println(ncounter + ", ");//DEBUGGING
	}
	
	/**
	 * Converts the given string of DNA to binary in
	 * the form of a long.
	 * 
	 * @param subsequence
	 * @return long
	 */
	private static long convertToBinary(String subsequence) {
		long key = 0;
		char chemical;
		StringBuilder binary = new StringBuilder();
		
		//Goes through each character in the subsequence
		for (int i = 0; i < subsequence.length(); i++) {
			chemical = subsequence.charAt(i);
			switch(chemical) {
				case 'A': //00
				case 'a':
					binary.append("00");
					break;
				case 'T': //11
				case 't':
					binary.append("11");
					break;
				case 'C': //01
				case 'c':
					binary.append("01");
					break;
				case 'G': //10
				case 'g':
					binary.append("10");
					break;
				default:
					System.err.print("ERROR: Invalid character in subsequence");
					break;	
			}
		}
		key = Long.parseLong(binary.toString(), 2);
		return key;
	}
	
	/**
	 * Returns true if the character is '/' or 'n' or 'N'
	 * 
	 * @param c
	 * @return boolean
	 */
	private static boolean isBreaker(char c) {
		boolean check = false;
		if (c == '/' || c == 'N' || c == 'n')
			check = true;
		return check;
	}
	
	/**
	 * 
	 * Returns true if the character is a lower or upper case
	 * letter of DNA
	 * 
	 * @param c
	 * @return boolean
	 */
	private static boolean isDNA(char c) {
		boolean check = false;
		if(c == 'A' || c == 'a' || c == 'T' || c == 't' || c == 'C' || c == 'c' || c == 'G' || c == 'g')
			check = true;
		return check;
	}
	
	/**
	 *  Prints the usage statement for the program
	 */
	private static void printUsage() {
		System.err.println("Usage: <0/1(no/with Cache)> <degree> <gbk file> <sequence length> [<cache size>] [<debug level>]");
	}
}
