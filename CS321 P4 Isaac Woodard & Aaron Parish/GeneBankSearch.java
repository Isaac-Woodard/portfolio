import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.Scanner;

/**
 * Searches a BTree file for the frequency of a certain subsequence
 * of DNA.
 * 
 * @author Isaac Woodard
 *
 */
public class GeneBankSearch {
	public static void main(String[] args) {// <0/1(no/with Cache)> <btree file> <query file> [<cache size>] [<debug level>]
		try {
			//Parse and process command line arguments	
			int cacheQ = Integer.parseInt(args[0]);
			boolean cache = false;
			int cacheSize = 0;
			//File bFile = new File(args[1]); //TODO: is this needed for RandomAccessFile?
			File qFile = new File(args[2]);
			int debugLevel = 0;
			if (cacheQ == 1) { //checks for if the user wants a cache
				cacheSize = Integer.parseInt(args[3]);
				cache = true;
			}
			if (cacheQ == 0 && args.length > 3) { //other check for if the user wants to debug
				debugLevel = Integer.parseInt(args[4]);
			}
			if (args.length > 4) { //checks for if the user wants to debug
				debugLevel = Integer.parseInt(args[5]);
			}
		
			//Search for each query
			RandomAccessFile rFile = new RandomAccessFile(args[1], "rw");
			if (rFile.length() == 0) { //check that the btree file isn't empty
				System.exit(0);
			}
			BTree bTree = new BTree(rFile, cacheSize);
			Scanner qScanner = new Scanner(qFile);
			while (qScanner.hasNextLine()) {
				String dna = qScanner.nextLine();
				long key = convertToBinary(dna);
				long frequency = bTree.search(key);
				if (frequency != 0) { //only print non-zero freqencies
					System.out.println(dna + ": " + frequency);
				}
			}
		} catch (NumberFormatException e) { //catches command line formatting problems
			printUsage();
		} catch (FileNotFoundException e1) {
			e1.printStackTrace();
		} catch (IOException e2) {
			e2.printStackTrace();
		}
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
	 * Converts the binary of a long into a sequence of dna.
	 * 
	 * @param key
	 * @param length of desired sequence
	 * @return sequence
	 */
	private static String convertToString(long key, int length) {
		StringBuilder sequence = new StringBuilder();
		String binary = Long.toBinaryString(key);
		String twoBits;
		String letter;
		
		//add a leading 0 if binary string length is odd
		if (binary.length() % 2 != 0) {
			binary = 0 + binary;
		}
		
		//Goes through the subsequence two bits at a time
		for (int i = 0; i < binary.length(); i+=2) {
			twoBits = binary.substring(i, i+2);
			switch(twoBits) {
				case "00": //A or a
					letter = "a";
					break;
				case "11": //T or t
					letter = "t";
					break;
				case "01": //C or c
					letter = "c";
					break;
				case "10": //G or g
					letter = "g";
					break;
				default:
					letter = "";
					System.err.println("ERROR: Invalid character in subsequence");
					//System.out.println(twoBits); //DEBUG
					break;	
			}
			sequence.append(letter);
		}
		
		//check if we need to add a leading "a"
		if (sequence.length() < length) {
			sequence.insert(0, 'a');
		}
		return sequence.toString();
	}
	
	/**
	 * Prints the usage statement for the program
	 */
	private static void printUsage() {
		System.err.println("Usage: <0/1(no/with Cache)> <btree file> <query file> [<cache size>] [<debug level>]");
	}
}
