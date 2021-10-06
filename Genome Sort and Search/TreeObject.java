/**
 * Holds a key value and a frequency counter for
 * a BTreeNode
 * @author Isaac Woodard
 *
 */
public class TreeObject {
	long key;
	long frequency;
	
	/**
	 * Creates the object with a specified key
	 * and frequency
	 */
	public TreeObject(long key, long frequency) {
		this.key = key;
		this.frequency = frequency;
	}

	/**
	 * Getter for the key
	 * @return key
	 */
	public long getKey() {
		return key;
	}

	/**
	 * Setter for the key
	 * @param key
	 */
	public void setKey(long key) {
		this.key = key;
	}
	
	/**
	 * Increments the frequency by 1
	 */
	public void incrementFreq() {
		frequency++;
	}
}
