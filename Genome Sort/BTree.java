import java.io.IOException;
import java.io.PrintWriter;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.text.DecimalFormat;

/**
 * Creates and manages a balanced-tree data structure
 *
 * @author Aaron Parish
 * @author Isaac Woodard
 *
 */
public class BTree {
	int t; // degree of BTree
	final int PAGE_SIZE = 4096; // this is the NODE_SIZE if the user specifies 0 as the degree
	final int NODE_SIZE;
	final int OPT_DEGREE = 45; //optimal degree for the tree
	long next; //the next empty position to create a node (also the end of the BTree file)
	RandomAccessFile rFile;
	BTreeNode root;
	Cache<BTreeNode> nodeCache;

	/**
 	* Constructor that creates a BTree with the given degree
 	*
 	* @param degree
 	* @param rFile
 	* @param cache
 	* @throws IOException
 	*/
	public BTree(int degree, RandomAccessFile rFile, int cacheSize) throws IOException {

    	nodeCache = new Cache<BTreeNode>(cacheSize);

    	if (degree == 0) {
        	NODE_SIZE = PAGE_SIZE;
        	this.t = 45; // the optimal degree
    	} else {
        	NODE_SIZE = findNodeSize(degree);
        	this.t = degree;
    	}
    	rFile.seek(0);
    	rFile.writeInt(t);
    	CreateBTree(rFile);
	}

	/**
 	* Constructor that creates a BTree from the given file
 	*
 	* @param rFile
 	* @param cache
 	* @param cacheSize
 	* @throws IOException
 	*/
	public BTree(RandomAccessFile rFile, int cacheSize) throws IOException {
			
		nodeCache = new Cache<BTreeNode>(cacheSize);

    	rFile.seek(0);
    	t = rFile.readInt();
    	if (t == 0) {
        	NODE_SIZE = PAGE_SIZE;
        	this.t = OPT_DEGREE; // the optimal degree
    	} else {
        	NODE_SIZE = findNodeSize(t);
    	}
    	this.rFile = rFile;
    	this.root = diskRead(rFile.readLong());
    	if (cacheSize != 0) { //add the root to the cache if it exists
    		nodeCache.addObject(root);
    	}
	}

	/**
 	* Helper constructor
 	*
 	* @param rFile
 	* @throws IOException
 	*/
	private void CreateBTree(RandomAccessFile rFile) throws IOException {
    	this.rFile = rFile;
    	next = 12;
    	BTreeNode x = new BTreeNode(); // Create empty root node x

    	x.pointer = allocateNode(); // Create pointer for node x
    	x.leaf = true;
    	x.numKeys = 0;
    	diskWrite(x); // Write changes to node x
    	root = x;
    	if(nodeCache.getSize() != 0) {
    		nodeCache.addObject(root);
    	}
	}

	/**
 	* Searches for a key in the tree and returns its frequency. Returns 0 if the
 	* key isn't found
	 * @throws IOException 
 	*/
	public long search(long key) throws IOException {
    	long frequency = BTreeSearch(root, key);

    	return frequency;
	}

	/**
 	* Helper method that searches for the specified key in the BTree. 
 	* Returns -1 if the key isn't found.
 	*
 	* @param current
 	* @param k
 	* @return
 	* @throws IOException
 	*/
	private long BTreeSearch(BTreeNode current, long k) throws IOException {

    	BTreeNode child = new BTreeNode();
    	int i = 0;

    	//search along the current node for the key
    	while (i < current.numKeys && k > current.treeArray[i].key) {
        	i++;
    	}
    	//check if we have found the key
    	if (i < current.numKeys && k == current.treeArray[i].key)
        	return (current.treeArray[i].frequency);
    	else if (current.leaf) { //return 0 frequency if we reach a leaf and haven't found the key
        	return 0;
    	} else { //if we aren't at a leaf, we read to a child and search it
        	child = diskRead(current.child[i]);
        	return BTreeSearch(child, k);
    	}
	}

	/**
 	* Helper method to insert into a node that isn't full
 	*
 	* @param x
 	* @param k
 	* @throws IOException
 	*/
	private void InsertNonFull(BTreeNode x, long k) throws IOException {
   	 
    	int i = x.numKeys;
    	if (i >= (2*t-1)) { //DEBUGGING
        	System.err.println("Error: tried to insert full node");
        	System.exit(0);
    	}
    	BTreeNode node = new BTreeNode();
   	 
    	//Check for duplicates (probably an inefficient implementation)
    	for (int j = 0; j < i; j++) {
        	if (k == x.treeArray[j].key) {
            	x.treeArray[j].incrementFreq();
            	diskWrite(x);
            	return;
        	}
    	}

    	if (x.leaf) {
        	while (i >= 1 && k < x.treeArray[i-1].key) {// from 0 to numKeys-1
            	x.treeArray[i].key = x.treeArray[i-1].key;
            	x.treeArray[i].frequency = x.treeArray[i-1].frequency;
            	i--;
        	}
        	x.treeArray[i].key = k;
        	x.treeArray[i].frequency = 1;
        	x.numKeys = (x.numKeys + 1);
        	diskWrite(x);
    	} else {
        	while (i >= 1 && k < x.treeArray[i-1].key) {// children from 0 to numKeys, keys from 0 to numKeys-1
            	i--;
        	}
    
        	node = diskRead(x.child[i]);
    
        	if (node.numKeys == ((2 * t) - 1)) {
            	SplitChild(x, i);
    
            	if (k > x.treeArray[i].key) {
                	i++;
            	}
        	}
        	InsertNonFull(diskRead(x.child[i]), k);
    	}
	}
    
	/**
 	* Inserts a key into the BTree. We first check if the root is
 	* full. If it is, we make a new root and make the old root a
 	* child of the new root. We then call InsertNonfull on the new
 	* root. If the original root is not full, we don't make a new
 	* root and instead call InsertNonFull on the original root.
 	* The method InsertNonFull should never be called on a root
 	* that is full.
 	*
 	* @param k
 	* @throws IOException
 	*/
	public void insert(long k) throws IOException {
    	BTreeNode currentRoot = root;
    	//Check if the root is full
    	if (currentRoot.numKeys == (2*t-1)){
        	//if root is full, we make a new node to be the root
        	BTreeNode newRoot = new BTreeNode();
        	newRoot.pointer = allocateNode();
        	newRoot.leaf = false;
        	newRoot.numKeys = 0;
        	newRoot.child[0] = currentRoot.pointer;
        	root = newRoot;
       	 
        	//since the old root is full, we need to split it
        	SplitChild(newRoot, 0);
       	 
        	//now we start our insertion at the new root
        	InsertNonFull(newRoot, k);
       	 
    	} else { //if root is not full, we start there
        	InsertNonFull(currentRoot, k);
    	}

	}

	/**
 	* Helper method that splits the full child of x at index
 	* i. Creates a new child of x at index i+1. The middle key
 	* of the child at i is moved to x, and the right half
 	* of the keys in that child are moved to the new child.
 	*
 	* @param parent
 	* @param c
 	* @throws IOException
 	*/
	private void SplitChild(BTreeNode parent, int c) throws IOException {
    	//Read in the full child and make the new child
    	BTreeNode oldChild = diskRead(parent.child[c]);
    	BTreeNode newChild = new BTreeNode();
    	newChild.pointer = allocateNode();
    	newChild.leaf = oldChild.leaf;
    	newChild.numKeys = t - 1;
   	 
    	//move second half of keys in old child to new child
    	for (int i = 0; i < t-1; i++) {
        	newChild.treeArray[i].key = oldChild.treeArray[i+t].key;
        	newChild.treeArray[i].frequency = oldChild.treeArray[i+t].frequency;
        	//remove values moved to newChild from oldChild
        	oldChild.treeArray[i+t].key = 0;
        	oldChild.treeArray[i+t].frequency = 0;
    	}
    	oldChild.numKeys = t-1;
   	 
    	//if old child is not a leaf, move second half of children from old child to new child
    	if (!oldChild.leaf) {
        	for (int i = 0; i < t; i++) {
            	newChild.child[i] = oldChild.child[i+t];
            	//Remove moved children from the old child
            	oldChild.child[i+t] = 0;
        	}
    	}
   	 
    	//make new child a child of parent at the index following the old child
    	for (int i = parent.numKeys; i > c; i--) { //we assume the parent is not full
        	parent.child[i+1] = parent.child[i];
    	}
    	parent.child[c+1] = newChild.pointer;
   	 
    	//make room for the middle key of old child and insert it into parent
    	for (int i = parent.numKeys-1; i > c-1; i--) { //again, we assume the parent is not full
        	parent.treeArray[i+1].key = parent.treeArray[i].key;
        	parent.treeArray[i+1].frequency = parent.treeArray[i].frequency;
    	}
    	parent.treeArray[c].key = oldChild.treeArray[t-1].key;
    	parent.treeArray[c].frequency = oldChild.treeArray[t-1].frequency;
    	oldChild.treeArray[t-1].key = 0; //remove value from old child
    	oldChild.treeArray[t-1].frequency = 0;
    	parent.numKeys = parent.numKeys + 1;
   	 
    	//write changes to disk
    	diskWrite(oldChild);
    	diskWrite(newChild);
    	diskWrite(parent);
   	 
	}

	/**
 	* Allocates space on file for a node.
 	* Returns the position to start a new
 	* node.
 	*
 	* @return next
 	*/
	private long allocateNode() {
    	long temp = next;
    	next += NODE_SIZE;
    	return temp;
	}

	/**
 	* Writes a BTreeNode to the disk. In particular, we write parent pointers,
 	* child pointers, keys and frequencies. Pointers are written as positions in
 	* the file.
 	*
 	* @param x
 	* @throws IOException
 	*/
	private void diskWrite(BTreeNode node) throws IOException {
   	 
    	rFile.seek(node.pointer); // Seek to position of node
       	rFile.write(toBytes(node)); // Passes node into toBytes()
   	 
    	// Node structure
    	// | numKeys | leaf | parent | child[i] | key[i] | frequency[i] |
	}
    
	/**
 	* Method for reading BTree nodes into a byte array
 	* @param node
 	* @return buffer.array()
 	*/
	public byte[] toBytes(BTreeNode node) {
   	 
    	int i;
    	byte leaf;
    	ByteBuffer buffer = ByteBuffer.allocate(NODE_SIZE);
    	buffer.putInt(node.numKeys); // Insert numKeys
   	 
    	if (node.leaf)
        	leaf = 1;
    	else
        	leaf = 0;
   	 
    	buffer.put(leaf);
   	 
    	buffer.putLong(node.parent);

    	for(i = 0; i < (2 * t); i++) {
        	buffer.putLong(node.child[i]); // Insert child pointer
    	}

    	for(i = 0; i < ((2 * t) - 1); i++){
        	buffer.putLong(node.treeArray[i].key); // Insert key
       	 
        	buffer.putLong(node.treeArray[i].frequency); // Insert frequency
    	}

    	return buffer.array();
	}

	/**
 	* Reads a BTreeNode from the disk. We pass an array of the longs to a BTreeNode
 	* in its constructor and it knows how to get its instance variables from the
 	* data. After creating the node, we store it into our node array or cache.
 	*
 	* @param i
 	* @return
 	* @throws IOException
 	*/
	private BTreeNode diskRead(long position) throws IOException {
		BTreeNode readNode;
		if (nodeCache.getSize() != 0) { //we first search the cache if it exists
			readNode = nodeCache.search(position);
			if (readNode != null) { //if we found the node, we return
				nodeCache.addObject(readNode);
				return readNode;
			}
		}
    	rFile.seek(position); // Seek to position of node on disk
       	byte[] byteArray = new byte[NODE_SIZE]; // Allocate byteArray to store data
       	rFile.read(byteArray); // Read entire node to byteArray
       	readNode = new BTreeNode(byteArray, position);
       	
       	//add node to cache if we had to read it and if it exists
       	if (nodeCache.getSize() != 0) {
       		nodeCache.addObject(readNode);
       	}

       	return readNode; // Create new node with byteArray data
	}

	/**
 	* Helper method for constructor that returns the size of a node on disk based
 	* on the degree of the bTree
 	*/
	private int findNodeSize(int degree) {
    	int size; // DISK_SIZE = Numkeys + isLeaf + Parent pointer + Child Pointers + Keys + Frequencies
                	// DISK_SIZE = x + (2t)y + (2t-1)z + (2t-1)z
                	// DISK_SIZE = 13 + 16t + 16t - 8 + 16t - 8
                	// DISK_SIZE = 48t - 3
    	size = 48 * degree - 3;
    	return size;
	}
    
	/**
 	* Prints an inOrder traversal of the BTree to the file with the given name
 	* @param name
 	* @param length
 	* @throws IOException
 	*/
	public void inOrder(String name, int length) throws IOException {
    	PrintWriter writer = new PrintWriter("dump", "UTF-8");    
   	 
    	traverse(writer, root, length);
           	 
    	writer.close();
	}
    
	/**
 	* Recursive in order traversal of the bTree
 	* @param writer
 	* @throws IOException
 	*/
	private void traverse(PrintWriter writer, BTreeNode node, int length) throws IOException {
    	if (node.leaf) {
        	for (int i=0; i < node.numKeys; i++) { //step through the keys of the leaf node
            	if (node.treeArray[i].key != 0) {
                	writer.println(node.treeArray[i].frequency + ", " + convertToString(node.treeArray[i].key, length)); //write to file <frequency>, <DNA string>
            	}
        	}
    	} else {
        	for (int i=0; i <= node.numKeys; i++) { //step through the children of the node
            	if (node.child[i] != 0) {
                	traverse(writer, diskRead(node.child[i]), length); //traverse one of the children
            	}
        	}
    	}
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
   	 
    	//check if we need to add leading "a"'s
    	while (sequence.length() < length) {
        	sequence.insert(0, 'a');
    	}
    	return sequence.toString();
	}
	
	/**
	 * This method must be called before closing the program.
	 * @throws IOException 
	 */
	public void writeRoot() throws IOException {
		rFile.seek(4);
		rFile.writeLong(root.pointer);
	}

	/**
 	* Class for the node of a balanced tree. Each node has a reference to its
 	* parent and children if any. It also holds a TreeObject for each key value and
 	* its frequency.
 	*
 	* @author Aaron Parish
 	* @author Isaac Woodard
 	*
 	*/
	public class BTreeNode {
   	 
    	TreeObject[] treeArray;
    	long[] child = new long[2 * t];
    	boolean leaf;
    	int numKeys;
    	long parent;
    	long pointer;

    	/**
     	* Constructor for the BTreeNode if the node isn't read from disk
     	*/
    	private BTreeNode() {
        	pointer = -1;
        	int i;
        	treeArray = new TreeObject[((2 * t) - 1)];
       	 
        	for(i = 0; i < ((2 * t) - 1); i++) { // Initialize all defaults to 0
            	treeArray[i] = new TreeObject(0, 0);
        	}
       	 
        	for(i = 0; i < (2 * t); i++) {
            	child[i] = 0;
        	}
       	 
        	numKeys = 0;
        	parent = 0;
        	leaf = true;
    	}

    	/**
     	* Constructor for the BTreeNode if the node is read from disk
     	*/
    	private BTreeNode(byte[] byteArray, long position) {
        	pointer = position;
       	 
        	treeArray = new TreeObject[((2 * t) - 1)];    
        	int i = 0;
        	ByteBuffer buffer = ByteBuffer.wrap(byteArray); // Using ByteBuffer on data from disk
        	this.numKeys = buffer.getInt(); // Read int - numKeys
        	this.leaf = 1 == buffer.get(); // Read boolean - leaf
        	this.parent = buffer.getLong(); // Read long - parent
       	 
            	for(i = 0; i < (2 * t); i++){   		
                	this.child[i] = buffer.getLong();
            	}
           	 
        	for(i = 0; i < ((2 * t) - 1); i++) { // Initialize all defaults to 0
           	treeArray[i] = new TreeObject(buffer.getLong(), buffer.getLong()); //Read key and frequency into Tree Object
        	}

        	// Node structure
        	// | numKeys | leaf | parent | child[i] | key[i] | frequency[i] |
    	}
    	
    	/**
    	 * Getter for the node's position
    	 * 
    	 * @return pointer
    	 */
    	public long getSelf() {
    		return pointer;
    	}
	}
	
	/**
	 * Creates a cache storage object. The cache stores frequently searched for objects. When an object is searched for
	 * it is moved to the top of the cache. When an object is searched for in the cache and not found, it is added to
	 * the cache and moved to the top. The cache will discard objects at the bottom (least used) rather than go over 
	 * capacity when adding new objects.
	 * 
	 * @author Isaac Woodard
	 * @param <T>
	 *
	 */
	public class Cache<T> {

		private int References;
		private int Hits;
		final int SIZE;
		private IUDoubleLinkedList<Object> list;
		
		/**
		 * Constructor for the Cache. Creates a cache to hold
		 * generic objects.
		 * 
		 * @param size
		 */
		public Cache(int size) {
			SIZE = size;
			References = 0;
			Hits = 0;
			list = new IUDoubleLinkedList<Object>();
		}
		
		/**
		 * Searches for and returns the specified object.
		 * Adds the object to the list at the top of the
		 * cache if it isn't found.
		 * 
		 * @param object
		 * @return
		 */
		public Object getObject(Object obj) {
			Object retval;
			References++;
			retval = null;
			
			if (list.contains(obj)) {
				Hits++;
				retval = list.remove(obj);
			}
			if (list.size() == SIZE) {
				list.removeLast();
			}
			list.addToFront(obj);
		
			return retval;
		}
		
		/**
		 * Searches the cache for a BTreeNode with the given 
		 * position. Returns null if the node isn't found.
		 * 
		 * @param position
		 * @param node
		 */
		public BTreeNode search(long position) {
			BTreeNode node;
			
			//search through the double linked list
			int i = 0;
			while (i < list.size()) {
				node = (BTreeNode) list.get(i);
				if (node.pointer == position) {
					return node;
				}
				i++;
			}
			node = null;
			
			return node;
		}
		
		/**
		 * Adds an object to the cache at the top.
		 * Removes the last object if the cache is
		 * full. Removes the object then adds it
		 * if the object is already in the list.
		 * 
		 * @param object
		 */
		public void addObject(Object obj) {
			if (list.contains(obj)) {
				list.remove(obj);
			}
			if (list.size() == SIZE) {
				list.removeLast();
			}
			list.addToFront(obj);
		}
		
		/**
		 * Removes an object from the cache and
		 * returns it.
		 * 
		 * @param object
		 */
		public Object removeObject(Object obj) {
			Object retval = list.remove(obj);
			return retval;
		}
		
		/**
		 * Removes all objects from the cache's memory.
		 */
		public void clearCache() {
			list = new IUDoubleLinkedList<Object>();
		}

		/**
		 * Gets the cache size.
		 * 
		 * @return cache size
		 */
		public int getSize() {
			return SIZE;
		}
		
		/**
		 * Gets the number of references.
		 * 
		 * @return References
		 */
		public int getReferences() {
			return References;
		}
		
		/**
		 * Gets the number of hits.
		 * 
		 * @return Hits
		 */
		public int getHits() {
			return Hits;
		}
		
		/**
		 * Gives the ratio of References to
		 * Hits.
		 * 
		 * @return hit ratio
		 */
		public double getHitRatio() {
			return (double)Hits/References;
		}
		/**
		 * Reports the cache size, number of references,
		 * number of hits, and hit ratio.
		 */
		public String toString() {
			StringBuilder str = new StringBuilder();
			DecimalFormat f = new DecimalFormat("#.00");	
			
			str.append("Size: " + SIZE + "\n");
			str.append("References: " + References + "\n");
			str.append("Hits: " + Hits + "\n");
			str.append("Hit Ratio: " + f.format(this.getHitRatio()) + "\n");
			
			return str.toString();
		}
	}

}
