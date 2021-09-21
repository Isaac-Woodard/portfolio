/**
 * This class creates a node to hold an element for use in a double linked list. Each
 * node has a reference to the node after it and the node before it.
 * 
 * @author Isaac Woodard & CS-221
 *
 * @param <T> type to store
 */
public class LinearNode2<T> {
	private T element;
	private LinearNode2<T> next;
	private LinearNode2<T> prev;
	
	/**
	 * This constructor creates an empty node with no next or prev reference.
	 */
	public LinearNode2() {
		element = null;
		next = prev = null;
	}
	
	/**
	 * This constructor creates a node with a given element and no next or prev reference.
	 * 
	 * @param element
	 */
	public LinearNode2(T element) {
		this.element = element;
		next = prev = null;
	}
	
	/**
	 * This constructor creates a node with a given element and a next and prev reference to
	 * the given nodes.
	 * 
	 * @param element
	 * @param nextNode
	 */
	public LinearNode2(T element, LinearNode2<T> nextNode, LinearNode2<T> prevNode) {
		this.element = element;
		next = nextNode;
		prev = prevNode;
	}
	
	/**
	 * Returns the element stored in the node.
	 * 
	 * @return the element
	 */
	public T getElement() {
		return element;
	}
	
	/**
	 * Sets the element stored in the node.
	 * 
	 * @param element
	 */
	public void setElement(T element) {
		this.element = element;
	}
	
	/**
	 * Returns the reference to the next node.
	 * 
	 * @return the next node reference
	 */
	public LinearNode2<T> getNext() {
		return next;
	}
	
	/**
	 * Sets the reference to the next node.
	 * 
	 * @param next
	 */
	public void setNext(LinearNode2<T> next) {
		this.next = next;
	}
	
	/**
	 * Returns the reference to the previous node.
	 * 
	 * @return the previous node reference
	 */
	public LinearNode2<T> getPrev() {
		return prev;
	}
	
	/**
	 * Sets the reference to the previous node.
	 * 
	 * @param prev
	 */
	public void setPrev(LinearNode2<T> prev) {
		this.prev = prev;
	}
	
	@Override
	public String toString() {
		return element.toString();
	}
}
