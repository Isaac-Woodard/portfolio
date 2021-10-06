import java.util.ConcurrentModificationException;
import java.util.Iterator;
import java.util.ListIterator;
import java.util.NoSuchElementException;

/**
 * This is a double linked list implementation of the IndexedUnsortedList
 * interface. It creates and manages a list of elements of type T.
 * 
 * @author Isaac Woodard & CS 221-2
 *
 * @param <T> type to store
 */
public class IUDoubleLinkedList<T> implements IndexedUnsortedList<T> {

	private LinearNode2<T> head;
	private LinearNode2<T> tail;
	private int size;
	private int modCount;

	/**
	 * The constructor for the class. Creates an empty list.
	 */
	public IUDoubleLinkedList() {
		head = tail = null;
		size = 0;
		modCount = 0;
	}

	@Override
	public void addToFront(T element) {
		LinearNode2<T> newNode = new LinearNode2<T>(element);
		if (isEmpty()) {
			head = tail = newNode;
		} else {
			newNode.setNext(head);
			head.setPrev(newNode);
			head = newNode;
		}
		size++;
		modCount++;
	}

	@Override
	public void addToRear(T element) {
		LinearNode2<T> newNode = new LinearNode2<T>(element);
		if (isEmpty()) {
			head = tail = newNode;
		} else {
			newNode.setPrev(tail);
			tail.setNext(newNode);
			tail = newNode;
		}
		size++;
		modCount++;
	}

	@Override
	public void add(T element) {
		addToRear(element);
	}

	@Override
	public void addAfter(T element, T target) {
		//consider writing this to track the current cursor value rather than foundTarget
		ListIterator<T> lit = listIterator();
		boolean foundTarget = false;
		while (lit.hasNext() && !foundTarget) {
				if (lit.next().equals(target)) {
					foundTarget = true;
				}
		}
		if (!foundTarget) {
			throw new NoSuchElementException();
		}
		lit.add(element);
	}

	@Override
	public void add(int index, T element) {
		if (index < 0 || index > size) {
			throw new IndexOutOfBoundsException();
		}

		if (index == 0) {
			addToFront(element);
		} else if (index == size) {
			addToRear(element);
		} else {
			LinearNode2<T> newNode = new LinearNode2<T>(element);
			LinearNode2<T> current = head;
			for (int i = 0; i < index; i++) {
				current = current.getNext(); // inserted thing will go in front
												// of current
			}
			newNode.setNext(current);
			newNode.setPrev(current.getPrev());
			current.setPrev(newNode);
			newNode.getPrev().setNext(newNode);

			size++;
			modCount++;
		}
	}

	@Override
	public T removeFirst() {
		if (size == 0) {
			throw new IllegalStateException();
		}
		
		T retVal = head.getElement();
		head = head.getNext();
		if (head == null) {
			tail = null;
		} else {
			head.setPrev(null);
		}
		
		size--;
		modCount++;
		return retVal;
	}

	@Override
	public T removeLast() {
		if (size == 0) {
			throw new IllegalStateException();
		}
		
		T retVal = tail.getElement();
		tail = tail.getPrev();
		if (tail == null) {
			head = null;
		} else {
			tail.setNext(null);
		}
		
		size--;
		modCount++;
		return retVal;
	}

	@Override
	public T remove(T element) {
		LinearNode2<T> current = head;
		while (current != null && !current.getElement().equals(element)) {
			current = current.getNext();
		}
		if (current == null) {
			throw new NoSuchElementException();
		}

		if (current == head) {
			head = head.getNext();
			if (head != null) {
				head.setPrev(null);
			}
		} else {
			current.getPrev().setNext(current.getNext());
		}
		if (current == tail) {
			tail = tail.getPrev();
			if (tail != null) {
				tail.setNext(null);
			}
		} else {
			current.getNext().setPrev(current.getPrev());
		}

		size--;
		modCount++;

		return current.getElement();
	}

	@Override
	public T remove(int index) {
		if (index < 0 || index >= size) {
			throw new IndexOutOfBoundsException();
		}

		ListIterator<T> lit = listIterator(index);
		T retVal = lit.next();
		lit.remove();

		return retVal;
	}

	@Override
	public void set(int index, T element) {
		if (index < 0 || index >= size) {
			throw new IndexOutOfBoundsException();
		}
		ListIterator<T> lit = listIterator(index);
		lit.next();
		lit.set(element);
	}

	@Override
	public T get(int index) {
		if (index < 0 || index >= size) {
			throw new IndexOutOfBoundsException();
		}
		ListIterator<T> lit = listIterator(index);
		return lit.next();
	}

	@Override
	public int indexOf(T element) {
		int index = -1;
		int i = 0;
		ListIterator<T> lit = listIterator();
		while (index < 0 && lit.hasNext()) {
			if (lit.next().equals(element)) {
				index = i;
			}
			i++;
		}
		return index;
	}

	@Override
	public T first() {
		if (isEmpty()) {
			throw new IllegalStateException();
		}
		return head.getElement();
	}

	@Override
	public T last() {
		if (isEmpty()) {
			throw new IllegalStateException();
		}
		return tail.getElement();
	}

	@Override
	public boolean contains(T target) {
		LinearNode2<T> current = head;
		while (current != null) {
			if (current.getElement().equals(target)) {
				return true;
			}
			current = current.getNext();
		}
		return false;
	}

	@Override
	public boolean isEmpty() {
		return (head == null);
	}

	@Override
	public int size() {
		return size;
	}

	@Override
	public Iterator<T> iterator() {
		return listIterator();
	}

	@Override
	public ListIterator<T> listIterator() {
		return new DLLIterator();
	}

	@Override
	public ListIterator<T> listIterator(int startingIndex) {
		return new DLLIterator(startingIndex);
	}
	
	@Override
	public String toString() {
		StringBuilder str = new StringBuilder();
		
		str.append("[");
		
		for (T element : this) {
			str.append(element);
			str.append(", ");
		}
		if (size > 0) {
			str.delete(str.length()-2, str.length());
		}
		
		str.append("]");
		
		return str.toString();
	}

	/**
	 * ListIterator for IUDoubleLinkedList. It can step forward or backward in
	 * the list as well as add and remove elements from the list.
	 * 
	 * @author Isaac Woodard & CS 221-2
	 *
	 */
	private class DLLIterator implements ListIterator<T> {
		private LinearNode2<T> nextNode;
		private int nextIndex;
		private LinearNode2<T> lastReturned;
		private int iterModCount;

		/**
		 * The basic constructor. Creates a list iterator starting at the
		 * beginning of the list.
		 */
		public DLLIterator() {
			this(0);
		}

		/**
		 * This constructor creates a list iterator starting before the
		 * specified index. If the startingIndex is equal to the list's
		 * size, the list iterator starts at the end of the list.
		 * 
		 * @param startingIndex
		 */
		public DLLIterator(int startingIndex) {
			if (startingIndex < 0 || startingIndex > size) {
				throw new IndexOutOfBoundsException();
			}
			nextNode = head;
			lastReturned = null;
			nextIndex = 0;
			iterModCount = modCount;
			for (int i = 0; i < startingIndex; i++) {
				nextNode = nextNode.getNext();
				nextIndex++;
			}
		}

		@Override
		public void add(T element) { // [previous] [new element] ^ [next]
			if (iterModCount != modCount) {
				throw new ConcurrentModificationException();
			}

			LinearNode2<T> newNode = new LinearNode2<T>(element);

			if (size == 0) {
				head = tail = newNode;
			} else if (nextNode == head) { // new head
				newNode.setNext(head);
				head.setPrev(newNode);
				head = newNode;
			} else if (nextNode == null) { // new tail
				newNode.setPrev(tail);
				tail.setNext(newNode);
				tail = newNode;
			} else {
				newNode.setNext(nextNode);
				newNode.setPrev(nextNode.getPrev());
				nextNode.setPrev(newNode);
				newNode.getPrev().setNext(newNode);
			}

			lastReturned = null;
			size++;
			modCount++;
			iterModCount++;
			nextIndex++;
		}

		@Override
		public boolean hasNext() {
			if (iterModCount != modCount) {
				throw new ConcurrentModificationException();
			}
			return (nextNode != null);
		}

		@Override
		public boolean hasPrevious() {
			if (iterModCount != modCount) {
				throw new ConcurrentModificationException();
			}
			return (nextNode != head);
		}

		@Override
		public T next() {
			if (!hasNext()) {
				throw new NoSuchElementException();
			}
			T retVal = nextNode.getElement();
			lastReturned = nextNode;
			nextNode = nextNode.getNext();

			nextIndex++;
			return retVal;
		}

		@Override
		public int nextIndex() {
			if (iterModCount != modCount) {
				throw new ConcurrentModificationException();
			}
			return nextIndex;
		}

		@Override
		public T previous() {
			if (!hasPrevious()) {
				throw new NoSuchElementException();
			}
			if (nextNode == null) {
				nextNode = tail;
			} else {
				nextNode = nextNode.getPrev();
			}

			nextIndex--;
			lastReturned = nextNode;
			return nextNode.getElement();
		}

		@Override
		public int previousIndex() {
			if (iterModCount != modCount) {
				throw new ConcurrentModificationException();
			}
			return nextIndex - 1;
		}

		@Override
		public void remove() {
			if (iterModCount != modCount) {
				throw new ConcurrentModificationException();
			}
			if (lastReturned == null) {
				throw new IllegalStateException();
			}

			if (lastReturned == head) {
				head = head.getNext();
				if (head != null) {
					head.setPrev(null);
				}
			} else {
				lastReturned.getPrev().setNext(lastReturned.getNext());
			}
			if (lastReturned == tail) {
				tail = tail.getPrev();
				if (tail != null) {
					tail.setNext(null);
				}
			} else {
				lastReturned.getNext().setPrev(lastReturned.getPrev());
			}

			if (nextNode != lastReturned) { // last move was next
				nextIndex--;
			} else { // last move was prev
				nextNode = lastReturned.getNext();
			}

			lastReturned = null;
			size--;
			modCount++;
			iterModCount++;
		}

		@Override
		public void set(T element) {
			if (iterModCount != modCount) {
				throw new ConcurrentModificationException();
			}
			if (lastReturned == null) {
				throw new IllegalStateException();
			}

			lastReturned.setElement(element);

			modCount++;
			iterModCount++;
		}

	}
}
