/*
 * Our min-heap implementation.
 *
 * Author (starter code): A. Tafliovich.
 */

#include "minheap.h"

#define ROOT_INDEX 1
#define NOTHING -1

/*************************************************************************
 ** Suggested helper functions -- to help designing your code
 *************************************************************************/
 
 /* Returns True if 'maybeIdx' is a valid index in minheap 'heap', and 'heap'
 * stores an element at that index. Returns False otherwise.
 */
bool isValidIndex(MinHeap* heap, int maybeIdx) {
	return (maybeIdx >= ROOT_INDEX && maybeIdx <= heap->size);
}

/* Returns node at index 'nodeIndex' in minheap 'heap'.
 * Precondition: 'nodeIndex' is a valid index in 'heap'
 *               'heap' is non-empty
 */
HeapNode nodeAt(MinHeap* heap, int nodeIndex) {
	return heap->arr[nodeIndex];
}

/* Returns priority of node at index 'nodeIndex' in minheap 'heap'.
 * Precondition: 'nodeIndex' is a valid index in 'heap'
 *               'heap' is non-empty
 */
int priorityAt(MinHeap* heap, int nodeIndex) {
	return heap->arr[nodeIndex].priority;
}

/* Returns ID of node at index 'nodeIndex' in minheap 'heap'.
 * Precondition: 'nodeIndex' is a valid index in 'heap'
 *               'heap' is non-empty
 */
int idAt(MinHeap* heap, int nodeIndex) {
	return heap->arr[nodeIndex].id;
}

/* Returns index of node with ID 'id' in minheap 'heap'.
 * Precondition: 'id' is a valid ID in 'heap'
 *               'heap' is non-empty
 */
int indexOf(MinHeap* heap, int id) {
	return heap->indexMap[id];
}

/* Returns the index of the left child of a node at index 'nodeIndex' in
 * minheap 'heap', if such exists.  Returns NOTHING if there is no such left
 * child.
 */
int leftIdx(MinHeap* heap, int nodeIndex) {
	if (!isValidIndex(heap, nodeIndex)) return NOTHING;	// verify index
	
	int left = 2 * nodeIndex;
	return isValidIndex(heap, left) ? left : NOTHING;
}

/* Returns the index of the right child of a node at index 'nodeIndex' in
 * minheap 'heap', if such exists.  Returns NOTHING if there is no such right
 * child.
 */
int rightIdx(MinHeap* heap, int nodeIndex) {
	if (!isValidIndex(heap, nodeIndex)) return NOTHING;	// verify index
	
	int right = (2 * nodeIndex) + 1;
	return isValidIndex(heap, right) ? right : NOTHING;
}

/* Returns the index of the parent of a node at index 'nodeIndex' in minheap
 * 'heap', if such exists.  Returns NOTHING if there is no such parent.
 */
int parentIdx(MinHeap* heap, int nodeIndex) {
	// Verify index and check if it is the root
	return (nodeIndex == 1 || !isValidIndex(heap, nodeIndex)) ? NOTHING : nodeIndex / 2;
}

/* Swaps contents of heap->arr[index1] and heap->arr[index2] if both 'index1'
 * and 'index2' are valid indices for minheap 'heap'. Has no effect
 * otherwise.
 */
void swap(MinHeap* heap, int index1, int index2) {
	if (isValidIndex(heap, index1) && isValidIndex(heap, index2)) {
		HeapNode temp = nodeAt(heap, index1);
		
		// Update indices in indexMap
		heap->indexMap[idAt(heap, index1)] = index2;
		heap->indexMap[idAt(heap, index2)] = index1;
		
		// Swap nodes in arr
		heap->arr[index1] = nodeAt(heap, index2);
		heap->arr[index2] = temp;
	}
}

/* Bubbles up the element newly inserted into minheap 'heap' at index
 * 'nodeIndex', if 'nodeIndex' is a valid index for heap. Has no effect
 * otherwise.
 */
void bubbleUp(MinHeap* heap, int nodeIndex) {
	if (heap == NULL || heap->size == 0) return;
	
	int parent = parentIdx(heap, nodeIndex);
	if (parent == NOTHING) return;
	
	int id = idAt(heap, nodeIndex);
	if (priorityAt(heap, parent) > priorityAt(heap, nodeIndex)) {
		swap(heap, parent, nodeIndex);
		bubbleUp(heap, indexOf(heap, id));
	}
}

/* Bubbles down the element newly inserted into minheap 'heap' at the root,
 * if it exists. Has no effect otherwise.
 */
void bubbleDown(MinHeap* heap) {
	if (heap != NULL && heap->size > 0) {
		int id = idAt(heap, ROOT_INDEX);
		int left = leftIdx(heap, ROOT_INDEX);
		int right = rightIdx(heap, ROOT_INDEX);
		
		// stop condition: no children or satisfied order of priorities
		// (left will exist if right exists, by nature of being nearly-complete)
		while (left != NOTHING) {
			int priority = priorityAt(heap, indexOf(heap, id));
			int priorityLeft = priorityAt(heap, left);
			
			// 1 child
			if (right == NOTHING) {
				if (priority > priorityLeft) 
					swap(heap, indexOf(heap, id), left);
				else break;
			}
			
			// 2 children
			else {
				int priorityRight = priorityAt(heap, right);
				if (priority > priorityLeft && priority > priorityRight) {
					if (priorityLeft < priorityRight) swap(heap, indexOf(heap, id), left);
					else swap(heap, indexOf(heap, id), right);
				}
				else if (priority > priorityLeft) swap(heap, indexOf(heap, id), left);
				else if (priority > priorityRight) swap(heap, indexOf(heap, id), right);
				else break;
			}
			
			left = leftIdx(heap, indexOf(heap, id));
			right = rightIdx(heap, indexOf(heap, id));
		}
	}
}

/*********************************************************************
 * Required functions
 ********************************************************************/

/* Returns the node with minimum priority in minheap 'heap'.
 * Precondition: heap is non-empty
 */
HeapNode getMin(MinHeap* heap) {
	return heap->arr[ROOT_INDEX];	// since heap non-empty by precond.
}

/* Removes and returns the node with minimum priority in minheap 'heap'.
 * Precondition: heap is non-empty
 */
HeapNode extractMin(MinHeap* heap) {
	// Swap root and bottom rightmost node
	swap(heap, ROOT_INDEX, heap->size);
	
	// Save and remove bottom rightmost node
	HeapNode save = nodeAt(heap, heap->size);
	heap->indexMap[idAt(heap, heap->size)] = 0;
	heap->arr[heap->size].priority = 0;
	heap->arr[heap->size].id = 0;
	heap->size--;
	
	// Bubble down newly swapped root node
	bubbleDown(heap);
	
	return save;
}

/* Inserts a new node with priority 'priority' and ID 'id' into minheap 'heap'.
 * Precondition: 'id' is unique within this minheap
 *               0 <= 'id' < heap->capacity
 *               heap->size < heap->capacity
 */
void insert(MinHeap* heap, int priority, int id) {
	HeapNode newNode;
	newNode.priority = priority;
	newNode.id = id;
	heap->arr[heap->size + 1] = newNode;		// insert into arr
	heap->indexMap[id] = heap->size + 1;		// insert into indexMap
	heap->size++;						// increment heap size
	
	// Bubble up newly inserted node
	bubbleUp(heap, heap->size);
}

/* Returns priority of the node with ID 'id' in 'heap'.
 * Precondition: 'id' is a valid node ID in 'heap'.
 */
int getPriority(MinHeap* heap, int id) {
	return priorityAt(heap, indexOf(heap, id));
}

/* Sets priority of node with ID 'id' in minheap 'heap' to 'newPriority', if
 * such a node exists in 'heap' and its priority is larger than
 * 'newPriority', and returns True. Has no effect and returns False, otherwise.
 * Note: this function bubbles up the node until the heap property is restored.
 */
bool decreasePriority(MinHeap* heap, int id, int newPriority) {
	HeapNode node = nodeAt(heap, indexOf(heap, id));
	
	if (node.priority <= newPriority) return false;
	
	heap->arr[indexOf(heap, id)].priority = newPriority;
	bubbleUp(heap, indexOf(heap, id));
	return true;
}

/* Returns a newly created empty minheap with initial capacity 'capacity'.
 * Precondition: capacity >= 0
 */
MinHeap* newHeap(int capacity) {
	MinHeap *new = malloc(sizeof(MinHeap));
	new->size = 0;
	new->capacity = capacity;
	new->arr = malloc(sizeof(HeapNode) * (capacity + 1));	// allocate for capacity and empty index 0
	new->indexMap = malloc(sizeof(int) * capacity);
	
	return new;
}

/* Frees all memory allocated for minheap 'heap'.
 */
void deleteHeap(MinHeap* heap) {
	if (heap == NULL) return;
	free(heap->arr);
	free(heap->indexMap);
	free(heap);
}

/*********************************************************************
 ** Helper function provided
 *********************************************************************/
void printHeap(MinHeap* heap) {
  printf("MinHeap with size: %d\n\tcapacity: %d\n\n", heap->size,
         heap->capacity);
  printf("index: priority [ID]\t ID: index\n");
  for (int i = 0; i < heap->capacity; i++)
    printf("%d: %d [%d]\t\t%d: %d\n", i, priorityAt(heap, i), idAt(heap, i), i,
           indexOf(heap, i));
  printf("%d: %d [%d]\t\t\n", heap->capacity, priorityAt(heap, heap->capacity),
         idAt(heap, heap->capacity));
  printf("\n\n");
}
