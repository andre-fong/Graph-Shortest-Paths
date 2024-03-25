/*
 * Graph algorithms.
 *
 * Author (of starter code): A. Tafliovich.
 */

#include <limits.h>

#include "graph.h"
#include "minheap.h"

#define NOTHING -1
#define DEBUG 0

// TODO: Try running algo w/ graph of ZERO vertices (0)

typedef struct records {
  int numVertices;    // total number of vertices in the graph
                      // vertex IDs are 0, 1, ..., numVertices-1
  MinHeap* heap;      // priority queue
  bool* finished;     // finished[id] is true iff vertex id is finished
                      //   i.e. no longer in the PQ
  int* predecessors;  // predecessors[id] is the predecessor of vertex id
  Edge* tree;         // keeps edges for the resulting tree
  int numTreeEdges;   // current number of edges in mst
} Records;

/*************************************************************************
 ** Suggested helper functions -- part of starter code
 *************************************************************************/
 
/* Returns true iff 'heap' is NULL or is empty. */
bool isEmpty(MinHeap* heap) {
	return (heap == NULL || heap->size == 0);
}

/* Creates, populates, and returns a MinHeap to be used by Prim's and
 * Dijkstra's algorithms on Graph 'graph' starting from vertex with ID
 * 'startVertex'.
 * Precondition: 'startVertex' is valid in 'graph'
 */
MinHeap* initHeap(Graph* graph, int startVertex) {
	MinHeap *new = newHeap(graph->numVertices);
	
	for (int i = 0; i < new->capacity; i++) {
		insert(new, INT_MAX, i);
	}
	
	decreasePriority(new, startVertex, 0);
	
	return new;
}

/* Creates, populates, and returns all records needed to run Prim's and
 * Dijkstra's algorithms on Graph 'graph' starting from vertex with ID
 * 'startVertex'.
 * Precondition: 'startVertex' is valid in 'graph'
 */
Records* initRecords(Graph* graph, int startVertex) {
	Records *newRecords = malloc(sizeof(Records));
	
	// Check malloc was successful
	if (newRecords == NULL) {
		fprintf(stderr, "Could not allocate memory for new record\n");
		exit(1);
	}
	newRecords->numVertices = graph->numVertices;
	
	MinHeap *heap = initHeap(graph, startVertex);
	newRecords->heap = heap;
	
	bool *finished = malloc(newRecords->numVertices * sizeof(bool));
	for (int i = 0; i < newRecords->numVertices; i++) {
		finished[i] = false;
	}
	newRecords->finished = finished;
	
	int *pred = malloc(newRecords->numVertices * sizeof(int));
	for (int i = 0; i < newRecords->numVertices; i++) {
		pred[i] = NOTHING;
	}
	pred[startVertex] = startVertex;	// set predecessor of start ind. to itself
	newRecords->predecessors = pred;
	
	Edge *tree = malloc(newRecords->numVertices * sizeof(Edge));
	newRecords->tree = tree;
	
	newRecords->numTreeEdges = 0;
	
	return newRecords;
}

/* Add a new edge to records at index ind. */
void addTreeEdge(Records* records, int ind, int fromVertex, int toVertex,
                 int weight) {
	Edge edge;
	edge.fromVertex = fromVertex;
	edge.toVertex = toVertex;
	edge.weight = weight;
	records->tree[ind] = edge;
	records->numTreeEdges++;	// to indicate we've added an edge to the tree
}

/* Creates and returns a path from 'vertex' to 'startVertex' from edges
 * in the distance tree 'distTree'.
 */
EdgeList* makePath(Edge* distTree, int vertex, int startVertex) {
	// Assuming correctness of Dijkstra's, distTree should be s.t.:
	// distTree[u.id] = an edge (u, v)
	
	// Base case: vertex == startVertex
	if (vertex == startVertex) return NULL;
	
	// Recursive case:
	int nextVertex = distTree[vertex].toVertex;
	EdgeList *next = makePath(distTree, nextVertex, startVertex);
	
	// Current weight is distance of current vertex in distTree - distance of next vertex in distTree
	int weight = distTree[vertex].weight - distTree[nextVertex].weight;
	Edge *edge = newEdge(vertex, nextVertex, weight);
	
	return newEdgeList(edge, next);
}

/* Frees all memory previously allocated for records. */
void deleteRecords(Records *records) {
	deleteHeap(records->heap);
	free(records->finished);
	free(records->predecessors);
	free(records->tree);
	free(records);
}

/* Allocates space for a new Edge array from records.
 * Prereqs: * should be called before deleteRecords(records).
 *          * records should follow a complete run of Prim's or Dijkstra's
 */
Edge *extractTreeFromRecords(Records *records) {
	Edge *tree = malloc(records->numTreeEdges * sizeof(Edge));
	for (int i = 0; i < records->numTreeEdges; i++) {
		tree[i] = records->tree[i];
	}
	return tree;
}
	
/* Get other vertex from edge in v.id's adjacency list */
int getIdFromEdge(Edge *edge, int id) {
	return (edge->fromVertex == id) ? edge->toVertex : edge->fromVertex;
}

/*************************************************************************
 ** Required functions
 *************************************************************************/

/* Runs Prim's algorithm on Graph 'graph' starting from vertex with ID
 * 'startVertex', and return the resulting MST: an array of Edges.
 * Returns NULL if 'startVertex' is not valid in 'graph'.
 * Precondition: 'graph' is connected.
 */
Edge* getMSTprim(Graph* graph, int startVertex) {
	// Check startVertex is valid in graph
	if (startVertex >= graph->numVertices || startVertex < 0) return NULL;
	
	// Initiate records for Prim's
	Records *records = initRecords(graph, startVertex);
	
	while (!isEmpty(records->heap)) {
		// Extract min and add to T
		HeapNode min = extractMin(records->heap);
		
		records->finished[min.id] = true;
		// Don't write (startVertex, startVertex, 0) to Edge array
		if (min.id != startVertex) {
			addTreeEdge(records, records->numTreeEdges, min.id, records->predecessors[min.id], min.priority);
		}
		
		// Iterate through min's adjacency list
		Vertex *v = graph->vertices[min.id];
		EdgeList *tr = v->adjList;
		
		while (tr != NULL) {
			int weight = tr->edge->weight;
			int nextVertexId = getIdFromEdge(tr->edge, min.id);
			
			// If next vertex in heap and current weight < priority of next vertex, decrease priority of next vertex
			if (!records->finished[nextVertexId] && weight < getPriority(records->heap, nextVertexId)) {
				decreasePriority(records->heap, nextVertexId, weight);
				records->predecessors[nextVertexId] = min.id;
			}
			
			tr = tr->next;
		}
	}
	
	Edge *tree = extractTreeFromRecords(records);
	deleteRecords(records);
	return tree;
}

/* Runs Dijkstra's algorithm on Graph 'graph' starting from vertex with ID
 * 'startVertex', and return the resulting distance tree: an array of edges.
 * Returns NULL if 'startVertex' is not valid in 'graph'.
 * Precondition: 'graph' is connected.
 */
Edge* getDistanceTreeDijkstra(Graph* graph, int startVertex) {
	// Check startVertex is valid in graph
	if (startVertex >= graph->numVertices || startVertex < 0) return NULL;
	
	// Initiate records for Dijkstra's
	Records *records = initRecords(graph, startVertex);
	
	while (!isEmpty(records->heap)) {
		// Extract min and add to T
		HeapNode min = extractMin(records->heap);
		records->finished[min.id] = true;
		// Write (u, u.pred) in index u for simplicity in later algo.
		addTreeEdge(records, min.id, min.id, records->predecessors[min.id], min.priority);
		
		// Iterate through min's adjacency list
		Vertex *v = graph->vertices[min.id];
		EdgeList *tr = v->adjList;
		
		while (tr != NULL) {
			int distance = min.priority + tr->edge->weight;
			int nextVertexId = getIdFromEdge(tr->edge, min.id);
			
			// If next vertex in heap and current dist < priority of next vertex, decrease priority of next vertex
			if (!records->finished[nextVertexId] && distance < getPriority(records->heap, nextVertexId)) {
				decreasePriority(records->heap, nextVertexId, distance);
				records->predecessors[nextVertexId] = min.id;
			}
			
			tr = tr->next;
		}
	}
	
	Edge *tree = extractTreeFromRecords(records);
	deleteRecords(records);
	return tree;
}

/* Creates and returns an array 'paths' of shortest paths from every vertex
 * in the graph to vertex 'startVertex', based on the information in the
 * distance tree 'distTree' produced by Dijkstra's algorithm on a graph with
 * 'numVertices' vertices and with the start vertex 'startVertex'.  paths[id]
 * is the list of edges of the form
 *   [(id -- id_1, w_0), (id_1 -- id_2, w_1), ..., (id_n -- start, w_n)]
 *   where w_0 + w_1 + ... + w_n = distance(id)
 * Returns NULL if 'startVertex' is not valid in 'distTree'.
 */
EdgeList** getShortestPaths(Edge* distTree, int numVertices, int startVertex) {
	// Check startVertex is valid in graph
	if (startVertex >= numVertices || startVertex < 0) return NULL;
	
	EdgeList **paths = malloc(numVertices * sizeof(EdgeList *));
	
	// Check malloc was successful
	if (paths == NULL) {
		fprintf(stderr, "Could not allocate memory for paths list\n");
		exit(1);
	}
	
	for (int i = 0; i < numVertices; i++) {
		paths[i] = makePath(distTree, i, startVertex);
	}
	
	return paths;
}

/*************************************************************************
 ** Provided helper functions -- part of starter code to help you debug!
 *************************************************************************/
void printRecords(Records* records) {
  if (records == NULL) return;

  int numVertices = records->numVertices;
  printf("Reporting on algorithm's records on %d vertices...\n", numVertices);

  printf("The PQ is:\n");
  printHeap(records->heap);

  printf("The finished array is:\n");
  for (int i = 0; i < numVertices; i++)
    printf("\t%d: %d\n", i, records->finished[i]);

  printf("The predecessors array is:\n");
  for (int i = 0; i < numVertices; i++)
    printf("\t%d: %d\n", i, records->predecessors[i]);

  printf("The TREE edges are:\n");
  for (int i = 0; i < records->numTreeEdges; i++) printEdge(&records->tree[i]);

  printf("... done.\n");
}
