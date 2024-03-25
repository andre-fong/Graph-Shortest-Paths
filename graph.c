/*
 * Our graph implementation.
 *
 * Author: A. Tafliovich.
 */

#include "graph.h"

/*********************************************************************
 ** Helper function provided in the starter code
 *********************************************************************/

void printEdge(Edge* edge) {
  if (edge == NULL)
    printf("NULL");
  else
    printf("(%d -- %d, %d)", edge->fromVertex, edge->toVertex, edge->weight);
}

void printEdgeList(EdgeList* head) {
  while (head != NULL) {
    printEdge(head->edge);
    printf(" --> ");
    head = head->next;
  }
  printf("NULL");
}

void printVertex(Vertex* vertex) {
  if (vertex == NULL) {
    printf("NULL");
  } else {
    printf("%d: ", vertex->id);
    printEdgeList(vertex->adjList);
  }
}

void printGraph(Graph* graph) {
  if (graph == NULL) {
    printf("NULL");
    return;
  }
  printf("Number of vertices: %d. Number of edges: %d.\n\n", graph->numVertices,
         graph->numEdges);

  for (int i = 0; i < graph->numVertices; i++) {
    printVertex(graph->vertices[i]);
    printf("\n");
  }
  printf("\n");
}

/*********************************************************************
 ** Required functions
 *********************************************************************/

/* Returns a newly created Edge from vertex with ID 'fromVertex' to vertex
 * with ID 'toVertex', with weight 'weight'.
 */
Edge* newEdge(int fromVertex, int toVertex, int weight) {
	Edge *new = malloc(sizeof(Edge));
	
	// Check malloc successful
	if (new == NULL) {
		fprintf(stderr, "Could not allocate memory for new edge\n");
		exit(1);
	}
	new->fromVertex = fromVertex;
	new->toVertex = toVertex;
	new->weight = weight;
	
	return new;
}

/* Returns a newly created EdgeList containing 'edge' and pointing to the next
 * EdgeList node 'next'.
 */
EdgeList* newEdgeList(Edge* edge, EdgeList* next) {
	EdgeList *newList = malloc(sizeof(EdgeList));
	
	// Check malloc successful
	if (newList == NULL) {
		fprintf(stderr, "Could not allocate memory for new edge list\n");
		exit(1);
	}
	newList->edge = edge;
	newList->next = next;
	
	return newList;
}

/* Returns a newly created Vertex with ID 'id', value 'value', and adjacency
 * list 'adjList'.
 * Precondition: 'id' is valid for this vertex
 */
Vertex* newVertex(int id, void* value, EdgeList* adjList) {
	Vertex *new = malloc(sizeof(Vertex));
	
	// Check malloc successful
	if (new == NULL) {
		fprintf(stderr, "Could not allocate memory for new vertex\n");
		exit(1);
	}
	new->id = id;
	new->value = value;
	new->adjList = adjList;
	
	return new;
}

/* Returns a newly created Graph with space for 'numVertices' vertices.
 * Precondition: numVertices >= 0
 */
Graph* newGraph(int numVertices) {
	Graph *new = malloc(sizeof(Graph));
	
	// Check malloc successful
	if (new == NULL) {
		fprintf(stderr, "Could not allocate memory for new graph\n");
		exit(1);
	}
	new->numVertices = numVertices;
	new->numEdges = 0;
	
	Vertex **vertices;
	// Necessary since calloc(0) may return NULL on some implementations
	if (numVertices == 0) {
		vertices = NULL;
	}
	else {
		vertices = (Vertex **)calloc(numVertices, sizeof(Vertex *));
		if (vertices == NULL) {
			fprintf(stderr, "Could not allocate memory for graph vertices list\n");
			exit(1);
		}
	}
	
	new->vertices = vertices;
	
	return new;
}

/* Frees memory allocated for EdgeList starting at 'head'.
 */
void deleteEdgeList(EdgeList* head) {
	EdgeList *pre = head;
	EdgeList *tr = NULL;
	while (pre != NULL) {
		tr = pre->next;
		free(pre->edge);	// free edge (edge lists have unique edges, ie. (a, b, 2) in a's list =/= (b, a, 2) in b's list)
		free(pre);			// then free edge list node
		pre = tr;
	}
}

/* Frees memory allocated for 'vertex' including its adjacency list.
 */
void deleteVertex(Vertex* vertex) {
	if (vertex == NULL) return;
	deleteEdgeList(vertex->adjList);
	free(vertex);
}

/* Frees memory allocated for 'graph'.
 */
void deleteGraph(Graph* graph) {
	for (int i = 0; i < graph->numVertices; i++) {
		deleteVertex(graph->vertices[i]);
	}
	free(graph->vertices);
	free(graph);
}
