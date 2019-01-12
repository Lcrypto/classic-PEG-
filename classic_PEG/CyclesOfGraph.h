#ifndef CYCLESOFGRAPH
#define CYCLESOFGRAPH


class NodesOfGraph;


class CyclesOfGraph {
 public:
  int M, N;
  int *(*H);
  int *cyclesTable;
  NodesOfGraph *nodesOfGraph;
  CyclesOfGraph(int mm, int n, int *(*h));
  ~CyclesOfGraph(void);
  void getCyclesTable(void);
  void printCyclesTable(void);
  int girth(void);
 private:
  int *tmp, *med, *tmpCycles;
};

#endif
