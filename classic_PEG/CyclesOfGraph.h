#ifndef CYCLESOFGRAPH
#define CYCLESOFGRAPH

#include <valarray>
#include <vector>
#include <memory>


class MatrixInt2D;
class NodesOfGraph;


class CyclesOfGraph {
public:
  CyclesOfGraph(int mm, int n, std::shared_ptr< MatrixInt2D > h);
  ~CyclesOfGraph();

  void getCyclesTable(void);
  void printCyclesTable(void);
  int FindGirth() const;

private:
  int M;
  int N;

  std::shared_ptr< MatrixInt2D > H;

  std::valarray< int > m_tmp;
  std::valarray< int > m_med;
  std::valarray< int > m_tmpCycles;

  std::valarray< int > m_cyclesTable;
  std::vector< NodesOfGraph > m_nodesOfGraph;
};

#endif
