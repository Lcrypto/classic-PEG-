#ifndef BIGGIRTH
#define BIGGIRTH

#include <vector>
#include <memory>
#include <valarray>


class MatrixInt2D;
class NodesInGraph;
class Random;

class BigGirth {
public:
  BigGirth(int m, int n, int quickEnc, int *symbolDegSequence, int *checkDegSequence, char *filename, int sglConcent, int tgtGirth, int verbose);
  ~BigGirth();

  void writeToFile_Hcompressed(void);
  void writeToFile_Hmatrix(void);
  void writeToFile(void);

  void loadH();

private:
  int selectParityConnect(int kthSymbol, int mthConnection, int & cycle, int quickEnc);
  void updateConnection(int kthSymbol);

private:
  int M, N;
  int K;
  int EXPAND_DEPTH;

  char *filename;

  std::shared_ptr< MatrixInt2D > H;
  std::valarray< int > localGirth;
  std::vector< NodesInGraph > nodesInGraph;
  std::unique_ptr< Random > myrandom;
};

#endif









