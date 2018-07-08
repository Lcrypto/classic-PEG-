#ifndef BIGGIRTH
#define BIGGIRTH

#include <stdlib.h>
#include <iostream> // C++ I/O library header
#include "Random.h"

using namespace std;

class NodesInGraph{
 public:
  int numOfConnectionParityBit;
  int *connectionParityBit;
  int numOfConnectionSymbolBit;
  int *connectionSymbolBit;
  int maxDegParity;

  NodesInGraph(void);
  ~NodesInGraph(void);
  void setNumOfConnectionSymbolBit(int deg);
  void initConnectionParityBit(void);
  void initConnectionParityBit(int deg);
};

class BigGirth {
 public:
  int M, N;
  int K;
  int EXPAND_DEPTH;
  char *filename;
  int *(*H);

  int *localGirth;
  
  NodesInGraph *nodesInGraph;
  Random *myrandom;

  BigGirth(int m, int n, int quickEnc, int *symbolDegSequence, int *checkDegSequence, char *filename, int sglConcent, int tgtGirth, int verbose);
  BigGirth(void);

  void writeToFile_Hcompressed(void);
  void writeToFile_Hmatrix(void);
  void writeToFile(void);

  void loadH(void);

  ~BigGirth(void);

 private:
  int selectParityConnect(int kthSymbol, int mthConnection, int & cycle, int quickEnc);
  void updateConnection(int kthSymbol);

};

#endif









