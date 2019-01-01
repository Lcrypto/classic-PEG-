#ifndef BIGGIRTH
#define BIGGIRTH

class NodesInGraph;
class Random;

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









