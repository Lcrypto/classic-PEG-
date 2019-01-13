#include "CyclesOfGraph.h"
#include "Utility.h"

#include <cstdlib>
#include <iostream>
#include <valarray>


using namespace std;

class NodesOfGraph{
public:
  NodesOfGraph();
  ~NodesOfGraph();

  int GetNumOfParityConnections() const { return int( parityConnections.size() ); }
  int GetNumOfSymbolConnections() const { return int( symbolConnections.size() ); }
  int GetNumOfSymbolMapping() const { return int( symbolMapping.size() ); }

  std::valarray< int > parityConnections;
  std::valarray< int > symbolConnections;
  std::valarray< int > symbolMapping;

  void setParityConnections(int num, int *values);
  void setSymbolConnections(int num, int *values);
  void setSymbolMapping(int num, int *values);
};


NodesOfGraph::NodesOfGraph() = default;

NodesOfGraph::~NodesOfGraph() = default;

void NodesOfGraph::setParityConnections(int num, int *values) 
{
  parityConnections.resize( num );
  std::memcpy( &parityConnections[0], values, num );
}

void NodesOfGraph::setSymbolConnections(int num, int *values) 
{
  symbolConnections.resize( num );
  std::memcpy( &symbolConnections[0], values, num );
}
void NodesOfGraph::setSymbolMapping(int num, int *values) 
{
  symbolMapping.resize( num );
  std::memcpy( &symbolMapping[0], values, num );
}


CyclesOfGraph::CyclesOfGraph(int mm, int n, std::shared_ptr< MatrixInt2D > h )
{
  int i, j, k, m, index;
  M=mm;
  N=n;
  H=h;

  m_tmp.resize( N );
  m_med.resize( N );
  m_tmpCycles.resize( N );
  m_cyclesTable.resize( N );
  m_nodesOfGraph.resize( N );

  //cout<<M<<" "<<N<<endl;
  /*
  for(i=0;i<M;i++){
    for(j=0;j<N;j++)
      cout<<(*H)[i][j]<<" ";
    cout<<endl;
  }
  */
  for(i=0;i<N;i++){
    index=0;
    for(j=0;j<M;j++){
      if((*H)[j][i]==1){
    m_tmp[index]=j;
    index++;
      }
    }
    m_nodesOfGraph[i].setSymbolConnections(index, &m_tmp[0]);
  }
  for(i=0;i<M;i++){
    index=0;
    for(j=0;j<N;j++){
      if((*H)[i][j]==1){
    m_tmp[index]=j;
    index++;
      }
    }
    m_nodesOfGraph[i].setParityConnections(index, &m_tmp[0]);
  }
  for(i=0;i<N;i++){
    index=0;
    for(j=0;j<m_nodesOfGraph[i].GetNumOfSymbolConnections();j++){
      for(k=0;k<m_nodesOfGraph[m_nodesOfGraph[i].symbolConnections[j]].GetNumOfParityConnections();k++){
    int t=0;
    for(m=0;m<index;m++){
      if(m_nodesOfGraph[m_nodesOfGraph[i].symbolConnections[j]].parityConnections[k]==m_tmp[m]){
        t=1; break;
      }
    }
    if(m_nodesOfGraph[m_nodesOfGraph[i].symbolConnections[j]].parityConnections[k]==i) t=1;
    if(t==0) {
      m_tmp[index]=m_nodesOfGraph[m_nodesOfGraph[i].symbolConnections[j]].parityConnections[k];
      index++;
    }
      }
    }
    m_nodesOfGraph[i].setSymbolMapping(index, &m_tmp[0]);
  }
}

CyclesOfGraph::~CyclesOfGraph() = default;
  
void CyclesOfGraph::getCyclesTable(void) {
  int i, j, k, m, n, t, imed;
  for(i=0;i<N;i++){
    //special handlement for nodes having only one or zero connections
    if(m_nodesOfGraph[i].GetNumOfSymbolConnections()<=1) {
      m_cyclesTable[i]=2*N; 
      continue;
    }
    for(j=0;j<m_nodesOfGraph[i].GetNumOfSymbolConnections()-1;j++){ //-1 because the graph is undirected
      for(k=0;k<m_nodesOfGraph[m_nodesOfGraph[i].symbolConnections[j]].GetNumOfParityConnections();k++){
    m_tmp[k]=m_nodesOfGraph[m_nodesOfGraph[i].symbolConnections[j]].parityConnections[k];
    //cout<<m_tmp[k]<<" ";
      }
      //cout<<endl;
      int cycles=2;
      int index=m_nodesOfGraph[m_nodesOfGraph[i].symbolConnections[j]].GetNumOfParityConnections();
    LOOP:
      imed=0;
      for(k=0;k<index;k++){
    if(m_tmp[k]==i) continue;
    //cout<<"k="<<k<<" "<<m_tmp[k]<<endl;
    for(m=0;m<m_nodesOfGraph[m_tmp[k]].GetNumOfSymbolConnections();m++){
      for(n=0;n<m_nodesOfGraph[i].GetNumOfSymbolConnections();n++){
        if((n!=j)&&(m_nodesOfGraph[m_tmp[k]].symbolConnections[m]==m_nodesOfGraph[i].symbolConnections[n])){
          cycles+=2;
          goto OUTLOOP;
        }
      }
    }
    for(m=0;m<m_nodesOfGraph[m_tmp[k]].GetNumOfSymbolMapping();m++){
      t=0;
      for(int l=0;l<imed;l++) {
        if(m_nodesOfGraph[m_tmp[k]].symbolMapping[m]==m_med[l]){
          t=1; break;
        }
      }
      if(t==0){
        m_med[imed]=m_nodesOfGraph[m_tmp[k]].symbolMapping[m];
        //cout<<m_med[imed]<<endl;
        imed++;
      }
    }
      }
      index=imed;//cout<<index<<" "<<endl;
      for(k=0;k<index;k++) {
    m_tmp[k]=m_med[k];//cout<<m_tmp[k]<<" ";
      }
      //cout<<"j="<<j<<endl;
      cycles+=2;
      if(cycles>=2*N) //dead lock 
    goto OUTLOOP;
      else
    goto LOOP;
    OUTLOOP:
      m_tmpCycles[j]=cycles;
    }
    //for(j=0;j<m_nodesOfGraph[i].GetNumOfSymbolConnections()-1;j++) cout<<m_tmpCycles[j]<<" ";
    //cout<<endl;
    m_cyclesTable[i]=m_tmpCycles[0];
    for(j=1;j<m_nodesOfGraph[i].GetNumOfSymbolConnections()-1;j++){
      if(m_cyclesTable[i]>m_tmpCycles[j])
    m_cyclesTable[i]=m_tmpCycles[j];
    }
    //OUTPUT cycles per symbol node
    //cout<<"i="<<i<<" "<<m_cyclesTable[i]<<endl;
  }
}
    
int CyclesOfGraph::FindGirth() const 
{
  int girth = 2 * N;
  
  for( int i = 0; i < N; ++i ) 
    if( girth > m_cyclesTable[i] ) 
      girth = m_cyclesTable[i];
  
  return girth;
}

void CyclesOfGraph::printCyclesTable(void){
  int i, temp[20];
  /*
  for(i=0;i<N;i++)
    cout<<m_cyclesTable[i]<<" ";
  cout<<endl;
  */
  for(i=0;i<20;i++) temp[i]=0;
  for(i=0;i<N;i++){
    if(m_cyclesTable[i]==4) temp[0]++;
    else if(m_cyclesTable[i]==6) temp[1]++;    
    else if(m_cyclesTable[i]==8) temp[2]++;
    else if(m_cyclesTable[i]==10) temp[3]++;
    else if(m_cyclesTable[i]==12) temp[4]++;
    else if(m_cyclesTable[i]==14) temp[5]++;
    else if(m_cyclesTable[i]==16) temp[6]++;
    else if(m_cyclesTable[i]==18) temp[7]++;
    else if(m_cyclesTable[i]==20) temp[8]++;
    else if(m_cyclesTable[i]==22) temp[9]++;
    else if(m_cyclesTable[i]==24) temp[10]++;
    else if(m_cyclesTable[i]==26) temp[11]++;
    else if(m_cyclesTable[i]==28) temp[12]++;
    else if(m_cyclesTable[i]==30) temp[13]++;
    else {
      cout<<"Wrong cycles calculation   "<<m_cyclesTable[i]<<endl;
      exit(-1);
    }
  }  
  cout<<endl;
  cout<<"Num of Nodes with local girth   4:  "<< temp[0]<<endl;
  cout<<"Num of Nodes with local girth   6:  "<< temp[1]<<endl;
  cout<<"Num of Nodes with local girth   8:  "<< temp[2]<<endl;
  cout<<"Num of Nodes with local girth   10:  "<< temp[3]<<endl;
  cout<<"Num of Nodes with local girth   12:  "<< temp[4]<<endl;
  cout<<"Num of Nodes with local girth   14:  "<< temp[5]<<endl;
  cout<<"Num of Nodes with local girth   16:  "<< temp[6]<<endl;
  cout<<"Num of Nodes with local girth   18:  "<< temp[7]<<endl;
  cout<<"Num of Nodes with local girth   20:  "<< temp[8]<<endl;
  cout<<"Num of Nodes with local girth   22:  "<< temp[9]<<endl;
  cout<<"Num of Nodes with local girth   24:  "<< temp[10]<<endl;
  cout<<"Num of Nodes with local girth   26:  "<< temp[11]<<endl;
  cout<<"Num of Nodes with local girth   28:  "<< temp[12]<<endl;
  cout<<"Num of Nodes with local girth   30:  "<< temp[13]<<endl;
}



