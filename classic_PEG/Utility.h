#ifndef UTILITY
#define UTILITY

#include <valarray>


class MatrixInt2D {
public:
  typedef int ValueType;

public:
  MatrixInt2D( std::size_t rowCount, std::size_t colCount ) 
    : m_rowCount( rowCount ), m_colCount( colCount ), m_data( 0, rowCount * colCount ) { }

  std::size_t GetColCount() const { return m_colCount; }
  std::size_t GetRowCount() const { return m_rowCount; }
  const std::valarray< ValueType > & GetData() const { return m_data; }

  ValueType* operator[] ( std::size_t row ) { return &m_data[ row * m_colCount ]; } 

private:
  std::valarray< ValueType > m_data;
  std::size_t m_rowCount;
  std::size_t m_colCount;
};

#endif
