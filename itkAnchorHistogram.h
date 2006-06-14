// histogram from the moving histogram operations
#ifndef __itkAnchorHistogram_h
#define __itkAnchorHistogram_h
#include "itkNumericTraits.h"

namespace itk {

namespace Function {
template <class TInputPixel, class TCompare>
class MorphologyHistogram
{
public:
  MorphologyHistogram()
    {
    if( useVectorBasedAlgorithm() )
      { initVector(); }
    }
  ~MorphologyHistogram(){}

  // define the method required by the functor and dispatch to the specialized methods

  void Init()
  {
    if( useVectorBasedAlgorithm() )
      { 
      initVector(); 
      }
    else
      {
      m_Map.clear();
      }
  }

  inline void AddBoundary()
    {
    if( useVectorBasedAlgorithm() )
      { AddBoundaryVector(); }
    else
      { AddBoundaryMap(); }
    }

  inline void RemoveBoundary()
    {
    if( useVectorBasedAlgorithm() )
      { RemoveBoundaryVector(); }
    else
      { RemoveBoundaryMap(); }
    }

  inline void AddPixel( const TInputPixel &p )
    {
    if( useVectorBasedAlgorithm() )
      { AddPixelVector( p ); }
    else
      { AddPixelMap( p ); }
    }

  inline void RemovePixel( const TInputPixel &p )
    {
    if( useVectorBasedAlgorithm() )
      { RemovePixelVector( p ); }
    else
      { RemovePixelMap( p ); }
    }

  inline TInputPixel GetValue()
    {
    if( useVectorBasedAlgorithm() )
      { return GetValueVector(); }
    else
      { return GetValueMap(); }
    }

  inline bool useVectorBasedAlgorithm()
    {
    // bool, short and char are acceptable for vector based algorithm: they do not require
    // too much memory. Other types are not usable with that algorithm
    return typeid(TInputPixel) == typeid(unsigned char)
        || typeid(TInputPixel) == typeid(signed char)
        || typeid(TInputPixel) == typeid(unsigned short)
        || typeid(TInputPixel) == typeid(signed short)
        || typeid(TInputPixel) == typeid(bool);
    }



  //
  // the map based algorithm
  //

  typedef typename std::map< TInputPixel, unsigned long, TCompare > MapType;

  inline void AddBoundaryMap()
    { m_Map[ m_Boundary ]++; }

  inline void RemoveBoundaryMap()
    { m_Map[ m_Boundary ]--; }

  inline void AddPixelMap( const TInputPixel &p )
    { m_Map[ p ]++; }

  inline void RemovePixelMap( const TInputPixel &p )
    { m_Map[ p ]--; }

  inline TInputPixel GetValueMap()
    {
    // clean the map
    typename MapType::iterator mapIt = m_Map.begin();
    while( mapIt != m_Map.end() )
      {
      if( mapIt->second == 0 )
        { 
        // this value must be removed from the histogram
        // The value must be stored and the iterator updated before removing the value
        // or the iterator is invalidated.
        TInputPixel toErase = mapIt->first;
        mapIt++;
        m_Map.erase( toErase );
        }
      else
        {
        mapIt++;
        // don't remove all the zero value found, just remove the one before the current maximum value
        // the histogram may become quite big on real type image, but it's an important increase of performances
        break;
        }
      }

    // and return the value
    return m_Map.begin()->first;
    }

  MapType m_Map;




  //
  // the vector based algorithm
  //

  inline void initVector()
    {
    // initialize members need for the vector based algorithm
    m_Vector.resize( static_cast<int>( NumericTraits< TInputPixel >::max() - NumericTraits< TInputPixel >::NonpositiveMin() + 1 ), 0 );
    if( m_Compare( NumericTraits< TInputPixel >::max(), NumericTraits< TInputPixel >::NonpositiveMin() ) )
      {
      m_CurrentValue = NumericTraits< TInputPixel >::NonpositiveMin();
      m_Direction = -1;
      }
    else
      {
      m_CurrentValue = NumericTraits< TInputPixel >::max();
      m_Direction = 1;
      }
    }
  

  inline void AddBoundaryVector()
    { AddPixelVector( m_Boundary ); }

  inline void RemoveBoundaryVector()
    { RemovePixelVector( m_Boundary ); }

  inline void AddPixelVector( const TInputPixel &p )
    {
    m_Vector[ static_cast<int>( p - NumericTraits< TInputPixel >::NonpositiveMin() ) ]++;
    if( m_Compare( p, m_CurrentValue ) )
      { m_CurrentValue = p; }
    }

  inline void RemovePixelVector( const TInputPixel &p )
    {
    m_Vector[ static_cast<int>( p - NumericTraits< TInputPixel >::NonpositiveMin() ) ]--;
    while( m_Vector[ static_cast<int>( m_CurrentValue - NumericTraits< TInputPixel >::NonpositiveMin() ) ] == 0 )
      { m_CurrentValue += m_Direction; }
    }

  inline TInputPixel GetValueVector()
    { return m_CurrentValue; }

  std::vector<unsigned long> m_Vector;
  TInputPixel m_CurrentValue;
  TCompare m_Compare;
  signed int m_Direction;



  // accessor for boundary value

  void SetBoundary( const TInputPixel & val )
    { m_Boundary = val; }

  TInputPixel m_Boundary;
};
} // end namespace Function
} // end namespace itk
#endif
