#ifndef __itkFlatStructuringElement_h
#define __itkFlatStructuringElement_h

#include "itkNeighborhood.h"
#include "itkSize.h"
#include "itkOffset.h"
#include <vector>
#include "itkVector.h"

namespace itk {

/** \class FlatStructuringElement
* \brief A class to support a variety of flat structuring elements, 
* including versions created by decomposition of lines.
**/

template<unsigned int VDimension>
class ITK_EXPORT FlatStructuringElement : public Neighborhood <bool, VDimension>
{
public:
  /** Standard class typedefs. */
  typedef FlatStructuringElement Self;
  //typedef Neighborhood<bool, VDimension> Superclass;

  /** External support for dimensionality. */
  itkStaticConstMacro(NeighborhoodDimension, unsigned int, VDimension);

  /** Radius typedef support. */
  typedef Size<VDimension> RadiusType;

  typedef Vector<float, VDimension> LType;
  typedef std::vector<LType> DecompType;

  /** Default destructor. */
  virtual ~FlatStructuringElement() {}

  /** Default consructor. */
  FlatStructuringElement() {m_Decomposition=false;}

  /** Various constructors */

  static FlatStructuringElement Box(RadiusType radius);
  
  // lines is the number of elements in the decomposition
  static FlatStructuringElement Poly(RadiusType radius, unsigned lines);

  bool Decomposable()
  {
    return m_Decomposition;
  }

  void PrintSelf(std::ostream &os, Indent indent) const;

  const DecompType &GetDecomp()
  {
    return(m_Lines);
  }

private:
  bool m_Decomposition;

  DecompType m_Lines;
  
  // dispatch between 2D and 3D
  struct DispatchBase {};
  template<unsigned int VDimension2>
  struct Dispatch : DispatchBase {};

  virtual FlatStructuringElement PolySub(const Dispatch<2> &, 
					 RadiusType radius, 
					 unsigned lines) const;

  virtual FlatStructuringElement PolySub(const Dispatch<3> &, 
					 RadiusType radius, 
					 unsigned lines) const;

  virtual FlatStructuringElement PolySub(const DispatchBase &, 
					 RadiusType radius, 
					 unsigned lines) const;


  bool checkParallel(LType NewVec, DecompType Lines);

  typedef struct {
    LType P1, P2, P3;
  } FacetType;

};
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkFlatStructuringElement.txx"
#endif



#endif
