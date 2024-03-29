#ifndef __itkAnchorOpenCloseLine_h
#define __itkAnchorOpenCloseLine_h

#include "itkAnchorHistogram.h"

//#define RAWHIST

namespace itk {



/** 
 * \class AnchorOpenCloseLine
 * \brief class to implement openings and closings using anchor
 * methods. This is the base class that must be instantiated with
 * appropriate definitions of greater, less and so on

**/
template<class TInputPix, class THistogramCompare,
	 class TFunction1, class TFunction2>
class ITK_EXPORT AnchorOpenCloseLine
{
public:
  /** Some convenient typedefs. */
  typedef TInputPix InputImagePixelType;
  AnchorOpenCloseLine();
  ~AnchorOpenCloseLine() {delete m_Histo;};
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Single-threaded version of GenerateData.  This filter delegates
   * to GrayscaleGeodesicErodeImageFilter. */
  void doLine(InputImagePixelType * buffer, unsigned bufflength);

  void SetSize(unsigned int size)
  {
    m_Size = size;
  }

private:
  unsigned int m_Size;
  TFunction1 m_TF1;
  TFunction2 m_TF2;

  typedef MorphologyHistogram<InputImagePixelType> Histogram;
  typedef MorphologyHistogramVec<InputImagePixelType,THistogramCompare> VHistogram;
  typedef MorphologyHistogramMap<InputImagePixelType,THistogramCompare> MHistogram;

  bool startLine(InputImagePixelType * buffer,
		 InputImagePixelType &Extreme,
		 Histogram &histo,
		 unsigned &outLeftP,
		 unsigned &outRightP);

  bool finishLine(InputImagePixelType * buffer,
		  InputImagePixelType &Extreme,
		  unsigned &outLeftP,
		  unsigned &outRightP);

  bool useVectorBasedHistogram()
  {
    // bool, short and char are acceptable for vector based algorithm: they do not require
    // too much memory. Other types are not usable with that algorithm
    return typeid(InputImagePixelType) == typeid(unsigned char)
        || typeid(InputImagePixelType) == typeid(signed char)
        || typeid(InputImagePixelType) == typeid(unsigned short)
        || typeid(InputImagePixelType) == typeid(signed short)
        || typeid(InputImagePixelType) == typeid(bool);
  }

  Histogram * m_Histo;

} ; // end of class


} // end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkAnchorOpenCloseLine.txx"
#endif

#endif

