#ifndef __itkAnchorErodeDilateLine_h
#define __itkAnchorErodeDilateLine_h

#include "itkAnchorHistogram.h"

namespace itk {

/** 
 * \class AnchorErodeDilateLine
 * \brief class to implement erosions and dilations using anchor
 * methods. This is the base class that must be instantiated with
 * appropriate definitions of greater, less and so on

**/
template<class TInputPix, class THistogramCompare,
	 class TFunction1, class TFunction2>
class ITK_EXPORT AnchorErodeDilateLine
{
public:
  /** Some convenient typedefs. */
  typedef TInputPix InputImagePixelType;

  void doLine(InputImagePixelType * buffer, InputImagePixelType * inbuffer, 
	      unsigned bufflength);

  void SetSize(unsigned int size)
  {
    m_Size = size;
  }
  //itkGetConstReferenceMacro(Size, unsigned int);

  //itkSetMacro(Direction, unsigned int);
  //itkGetConstReferenceMacro(Direction, unsigned int);

  void PrintSelf(std::ostream &os, Indent indent) const;
  AnchorErodeDilateLine();
  ~AnchorErodeDilateLine() {delete m_Histo;};


private:
  unsigned int m_Size;
  TFunction1 m_TF1;
  TFunction2 m_TF2;

  bool m_useVec;

  typedef MorphologyHistogram<InputImagePixelType> Histogram;
  typedef MorphologyHistogramVec<InputImagePixelType,THistogramCompare> VHistogram;
  typedef MorphologyHistogramMap<InputImagePixelType,THistogramCompare> MHistogram;

  bool startLine(InputImagePixelType * buffer,
		 InputImagePixelType * inbuffer,
		 InputImagePixelType &Extreme,
		 Histogram &histo,
		 unsigned &outLeftP,
		 unsigned &outRightP,
		 unsigned &inLeftP,
		 unsigned &inRightP,
		 unsigned middle);

  bool finishLine(InputImagePixelType * buffer,
		  InputImagePixelType * inbuffer,
		  InputImagePixelType &Extreme,
		  Histogram &histo,
		  unsigned &outLeftP,
		  unsigned &outRightP,
		  unsigned &inLeftP,
		  unsigned &inRightP,
		  unsigned middle);

  bool useVectorBasedHistogram()
  {
    return(false);
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
#include "itkAnchorErodeDilateLine.txx"
#endif

#endif

