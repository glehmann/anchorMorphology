#ifndef __itkAnchorErodeDilateImageFilter_h
#define __itkAnchorErodeDilateImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkProgressReporter.h"
#include "itkAnchorHistogram.h"

//#define RAWHIST

namespace itk {



/** 
 * \class AnchorErodeDilateImageFilter
 * \brief class to implement erosions and dilations using anchor
 * methods. This is the base class that must be instantiated with
 * appropriate definitions of greater, less and so on

**/
template<class TInputImage, class TOutputImage, class THistogramCompare,
	 class TFunction1, class TFunction2>
class ITK_EXPORT AnchorErodeDilateImageFilter :
    public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef AnchorErodeDilateImageFilter Self;
  typedef ImageToImageFilter<TInputImage, TOutputImage>
  Superclass;
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Some convenient typedefs. */
  typedef TInputImage InputImageType;
  typedef TOutputImage OutputImageType;
  typedef typename InputImageType::Pointer         InputImagePointer;
  typedef typename InputImageType::ConstPointer    InputImageConstPointer;
  typedef typename InputImageType::RegionType      InputImageRegionType;
  typedef typename InputImageType::PixelType       InputImagePixelType;
  typedef typename OutputImageType::Pointer        OutputImagePointer;
  typedef typename OutputImageType::ConstPointer   OutputImageConstPointer;
  typedef typename OutputImageType::RegionType     OutputImageRegionType;
  typedef typename OutputImageType::PixelType      OutputImagePixelType;

  /** ImageDimension constants */
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);
  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      TOutputImage::ImageDimension);

  /** Standard New method. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(AnchorErodeDilateImageFilter,
               ImageToImageFilter);

  itkSetMacro(Size, unsigned int);
  itkGetConstReferenceMacro(Size, unsigned int);
  itkBooleanMacro(Size);

  itkSetMacro(Direction, unsigned int);
  itkGetConstReferenceMacro(Direction, unsigned int);
  itkBooleanMacro(Direction);


protected:
  AnchorErodeDilateImageFilter();
  ~AnchorErodeDilateImageFilter() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Single-threaded version of GenerateData.  This filter delegates
   * to GrayscaleGeodesicErodeImageFilter. */
  void GenerateData();


private:
  AnchorErodeDilateImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  unsigned int m_Size;
  unsigned int m_Direction;
  TFunction1 m_TF1;
  TFunction2 m_TF2;

  typedef MorphologyHistogram<InputImagePixelType> Histogram;
  typedef MorphologyHistogramVec<InputImagePixelType,THistogramCompare> VHistogram;
  typedef MorphologyHistogramMap<InputImagePixelType,THistogramCompare> MHistogram;

  bool startLine(OutputImagePixelType * buffer,
		 OutputImagePixelType * inbuffer,
		 OutputImagePixelType &Extreme,
#ifdef RAWHIST
		 unsigned int *  histo,
#else
		 Histogram &histo,
#endif
		 unsigned &outLeftP,
		 unsigned &outRightP,
		 unsigned &inLeftP,
		 unsigned &inRightP,
		 unsigned middle);

  bool finishLine(OutputImagePixelType * buffer,
		  OutputImagePixelType * inbuffer,
		  OutputImagePixelType &Extreme,
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


  void check(unsigned val, unsigned start, unsigned end, std::string text=std::string(""))
  {
    if ((val < start) || (val >= end))
      {
      std::cout << "******************" << text << " " << val << std::endl;
      }
  }
} ; // end of class


} // end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkAnchorErodeDilateImageFilter.txx"
#endif

#endif

