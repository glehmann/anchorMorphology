#ifndef __itkAnchorOpenCloseImageFilter_h
#define __itkAnchorOpenCloseImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkProgressReporter.h"

//#define RAWHIST

namespace itk {



/** 
 * \class AnchorOpenCloseImageFilter
 * \brief class to implement openings and closings using anchor
 * methods. This is the base class that must be instantiated with
 * appropriate definitions of greater, less and so on

**/
template<class TInputImage, class TOutputImage, class THistogram,
	 class TFunction1, class TFunction2>
class ITK_EXPORT AnchorOpenCloseImageFilter :
    public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef AnchorOpenCloseImageFilter Self;
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
  itkTypeMacro(AnchorOpenCloseImageFilter,
               ImageToImageFilter);

  itkSetMacro(Size, unsigned int);
  itkGetConstReferenceMacro(Size, unsigned int);
  itkBooleanMacro(Size);

  itkSetMacro(Direction, unsigned int);
  itkGetConstReferenceMacro(Direction, unsigned int);
  itkBooleanMacro(Direction);


protected:
  AnchorOpenCloseImageFilter();
  ~AnchorOpenCloseImageFilter() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Single-threaded version of GenerateData.  This filter delegates
   * to GrayscaleGeodesicErodeImageFilter. */
  void GenerateData();


private:
  AnchorOpenCloseImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  unsigned int m_Size;
  unsigned int m_Direction;
  TFunction1 m_TF1;
  TFunction2 m_TF2;

  bool startLine(OutputImagePixelType * buffer,
		 OutputImagePixelType &Extreme,
#ifdef RAWHIST
		 unsigned int *  histo,
#else
		 THistogram &histo,
#endif
		 unsigned &outLeftP,
		 unsigned &outRightP);

  bool finishLine(OutputImagePixelType * buffer,
		  OutputImagePixelType &Extreme,
		  unsigned &outLeftP,
		  unsigned &outRightP);



} ; // end of class


} // end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkAnchorOpenCloseImageFilter.txx"
#endif

#endif

