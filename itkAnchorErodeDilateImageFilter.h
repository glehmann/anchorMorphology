#ifndef __itkAnchorErodeDilateImageFilter_h
#define __itkAnchorErodeDilateImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkProgressReporter.h"
#include "itkAnchorErodeDilateLine.h"
#include "itkBresenhamLine.h"

#define ANCHOR_ALGORITHM

namespace itk {

/** 
 * \class AnchorErodeDilateImageFilter
 * \brief class to implement erosions and dilations using anchor
 * methods. This is the base class that must be instantiated with
 * appropriate definitions of greater, less and so on

**/
template<class TImage, class TKernel, 
	 class TFunction1, class TFunction2>
class ITK_EXPORT AnchorErodeDilateImageFilter :
    public ImageToImageFilter<TImage, TImage>
{
public:
  /** Standard class typedefs. */
  typedef AnchorErodeDilateImageFilter Self;
  typedef ImageToImageFilter<TImage, TImage>
  Superclass;
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Some convenient typedefs. */
  /** Kernel typedef. */
  typedef TKernel KernelType;

  typedef TImage InputImageType;
  typedef typename InputImageType::Pointer         InputImagePointer;
  typedef typename InputImageType::ConstPointer    InputImageConstPointer;
  typedef typename InputImageType::RegionType      InputImageRegionType;
  typedef typename InputImageType::PixelType       InputImagePixelType;
  typedef typename TImage::IndexType         IndexType;
  typedef typename TImage::SizeType          SizeType;

  /** ImageDimension constants */
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TImage::ImageDimension);
  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      TImage::ImageDimension);

  /** Standard New method. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(AnchorErodeDilateImageFilter,
               ImageToImageFilter);

  void SetKernel( const KernelType& kernel )
  {
    m_Kernel=kernel;
    m_KernelSet = true;
  }



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

  TKernel m_Kernel;
  bool m_KernelSet;
  typedef BresenhamLine<TImage::ImageDimension> BresType;

#ifdef ANCHOR_ALGORITHM
  // the class that operates on lines
  typedef AnchorErodeDilateLine<InputImagePixelType, TFunction1, TFunction2> AnchorLineType;

  AnchorLineType AnchorLine;
#endif

} ; // end of class


} // end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkAnchorErodeDilateImageFilter.txx"
#endif

#endif

