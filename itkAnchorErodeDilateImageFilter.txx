
#ifndef __itkAnchorErodeDilateImageFilter_txx
#define __itkAnchorErodeDilateImageFilter_txx

#include "itkAnchorErodeDilateImageFilter.h"

//#include "itkNeighborhoodAlgorithm.h"

#include "itkAnchorUtilities.h"

namespace itk {

template <class TImage, class TKernel, class TFunction1, class TFunction2>
AnchorErodeDilateImageFilter<TImage, TKernel, TFunction1, TFunction2>
::AnchorErodeDilateImageFilter()
{
  m_KernelSet = false;
}

template <class TImage, class TKernel, class TFunction1, class TFunction2>
void
AnchorErodeDilateImageFilter<TImage, TKernel, TFunction1, TFunction2>
::GenerateData()
{

  // check that we are using a decomposable kernel
  if (!m_Kernel.GetDecomposable())
    {
    itkWarningMacro("Anchor morphology only works with decomposable structuring elements");
    return;
    }
  if (!m_KernelSet)
    {
    itkWarningMacro("No kernel set - quitting");
    return;
    }
  // TFunction1 will be < for erosions
  // TFunction2 will be <=


  // the initial version will adopt the methodology of loading a line
  // at a time into a buffer vector, carrying out the opening or
  // closing, and then copy the result to the output. Hopefully this
  // will improve cache performance when working along non raster
  // directions.

  // Allocate the output
  this->AllocateOutputs();
  InputImagePointer output = this->GetOutput();
  InputImageConstPointer input = this->GetInput();

  // get the region size
  InputImageRegionType OReg = output->GetRequestedRegion();
  // maximum buffer length is sum of dimensions
  unsigned int bufflength = 0;
  for (unsigned i = 0; i<TImage::ImageDimension; i++)
    {
    bufflength += OReg.GetSize()[i];
    }
  

  InputImagePixelType * buffer = new InputImagePixelType[bufflength];
  InputImagePixelType * inbuffer = new InputImagePixelType[bufflength];

  // iterate over all the structuring elements
  typename KernelType::DecompType decomposition = m_Kernel.GetLines();
  BresType BresLine;
  ProgressReporter progress(this, 0, decomposition.size());

  std::cout << decomposition.size() << " lines will be used" << std::endl;

  for (unsigned i = 0; i < decomposition.size(); i++)
    {
    typename KernelType::LType ThisLine = decomposition[i];
    typename BresType::OffsetArray TheseOffsets = BresLine.buildLine(ThisLine, bufflength);
    unsigned int SELength = getLinePixels<typename KernelType::LType>(ThisLine);
    // want lines to be odd
    if (!(SELength%2))
      ++SELength;

    std::cout << ThisLine << SELength << std::endl;
    AnchorLine.SetSize(SELength);

    InputImageRegionType BigFace = mkEnlargedFace<InputImageType, typename KernelType::LType>(input, OReg, ThisLine);
    doFace<TImage, BresType, AnchorLineType, typename KernelType::LType>(input, output, ThisLine, AnchorLine, 
									   TheseOffsets, inbuffer, buffer, OReg, BigFace);

    // after the first pass the input will be taken from the output
    input = this->GetOutput();
    progress.CompletedPixel();
    }


  delete [] buffer;
  delete [] inbuffer;
}


template<class TImage, class TKernel, class TFunction1, class TFunction2>
void
AnchorErodeDilateImageFilter<TImage, TKernel, TFunction1, TFunction2>
::PrintSelf(std::ostream &os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}


} // end namespace itk

#endif
