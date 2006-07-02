#ifndef __itkAnchorOpenCloseImageFilter_txx
#define __itkAnchorOpenCloseImageFilter_txx

#include "itkAnchorOpenCloseImageFilter.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkAnchorUtilities.h"

namespace itk {

template <class TImage, class TKernel, class LessThan, class GreaterThan, class LessEqual, class GreaterEqual>
AnchorOpenCloseImageFilter<TImage, TKernel, LessThan, GreaterThan, LessEqual, GreaterEqual>
::AnchorOpenCloseImageFilter()
{
  m_KernelSet = false;
}

template <class TImage, class TKernel, class LessThan, class GreaterThan, class LessEqual, class GreaterEqual>
void
AnchorOpenCloseImageFilter<TImage, TKernel, LessThan, GreaterThan, LessEqual, GreaterEqual>
::GenerateData()
{

  // check that we are using a decomposable kernel
  if (!m_Kernel.Decomposable())
    {
    itkWarningMacro("Anchor morphology only works with decomposable structuring elements");
    return;
    }
  if (!m_KernelSet)
    {
    itkWarningMacro("No kernel set - quitting");
    return;
    }

  
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

  //unsigned linecount = OReg.GetNumberOfPixels()/bufflength;

  InputImagePixelType * inbuffer = new InputImagePixelType[bufflength];
  InputImagePixelType * outbuffer = new InputImagePixelType[bufflength];

  // iterate over all the structuring elements
  typename KernelType::DecompType decomposition = m_Kernel.GetDecomp();
  BresType BresLine;
  ProgressReporter progress(this, 0, decomposition.size()*2);

  // first stage -- all of the erosions if we are doing an opening
  for (unsigned i = 0; i < decomposition.size() - 1; i++)
    {
    typename KernelType::LType ThisLine = decomposition[i];
    typename BresType::OffsetArray TheseOffsets = BresLine.buildLine(ThisLine, bufflength);
    unsigned int SELength = getLinePixels<typename KernelType::LType>(ThisLine);
    // want lines to be odd
    if (!(SELength%2))
      ++SELength;
    AnchorLineErode.SetSize(SELength);
    typedef typename std::list<InputImageRegionType> FaceListType;
    typedef typename FaceListType::iterator FaceListIterator;
    FaceListType faceList;
    FaceListIterator fit;
    // Now figure out which faces of the image we should be starting
    // from with this line
    faceList = mkFaceList<InputImageRegionType, typename KernelType::LType>(OReg, ThisLine);

    for ( fit = faceList.begin(); fit != faceList.end(); ++fit)
      {
      //std::cout << ThisLine << std::endl;
      //std::cout << (*fit) << std::endl;
      doFace<TImage, BresType, AnchorLineErodeType>(input, output, AnchorLineErode,
						    TheseOffsets, inbuffer, outbuffer, 
						    OReg, *fit);
      }
    //std::cout << "-------------------" << std::endl;
    // after the first pass the input will be taken from the output
    input = this->GetOutput();
    progress.CompletedPixel();
    }
  // now do the opening in the middle of the chain
  {
  unsigned i = decomposition.size() - 1;
  typename KernelType::LType ThisLine = decomposition[i];
  typename BresType::OffsetArray TheseOffsets = BresLine.buildLine(ThisLine, bufflength);
  unsigned int SELength = getLinePixels<typename KernelType::LType>(ThisLine);
  // want lines to be odd
  if (!(SELength%2))
    ++SELength;

  AnchorLineOpen.SetSize(SELength);

  typedef typename std::list<InputImageRegionType> FaceListType;
  typedef typename FaceListType::iterator FaceListIterator;
  FaceListType faceList;
  FaceListIterator fit;
  // Now figure out which faces of the image we should be starting
  // from with this line
  faceList = mkFaceList<InputImageRegionType, typename KernelType::LType>(OReg, ThisLine);

  for ( fit = faceList.begin(); fit != faceList.end(); ++fit)
    {
    doFaceOpen(input, output,
	       TheseOffsets, outbuffer, 
	       OReg, *fit);
    }
  // equivalent to two passes
  progress.CompletedPixel();
  progress.CompletedPixel();  
  }

  // Now for the rest of the dilations -- note that i needs to be signed
  for (int i = decomposition.size() - 2; i >= 0; --i)
    {
    typename KernelType::LType ThisLine = decomposition[i];
    typename BresType::OffsetArray TheseOffsets = BresLine.buildLine(ThisLine, bufflength);
    unsigned int SELength = getLinePixels<typename KernelType::LType>(ThisLine);
    // want lines to be odd
    if (!(SELength%2))
      ++SELength;
  
    AnchorLineDilate.SetSize(SELength);
    typedef typename std::list<InputImageRegionType> FaceListType;
    typedef typename FaceListType::iterator FaceListIterator;
    FaceListType faceList;
    FaceListIterator fit;
    // Now figure out which faces of the image we should be starting
    // from with this line
    faceList = mkFaceList<InputImageRegionType, typename KernelType::LType>(OReg, ThisLine);

    for ( fit = faceList.begin(); fit != faceList.end(); ++fit)
      {

      doFace<TImage, BresType, AnchorLineDilateType>(input, output, AnchorLineDilate,
						     TheseOffsets, inbuffer, outbuffer, 
						     OReg, *fit);

      }
    progress.CompletedPixel();
    }

  delete [] inbuffer;
  delete [] outbuffer;
}

template<class TImage, class TKernel, class LessThan, class GreaterThan, class LessEqual, class GreaterEqual>
void
AnchorOpenCloseImageFilter<TImage, TKernel, LessThan, GreaterThan, LessEqual, GreaterEqual>
::doFaceOpen(InputImageConstPointer input,
	     InputImagePointer output,
	     typename BresType::OffsetArray LineOffsets,
	     InputImagePixelType * outbuffer,	      
	     const InputImageRegionType AllImage, 
	     const InputImageRegionType face)
{
  // iterate over the face
  typedef ImageRegionConstIteratorWithIndex<InputImageType> ItType;
  ItType it(input, face);
  it.GoToBegin();
  while (!it.IsAtEnd()) 
    {
    typename TImage::IndexType Ind = it.GetIndex();
    unsigned int len;
    fillLineBuffer<TImage, BresType>(input, Ind, LineOffsets, AllImage, outbuffer, len);
    AnchorLineOpen.doLine(outbuffer,len);
    copyLineToImage<TImage, BresType>(output, Ind, LineOffsets, outbuffer, len);
    ++it;
    }
}

template<class TImage, class TKernel, class LessThan, class GreaterThan, class LessEqual, class GreaterEqual>
void
AnchorOpenCloseImageFilter<TImage, TKernel, LessThan, GreaterThan, LessEqual, GreaterEqual>
::PrintSelf(std::ostream &os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}


} // end namespace itk

#endif
