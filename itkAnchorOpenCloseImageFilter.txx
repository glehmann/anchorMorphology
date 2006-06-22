#ifndef __itkAnchorOpenCloseImageFilter_txx
#define __itkAnchorOpenCloseImageFilter_txx

#include "itkAnchorOpenCloseImageFilter.h"
#include "itkImageLinearConstIteratorWithIndex.h"
#include "itkImageLinearIteratorWithIndex.h"

namespace itk {

template <class TInputImage, class TOutputImage, class THistogramCompare, class TFunction1, class TFunction2>
AnchorOpenCloseImageFilter<TInputImage, TOutputImage, THistogramCompare, TFunction1, TFunction2>
::AnchorOpenCloseImageFilter()
{
  m_Direction = 0;
  m_Size=2;
  AnchorLine.SetSize(m_Size);
}

template <class TInputImage, class TOutputImage, class THistogramCompare, class TFunction1, class TFunction2>
void
AnchorOpenCloseImageFilter<TInputImage, TOutputImage, THistogramCompare, TFunction1, TFunction2>
::GenerateData()
{
  // TFunction1 will be >= for openings
  // TFunction2 will be <=
  // TFunction3 will be >


  // the initial version will adopt the methodology of loading a line
  // at a time into a buffer vector, carrying out the opening or
  // closing, and then copy the result to the output. Hopefully this
  // will improve cache performance when working along non raster
  // directions.
  AnchorLine.SetSize(m_Size);
  
  // Allocate the output
  this->AllocateOutputs();
  typename TOutputImage::Pointer output = this->GetOutput();
  typename TInputImage::ConstPointer input = this->GetInput();
  
  // get the region size
  OutputImageRegionType OReg = output->GetRequestedRegion();
  unsigned int bufflength = OReg.GetSize()[m_Direction];
  unsigned linecount = OReg.GetNumberOfPixels()/bufflength;
  ProgressReporter progress(this, 0, linecount);

  InputImagePixelType * buffer = new OutputImagePixelType[bufflength];

  // set up a line iterator
  typedef itk::ImageLinearConstIteratorWithIndex<InputImageType> InputLineIteratorType;
  typedef itk::ImageLinearIteratorWithIndex<OutputImageType> OutputLineIteratorType;

  InputLineIteratorType inLineIt(input, output->GetRequestedRegion());
  inLineIt.SetDirection(m_Direction);
  OutputLineIteratorType outLineIt(output, output->GetRequestedRegion());
  outLineIt.SetDirection(m_Direction);

  for (inLineIt.GoToBegin(), outLineIt.GoToBegin(); ! inLineIt.IsAtEnd(); inLineIt.NextLine(),outLineIt.NextLine())
    {
    // copy the line to the buffer
    inLineIt.GoToBeginOfLine();
    unsigned pos = 0;
    while (! inLineIt.IsAtEndOfLine())
      {
      InputImagePixelType PVal = inLineIt.Get();
      buffer[pos]=PVal;
      ++pos; 
      ++inLineIt;
      }
    // done copying
    // start the real work - everything here will be done with index
    AnchorLine.doLine(buffer, bufflength);
    // copy the buffer to output
    inLineIt.GoToBeginOfLine();
    pos = 0;
    while (! outLineIt.IsAtEndOfLine())
      {
      outLineIt.Set(static_cast<OutputImagePixelType>(buffer[pos]));
      ++pos;
      ++outLineIt;
      }
    // done with this line
    progress.CompletedPixel();

    }

  delete [] buffer;
}


template<class TInputImage, class TOutputImage, class THistogramCompare, class TFunction1, class TFunction2>
void
AnchorOpenCloseImageFilter<TInputImage, TOutputImage, THistogramCompare, TFunction1, TFunction2>
::PrintSelf(std::ostream &os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "Direction: " << m_Direction << std::endl;
  os << indent << "Size: " << m_Size << std::endl;
}


} // end namespace itk

#endif
