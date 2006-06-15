#ifndef __itkAnchorErodeDilateImageFilter_txx
#define __itkAnchorErodeDilateImageFilter_txx

#include "itkAnchorErodeDilateImageFilter.h"
#include "itkImageLinearConstIteratorWithIndex.h"
#include "itkImageLinearIteratorWithIndex.h"

namespace itk {

template <class TInputImage, class TOutputImage, class THistogram, class TFunction1, class TFunction2>
AnchorErodeDilateImageFilter<TInputImage, TOutputImage, THistogram, TFunction1, TFunction2>
::AnchorErodeDilateImageFilter()
{
  m_Direction = 0;
  m_Size=2;
}

template <class TInputImage, class TOutputImage, class THistogram, class TFunction1, class TFunction2>
void
AnchorErodeDilateImageFilter<TInputImage, TOutputImage, THistogram, TFunction1, TFunction2>
::GenerateData()
{
  // TFunction1 will be < for erosions
  // TFunction2 will be <=


  // the initial version will adopt the methodology of loading a line
  // at a time into a buffer vector, carrying out the opening or
  // closing, and then copy the result to the output. Hopefully this
  // will improve cache performance when working along non raster
  // directions.

  // Allocate the output
  this->AllocateOutputs();
  typename TOutputImage::Pointer output = this->GetOutput();
  typename TInputImage::ConstPointer input = this->GetInput();
  
  // get the region size
  OutputImageRegionType OReg = output->GetRequestedRegion();
  unsigned int bufflength = OReg.GetSize()[m_Direction];
  unsigned linecount = OReg.GetNumberOfPixels()/bufflength;
  ProgressReporter progress(this, 0, linecount);

  OutputImagePixelType * buffer = new OutputImagePixelType[bufflength];
  OutputImagePixelType * inbuffer = new OutputImagePixelType[bufflength];

  // set up a line iterator
  typedef itk::ImageLinearConstIteratorWithIndex<InputImageType> InputLineIteratorType;
  typedef itk::ImageLinearIteratorWithIndex<OutputImageType> OutputLineIteratorType;

  InputLineIteratorType inLineIt(input, output->GetRequestedRegion());
  inLineIt.SetDirection(m_Direction);
  OutputLineIteratorType outLineIt(output, output->GetRequestedRegion());
  outLineIt.SetDirection(m_Direction);

  unsigned int middle = m_Size/2;

#ifdef RAWHIST
  unsigned int * histo = new unsigned int[256];
#else
  // create a histogram
  Histogram *histo;
  if (useVectorBasedHistogram())
    {
    histo = new VHistogram;
    } 
  else
    {
    histo = new MHistogram;
    }
#endif

  for (inLineIt.GoToBegin(), outLineIt.GoToBegin(); ! inLineIt.IsAtEnd(); inLineIt.NextLine(),outLineIt.NextLine())
    {
    // copy the line to the buffer
    inLineIt.GoToBeginOfLine();
    unsigned pos = 0;
    while (! inLineIt.IsAtEndOfLine())
      {
      OutputImagePixelType PVal = static_cast<OutputImagePixelType>(inLineIt.Get());
//      buffer[pos]=PVal;
      inbuffer[pos]=PVal;
      ++pos; 
      ++inLineIt;
      }
    // done copying
    // start the real work - everything here will be done with index
    // arithmetic rather than pointer arithmetic
    unsigned outLeftP = 0, outRightP = bufflength - 1;
    unsigned inLeftP = 0, inRightP = bufflength - 1;
    OutputImagePixelType Extreme;
    histo->Reset();

    // Left border, first half of structuring element
    Extreme = inbuffer[inLeftP];
    histo->AddPixel(Extreme);
    for (unsigned i = 0; i < middle; i++)
      {
      ++inLeftP;
      histo->AddPixel(inbuffer[inLeftP]);
      if (m_TF1(inbuffer[inLeftP], Extreme))
	{
	Extreme = inbuffer[inLeftP];
	}
      }
    buffer[outLeftP] = Extreme;
    
    // Second half of SE
    for (unsigned i = 0; i < m_Size - middle - 1; i++)
      {
      ++inLeftP;
      ++outLeftP;
      histo->AddPixel(inbuffer[inLeftP]);
      if (m_TF1(inbuffer[inLeftP], Extreme))
	{
	Extreme = inbuffer[inLeftP];
	}
      buffer[outLeftP] = Extreme;
      }
    // Use the histogram until we find a new minimum 
    while ((inLeftP < inRightP) && m_TF2(Extreme, inbuffer[inLeftP + 1]))
      {
      ++inLeftP;
      ++outLeftP;
      histo->RemovePixel(inbuffer[inLeftP - m_Size]);
      histo->AddPixel(inbuffer[inLeftP]);
      Extreme = histo->GetValue();
      buffer[outLeftP] = Extreme;
      }
    Extreme = buffer[outLeftP];

    while (startLine(buffer, inbuffer, Extreme, *histo, outLeftP, outRightP, inLeftP, inRightP, middle)){}

    finishLine(buffer, inbuffer, Extreme, *histo, outLeftP, outRightP, inLeftP, inRightP, middle);

    // copy the buffer to output
    inLineIt.GoToBeginOfLine();
    pos = 0;
    while (! outLineIt.IsAtEndOfLine())
      {
      outLineIt.Set(buffer[pos]);
      ++pos;
      ++outLineIt;
      }
    // done with this line
    progress.CompletedPixel();
    }

  delete [] buffer;
  delete [] inbuffer;

}

template<class TInputImage, class TOutputImage, class THistogram, class TFunction1, class TFunction2>
bool
AnchorErodeDilateImageFilter<TInputImage, TOutputImage, THistogram, TFunction1, TFunction2>
::startLine(OutputImagePixelType * buffer,
	    OutputImagePixelType * inbuffer,
	    OutputImagePixelType &Extreme,
#ifdef RAWHIST
	    unsigned int * histo,	    
#else
	    Histogram &histo,
#endif
	    unsigned &outLeftP,
	    unsigned &outRightP,
	    unsigned &inLeftP,
	    unsigned &inRightP,
	    unsigned middle)
{
  // This returns true to indicate return to startLine label in pseudo
  // code, and false to indicate finshLine
  unsigned currentP = inLeftP + 1;
  unsigned sentinel;
  
  while ((currentP < inRightP) && m_TF2(inbuffer[currentP], Extreme))
    {
    Extreme = inbuffer[currentP];
    ++outLeftP;
    buffer[outLeftP] = Extreme;
    ++currentP;
    }
  inLeftP = currentP - 1;

  sentinel = inLeftP + m_Size;
  if (sentinel > inRightP)
    {
    // finish
    return (false);
    }
  ++outLeftP;
  buffer[outLeftP] = Extreme;

  // ran m_Size pixels ahead
  ++currentP;
  while (currentP < sentinel)
    {
    if (m_TF2(inbuffer[currentP], Extreme))
      {
      Extreme = inbuffer[currentP];
      ++outLeftP;
      buffer[outLeftP] = Extreme;
      inLeftP = currentP;
      return (true);
      }
    ++currentP;
    ++outLeftP;
    buffer[outLeftP] = Extreme;
    }
  // We didn't find a smaller (for erosion) value in the segment of
  // reach of inLeftP. currentP is the first position outside the
  // reach of inLeftP
  if (m_TF2(inbuffer[currentP], Extreme))
    {
    Extreme = inbuffer[currentP];
    ++outLeftP;
    buffer[outLeftP] = Extreme;
    inLeftP = currentP;
    return (true);
    }
  else
    {
    // Now we need a histogram
    // Initialise it
    histo.Reset();
    ++outLeftP;
    ++inLeftP;
    for (unsigned aux = inLeftP; aux <= currentP; ++aux)
      {
      histo.AddPixel(inbuffer[aux]);
      }
    Extreme = histo.GetValue();
    buffer[outLeftP] = Extreme;
    }

  while (currentP < inRightP)
    {
    ++currentP;
    if (m_TF2(inbuffer[currentP], Extreme))
      {
      // Found a new extrem
      Extreme = inbuffer[currentP];
      ++outLeftP;
      buffer[outLeftP] = Extreme;
      inLeftP = currentP;
      return(true);
      }
    else
      {
      // update histogram
      histo.AddPixel(inbuffer[currentP]);
      histo.RemovePixel(inbuffer[inLeftP]);
      // find extreme
      Extreme = histo.GetValue();
      ++inLeftP;
      ++outLeftP;
      buffer[outLeftP] = Extreme;
      }
    }
  return(false);
}

template<class TInputImage, class TOutputImage, class THistogram, class TFunction1, class TFunction2>
bool
AnchorErodeDilateImageFilter<TInputImage, TOutputImage, THistogram, TFunction1, TFunction2>
::finishLine(OutputImagePixelType * buffer,
	     OutputImagePixelType * inbuffer,
	     OutputImagePixelType &Extreme,
	     Histogram &histo,
	     unsigned &outLeftP,
	     unsigned &outRightP,
	     unsigned &inLeftP,
	     unsigned &inRightP,
	     unsigned middle)
{
  // Handles the right border.
  // First half of the structuring element
  histo.Reset();
  Extreme = inbuffer[inRightP];
  histo.AddPixel(Extreme);

  for (unsigned i = 0; i < middle; i++)
    {
    --inRightP; 
    histo.AddPixel(inbuffer[inRightP]);
    if (m_TF1(inbuffer[inRightP], Extreme))
      {
      Extreme = inbuffer[inRightP];
      }
    }
  buffer[outRightP] = Extreme;
  // second half of SE
  for (unsigned i = 0; (i<m_Size - middle - 1) && (outLeftP < outRightP); i++)
    {
    --inRightP;
    --outRightP;
    histo.AddPixel(inbuffer[inRightP]);
    if (m_TF1(inbuffer[inRightP], Extreme))
      {
      Extreme = inbuffer[inRightP];
      }
    buffer[outRightP] = Extreme;
    }

  while (outLeftP < outRightP)
    {
    --inRightP;
    --outRightP;
    histo.RemovePixel(inbuffer[inRightP + m_Size]);
    histo.AddPixel(inbuffer[inRightP]);
    if (m_TF1(inbuffer[inRightP], Extreme))
      {
      Extreme = inbuffer[inRightP];
      }
    Extreme = histo.GetValue();
    buffer[outRightP] = Extreme;
    }
  
  
}

template<class TInputImage, class TOutputImage, class THistogram, class TFunction1, class TFunction2>
void
AnchorErodeDilateImageFilter<TInputImage, TOutputImage, THistogram, TFunction1, TFunction2>
::PrintSelf(std::ostream &os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "Direction: " << m_Direction << std::endl;
  os << indent << "Size: " << m_Size << std::endl;
}


} // end namespace itk

#endif
