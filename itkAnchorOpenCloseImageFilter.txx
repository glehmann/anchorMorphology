#ifndef __itkAnchorOpenCloseImageFilter_txx
#define __itkAnchorOpenCloseImageFilter_txx

#include "itkAnchorOpenCloseImageFilter.h"
#include "itkImageLinearConstIteratorWithIndex.h"
#include "itkImageLinearIteratorWithIndex.h"

namespace itk {

template <class TInputImage, class TOutputImage, class THistogram, class TFunction1, class TFunction2>
AnchorOpenCloseImageFilter<TInputImage, TOutputImage, THistogram, TFunction1, TFunction2>
::AnchorOpenCloseImageFilter()
{
  m_Direction = 0;
  m_Size=2;
}

template <class TInputImage, class TOutputImage, class THistogram, class TFunction1, class TFunction2>
void
AnchorOpenCloseImageFilter<TInputImage, TOutputImage, THistogram, TFunction1, TFunction2>
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

  // set up a line iterator
  typedef itk::ImageLinearConstIteratorWithIndex<InputImageType> InputLineIteratorType;
  typedef itk::ImageLinearIteratorWithIndex<OutputImageType> OutputLineIteratorType;

  InputLineIteratorType inLineIt(input, output->GetRequestedRegion());
  inLineIt.SetDirection(m_Direction);
  OutputLineIteratorType outLineIt(output, output->GetRequestedRegion());
  outLineIt.SetDirection(m_Direction);

#ifdef RAWHIST
  unsigned int * histo = new unsigned int[256];
#else
  // create a histogram
  THistogram histo;
#endif

  for (inLineIt.GoToBegin(), outLineIt.GoToBegin(); ! inLineIt.IsAtEnd(); inLineIt.NextLine(),outLineIt.NextLine())
    {
    // copy the line to the buffer
    inLineIt.GoToBeginOfLine();
    unsigned pos = 0;
    while (! inLineIt.IsAtEndOfLine())
      {
      OutputImagePixelType PVal = static_cast<OutputImagePixelType>(inLineIt.Get());
      buffer[pos]=PVal;
      ++pos; 
      ++inLineIt;
      }
    // done copying
    // start the real work - everything here will be done with index
    // arithmetic rather than pointer arithmetic
    unsigned outLeftP = 0, outRightP = bufflength - 1;
    // left side
    while ((outLeftP < outRightP) && m_TF1(buffer[outLeftP], buffer[outLeftP+1]))
      {
      ++outLeftP;
      }
    while ((outLeftP < outRightP) && m_TF2(buffer[outRightP-1], buffer[outRightP]))
      {
      --outRightP;
      }
    OutputImagePixelType Extreme;
    while (startLine(buffer, Extreme, histo, outLeftP, outRightP)){}

    finishLine(buffer, Extreme, outLeftP, outRightP);
    std::cout << "------------------------------" << std::endl;

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
//  delete [] histo;
}

template<class TInputImage, class TOutputImage, class THistogram, class TFunction1, class TFunction2>
bool
AnchorOpenCloseImageFilter<TInputImage, TOutputImage, THistogram, TFunction1, TFunction2>
::startLine(OutputImagePixelType * buffer,
	    OutputImagePixelType &Extreme,
#ifdef RAWHIST
	    unsigned int * histo,	    
#else
	    THistogram &histo,
#endif
	    unsigned &outLeftP,
	    unsigned &outRightP)
{
  // This returns true to indicate return to startLine label in pseudo
  // code, and false to indicate finshLine
  Extreme = buffer[outLeftP];
  unsigned currentP = outLeftP + 1;
  unsigned sentinel, endP;
  
  while ((currentP < outRightP) && m_TF2(buffer[currentP], Extreme))
    {
    Extreme = buffer[currentP];
    ++outLeftP;
    ++currentP;
    }

  sentinel = outLeftP + m_Size;
  if (sentinel > outRightP)
    {
    // finish
    return (false);
    }
  ++currentP;
  // ran m_Size pixels ahead
  while (currentP < sentinel)
    {
    if (m_TF2(buffer[currentP], Extreme))
      {
#if 0
      endP = currentP;
      ++outLeftP;
      while (outLeftP < endP)
	{
	buffer[outLeftP] = Extreme;
	++outLeftP;
	}
      outLeftP = currentP;
#else
      endP = currentP;
      for (unsigned PP = outLeftP + 1; PP < endP; ++PP)
	{
	buffer[PP] = Extreme;
	}
      outLeftP = currentP;
#endif
      return (true);
      }
    ++currentP;
    }
  // We didn't find a smaller (for opening) value in the segment of
  // reach of outLeftP. currentP is the first position outside the
  // reach of outLeftP
  if (m_TF2(buffer[currentP], Extreme))
    {
#if 0
    endP = currentP;
    ++outLeftP;
    while(outLeftP < endP)
      {
      buffer[outLeftP] = Extreme;
      ++outLeftP;
      }
    outLeftP = currentP;
#else
    endP = currentP;
    for (unsigned PP = outLeftP + 1; PP < endP; ++PP)
      {
      buffer[PP] = Extreme;
      }
    outLeftP = currentP;
#endif
    return(true);
    }
  else
    {
    // Now we need a histogram
    // Initialise it
#ifdef RAWHIST
    std::fill(histo, &(histo[256]), 0);
#else
    histo.Init();
#endif
    ++outLeftP;
    for (unsigned aux = outLeftP; aux <= currentP; ++aux)
      {
      std::cout << "Adding " << (int)buffer[aux] << std::endl;
#ifdef RAWHIST
      histo[buffer[aux]]++;
#else
      histo.AddPixel(buffer[aux]);
#endif
      }
    // find the minimum value. The version
    // in the paper assumes integer pixel types and initializes the
    // search to the current extreme. Hopefully the latter is an
    // optimization step.
#ifdef RAWHIST
    ++Extreme;
    while (histo[Extreme] <= 0) {++Extreme;}
    std::cout << "A) Extreme from hist = " << (int) Extreme << std::endl;
    histo[buffer[outLeftP]]--;
    std::cout << "Removing " << (int)buffer[outLeftP] << std::endl;
    buffer[outLeftP] = Extreme;
    histo[Extreme]++;
    std::cout << "Adding " << (int)Extreme << std::endl;
#else
    Extreme = histo.GetValue();
    std::cout << "A) Extreme from hist = " << (int) Extreme << std::endl;
    histo.RemovePixel(buffer[outLeftP]);
    std::cout << "Removing " << (int)buffer[outLeftP] << std::endl;
    buffer[outLeftP] = Extreme;
    histo.AddPixel(Extreme);
    std::cout << "Adding " << (int)Extreme << std::endl;
#endif
    }

  while (currentP < outRightP)
    {
    ++currentP;
    if (m_TF2(buffer[currentP], Extreme))
      {
      // Found a new extrem
      endP = currentP;
#if 0
      ++outLeftP;
      while (outLeftP < endP)
	{
	buffer[outLeftP] = Extreme;
	++outLeftP;
	}
#else
      for (unsigned PP = outLeftP + 1; PP < endP; PP++)
	{
	buffer[PP]=Extreme;
	}
#endif
      outLeftP = currentP;
      return(true);
      }
    else
      {
      /* histogram update */
#ifdef RAWHIST
      histo[buffer[currentP]]++;
      std::cout << "Adding " << (int)buffer[currentP] << std::endl;
      histo[buffer[outLeftP]]--;
      std::cout << "Removing " << (int)buffer[outLeftP] << std::endl;
      while (histo[Extreme] <= 0) {Extreme++;}
      ++outLeftP;
      histo[buffer[outLeftP]]--;
      std::cout << "Removing " << (int)buffer[outLeftP] << std::endl;
      buffer[outLeftP] = Extreme;
      histo[Extreme]++;
      std::cout << "B) Extreme from hist = " << (int) Extreme << std::endl;
#else
      histo.AddPixel(buffer[currentP]);
      histo.RemovePixel(buffer[outLeftP]);
      Extreme = histo.GetValue();
      ++outLeftP;
      histo.RemovePixel(buffer[outLeftP]);
      buffer[outLeftP] = Extreme;
      histo.AddPixel(Extreme);
      std::cout << "B) Extreme from hist = " << (int) Extreme << std::endl;
#endif
      }
    }
  // Finish the line
  while (outLeftP < outRightP)
    {
#ifdef RAWHIST
    histo[buffer[outLeftP]]--;
    while (histo[Extreme] <= 0) { Extreme++;}
 
    ++outLeftP;
    histo[buffer[outLeftP]]--;
    
    buffer[outLeftP] = Extreme;
    histo[Extreme]++;
    std::cout << "C) Extreme from hist = " << (int) Extreme << std::endl;

#else
    histo.RemovePixel(buffer[outLeftP]);
    Extreme = histo.GetValue();
    ++outLeftP;
    histo.RemovePixel(buffer[outLeftP]);
    buffer[outLeftP] = Extreme;
    histo.AddPixel(Extreme);
    std::cout << "C) Extreme from hist = " << (int) Extreme << std::endl;
#endif
    }
  return(false);
}

template<class TInputImage, class TOutputImage, class THistogram, class TFunction1, class TFunction2>
bool
AnchorOpenCloseImageFilter<TInputImage, TOutputImage, THistogram, TFunction1, TFunction2>
::finishLine(OutputImagePixelType * buffer,
	     OutputImagePixelType &Extreme,
	     unsigned &outLeftP,
	     unsigned &outRightP)
{
  while (outLeftP < outRightP)
    {
    if (m_TF2(buffer[outLeftP], buffer[outRightP]))
      {
      Extreme = buffer[outRightP];
      --outRightP;
      if (!m_TF2(buffer[outRightP], Extreme))
	{
	buffer[outRightP] = Extreme;
	}
      }
    else
      {
      Extreme = buffer[outLeftP];
      ++outLeftP;
      if (!m_TF2(buffer[outLeftP], Extreme))
	{
	buffer[outLeftP] = Extreme;
	}
      }
    }
}

template<class TInputImage, class TOutputImage, class THistogram, class TFunction1, class TFunction2>
void
AnchorOpenCloseImageFilter<TInputImage, TOutputImage, THistogram, TFunction1, TFunction2>
::PrintSelf(std::ostream &os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

}


} // end namespace itk

#endif
