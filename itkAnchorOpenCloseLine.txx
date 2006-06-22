#ifndef __itkAnchorOpenCloseLine_txx
#define __itkAnchorOpenCloseLine_txx

#include "itkAnchorOpenCloseLine.h"

namespace itk {

template <class TInputPix, class THistogramCompare, class TFunction1, class TFunction2>
AnchorOpenCloseLine<TInputPix, THistogramCompare, TFunction1, TFunction2>
::AnchorOpenCloseLine()
{
  m_Size=2;
  if (useVectorBasedHistogram())
    {
    m_Histo = new VHistogram;
    } 
  else
    {
    m_Histo = new MHistogram;
    }
}

template <class TInputPix, class THistogramCompare, class TFunction1, class TFunction2>
void
AnchorOpenCloseLine<TInputPix, THistogramCompare, TFunction1, TFunction2>
::doLine(InputImagePixelType * buffer, unsigned bufflength)
{
  // TFunction1 will be >= for openings
  // TFunction2 will be <=
  // TFunction3 will be >


  // the initial version will adopt the methodology of loading a line
  // at a time into a buffer vector, carrying out the opening or
  // closing, and then copy the result to the output. Hopefully this
  // will improve cache performance when working along non raster
  // directions.


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
  InputImagePixelType Extreme;
  while (startLine(buffer, Extreme, *m_Histo, outLeftP, outRightP)){}
  
  finishLine(buffer, Extreme, outLeftP, outRightP);
}

template<class TInputPix, class THistogramCompare, class TFunction1, class TFunction2>
bool
AnchorOpenCloseLine<TInputPix, THistogramCompare, TFunction1, TFunction2>
::startLine(InputImagePixelType * buffer,
	    InputImagePixelType &Extreme,
	    Histogram &histo,
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
      endP = currentP;
      for (unsigned PP = outLeftP + 1; PP < endP; ++PP)
	{
	buffer[PP] = Extreme;
	}
      outLeftP = currentP;
      return (true);
      }
    ++currentP;
    }
  // We didn't find a smaller (for opening) value in the segment of
  // reach of outLeftP. currentP is the first position outside the
  // reach of outLeftP
  if (m_TF2(buffer[currentP], Extreme))
    {
    endP = currentP;
    for (unsigned PP = outLeftP + 1; PP < endP; ++PP)
      {
      buffer[PP] = Extreme;
      }
    outLeftP = currentP;
    return(true);
    }
  else
    {
    // Now we need a histogram
    // Initialise it
    histo.Reset();
    ++outLeftP;
    for (unsigned aux = outLeftP; aux <= currentP; ++aux)
      {
      histo.AddPixel(buffer[aux]);
      }
    // find the minimum value. The version
    // in the paper assumes integer pixel types and initializes the
    // search to the current extreme. Hopefully the latter is an
    // optimization step.
    Extreme = histo.GetValue();
    histo.RemovePixel(buffer[outLeftP]);
    buffer[outLeftP] = Extreme;
    histo.AddPixel(Extreme);
    }

  while (currentP < outRightP)
    {
    ++currentP;
    if (m_TF2(buffer[currentP], Extreme))
      {
      // Found a new extrem
      endP = currentP;
      for (unsigned PP = outLeftP + 1; PP < endP; PP++)
	{
	buffer[PP]=Extreme;
	}
      outLeftP = currentP;
      return(true);
      }
    else
      {
      /* histogram update */
      histo.AddPixel(buffer[currentP]);
      histo.RemovePixel(buffer[outLeftP]);
      Extreme = histo.GetValue();
      ++outLeftP;
      histo.RemovePixel(buffer[outLeftP]);
      buffer[outLeftP] = Extreme;
      histo.AddPixel(Extreme);
      }
    }
  // Finish the line
  while (outLeftP < outRightP)
    {
    histo.RemovePixel(buffer[outLeftP]);
    Extreme = histo.GetValue();
    ++outLeftP;
    histo.RemovePixel(buffer[outLeftP]);
    buffer[outLeftP] = Extreme;
    histo.AddPixel(Extreme);
    }
  return(false);
}

template<class TInputPix, class THistogramCompare, class TFunction1, class TFunction2>
bool
AnchorOpenCloseLine<TInputPix,  THistogramCompare, TFunction1, TFunction2>
::finishLine(InputImagePixelType * buffer,
	     InputImagePixelType &Extreme,
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

template<class TInputPix, class THistogramCompare, class TFunction1, class TFunction2>
void
AnchorOpenCloseLine<TInputPix, THistogramCompare, TFunction1, TFunction2>
::PrintSelf(std::ostream &os, Indent indent) const
{
  os << indent << "Size: " << m_Size << std::endl;
}


} // end namespace itk

#endif
