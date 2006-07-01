#ifndef __itkAnchorErodeDilateLine_txx
#define __itkAnchorErodeDilateLine_txx

#include "itkAnchorErodeDilateLine.h"

namespace itk {

template <class TInputPix, class TFunction1, class TFunction2>
AnchorErodeDilateLine<TInputPix, TFunction1, TFunction2>
::AnchorErodeDilateLine()
{
  m_Size=2;
  // create a histogram
  if (useVectorBasedHistogram())
    {
    m_Histo = new VHistogram;
    } 
  else
    {
    m_Histo = new MHistogram;
    }
}

template <class TInputPix, class TFunction1, class TFunction2>
void
AnchorErodeDilateLine<TInputPix, TFunction1, TFunction2>
::doLine(InputImagePixelType * buffer, InputImagePixelType * inbuffer, unsigned bufflength)
{
  // TFunction1 will be < for erosions
  // TFunction2 will be <=


  // the initial version will adopt the methodology of loading a line
  // at a time into a buffer vector, carrying out the opening or
  // closing, and then copy the result to the output. Hopefully this
  // will improve cache performance when working along non raster
  // directions.
  if (bufflength <= m_Size)
    {
    // No point doing anything fancy - just look for the extreme value
    // This is important when operating near the corner of images with
    // angled structuring elements 
    InputImagePixelType Extreme = inbuffer[0];
    for (unsigned i = 0;i < bufflength;i++) 
      {
      if (m_TF1(Extreme, inbuffer[i]))
	Extreme = inbuffer[i];
      }

    for (unsigned i = 0;i < bufflength;i++) 
      {
      buffer[i] = Extreme;
      }
    return;
    }

  unsigned int middle = m_Size/2;

  unsigned outLeftP = 0, outRightP = bufflength - 1;
  unsigned inLeftP = 0, inRightP = bufflength - 1;
  InputImagePixelType Extreme;
  m_Histo->Reset();

  // Left border, first half of structuring element
  Extreme = inbuffer[inLeftP];
  m_Histo->AddPixel(Extreme);
  for (unsigned i = 0; (i < middle); i++)
    {
    ++inLeftP;
    m_Histo->AddPixel(inbuffer[inLeftP]);
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
    m_Histo->AddPixel(inbuffer[inLeftP]);
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
    if ((inLeftP - m_Size) < 0) 
      {
      }
    m_Histo->RemovePixel(inbuffer[inLeftP - m_Size]);
    m_Histo->AddPixel(inbuffer[inLeftP]);
    Extreme = m_Histo->GetValue();
    buffer[outLeftP] = Extreme;
    }
  Extreme = buffer[outLeftP];

  while (startLine(buffer, inbuffer, Extreme, *m_Histo, outLeftP, outRightP, inLeftP, inRightP, middle)){}

  finishLine(buffer, inbuffer, Extreme, *m_Histo, outLeftP, outRightP, inLeftP, inRightP, middle);

}

template<class TInputPix, class TFunction1, class TFunction2>
bool
AnchorErodeDilateLine<TInputPix, TFunction1, TFunction2>
::startLine(InputImagePixelType * buffer,
	    InputImagePixelType * inbuffer,
	    InputImagePixelType &Extreme,
	    Histogram &histo,
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

template<class TInputPix, class TFunction1, class TFunction2>
bool
AnchorErodeDilateLine<TInputPix, TFunction1, TFunction2>
::finishLine(InputImagePixelType * buffer,
	     InputImagePixelType * inbuffer,
	     InputImagePixelType &Extreme,
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

template<class TInputPix, class TFunction1, class TFunction2>
void
AnchorErodeDilateLine<TInputPix, TFunction1, TFunction2>
::PrintSelf(std::ostream &os, Indent indent) const
{
  os << indent << "Size: " << m_Size << std::endl;
}


} // end namespace itk

#endif
