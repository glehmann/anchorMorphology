#ifndef __itkAnchorUtilities_txx
#define __itkAnchorUtilities_txx

#include "itkAnchorUtilities.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionConstIterator.h"
#include "itkNeighborhoodAlgorithm.h"

namespace itk {

/**
 * \class AnchorUtilities 
 * \brief functionality in common for anchor openings/closings and
 * erosions/dilation
 *
**/

template <class TRegion, class TLine>
bool needToDoFace(const TRegion AllImage,
		  const TRegion face,
		  const TLine line)
{
  // can't use the continuous IsInside (even if I could get it to
  // work) because on the edge doesn't count as inside for this test

  // If the component of the vector orthogonal to the face doesn't go
  // inside the image then we can ignore the face
  
  // find the small dimension of the face - should only be one
  typename TRegion::SizeType ISz = AllImage.GetSize();
  typename TRegion::IndexType ISt = AllImage.GetIndex();

  typename TRegion::SizeType FSz = face.GetSize();
  typename TRegion::IndexType FSt = face.GetIndex();

  unsigned smallDim;
  for (unsigned i = 0; i < AllImage.GetImageDimension(); i++)
    {
    if (FSz[i] == 1)
      {
      smallDim = i;
      break;
      }
    }
  long startI = ISt[smallDim];
  long endI = (ISt[smallDim] + ISz[smallDim] - 1);
  long facePos = FSt[smallDim] + FSz[smallDim] - 1;
  if (facePos == startI) 
    {
    // at the start of dimension - vector must be positive
    if (line[smallDim] > 0.000001) return true;  
    // some small angle that we consider to be zero - should be more rigorous
    }
  else
    {    
    // at the end of dimension - vector must be positive
    if (line[smallDim] < -0.000001) return true;  
    }
  return (false);
  
}
template <class TImage, class TBres, class TLine>
int computeStartEnd(const typename TImage::IndexType StartIndex,
		    const TLine line,
		    const float tol,
		    const typename TBres::OffsetArray LineOffsets,
		    const typename TImage::RegionType AllImage, 
		    unsigned &start,
		    unsigned &end)
{
  // compute intersection between ray and box
  typename TImage::IndexType ImStart = AllImage.GetIndex();
  typename TImage::SizeType ImSize = AllImage.GetSize();
  float Tfar = NumericTraits<float>::max();
  float Tnear  = NumericTraits<float>::NonpositiveMin();
  float domdir  = NumericTraits<float>::NonpositiveMin();
  int sPos, ePos;
  unsigned perpdir;
  for (unsigned i = 0; i< TImage::RegionType::ImageDimension; i++)
    {
    if (fabs(line[i]) > domdir)
      {
      domdir = fabs(line[i]);
      perpdir = i;
      }
    if (fabs(line[i]) > tol)
      {
      int P1 = ImStart[i] - StartIndex[i];
      int P2 = ImStart[i] + ImSize[i] - 1 - StartIndex[i];
      float T1 = ((float)(P1))/line[i];
      float T2 = ((float)(P2))/line[i];
      
      if (T1 > T2)
	{
	// T1 is meant to be the near face
	std::swap(T1, T2);
	std::swap(P1, P2);
	}
      // want the farthest Tnear and nearest TFar
      if (T1 > Tnear) 
	{
	Tnear = T1;  
	sPos = abs(P1);
	}
      if (T2 < Tfar) 
	{
	Tfar = T2;
	ePos = abs(P2);
	}
      }
    else
      {
      // parallel to an axis - check for intersection at all
      if ((StartIndex[i] < ImStart[i]) || (StartIndex[i] > ImStart[i] + ImSize[i] - 1))
	{
	// no intersection
//	std::cout << StartIndex << "No intersection - parallel" << std::endl;
	start=end=0;
	return(0);
	}
      }
    }
  sPos = (int)(Tnear*fabs(line[perpdir]) + 0.5);
  ePos = (int)(Tfar*fabs(line[perpdir]) + 0.5);
  
  //std::cout << Tnear << " " << Tfar << std::endl;
  if (Tfar < Tnear) // seems to need some margin
    {
    // in theory, no intersection, but search between them
    bool intersection = false;
    unsigned inside;
    if (Tnear - Tfar < 10)
      {
//      std::cout << "Searching " << Tnear << " " << Tfar << std::endl;
      for (unsigned i = ePos; i<= sPos; i++)
	{
	if (AllImage.IsInside(StartIndex + LineOffsets[i]))
	  {
	  inside = i;
	  intersection=true;
	  break;
	  }
	}
      }
    if (intersection)
      {
//      std::cout << "Found intersection after all" << std::endl;
      sPos = ePos = inside;
      while (AllImage.IsInside(StartIndex + LineOffsets[ePos + 1])) ++ePos;
      while (AllImage.IsInside(StartIndex + LineOffsets[sPos - 1])) --sPos;
      start = sPos;
      end = ePos;
      }
    else
      {
//      std::cout << StartIndex << "No intersection" << std::endl;
      start=end=0;
      return(0);
      }
    }
  else
    {
    
    if (AllImage.IsInside(StartIndex + LineOffsets[sPos]))
      {
      for (;;)
	{
	if (!AllImage.IsInside(StartIndex + LineOffsets[sPos - 1])) break;
	else --sPos;
	}
      }
    else
      {
      for(;;)
	{
	++sPos;
	if (!AllImage.IsInside(StartIndex + LineOffsets[sPos])) ++sPos;
	else break;
	}
      }
    if (AllImage.IsInside(StartIndex + LineOffsets[ePos]))
      {
      for(;;)
	{
	if (!AllImage.IsInside(StartIndex + LineOffsets[ePos + 1])) break;
	else ++ePos;
	}
      }
    else
      {
      for (;;)
	{
	--ePos;
	if (!AllImage.IsInside(StartIndex + LineOffsets[ePos])) --ePos;
	else break;
	}
      }
    }
  start = sPos;
  end = ePos;
  return (1);
}

#ifdef ANCHOR_ALGORITHM

template <class TImage, class TBres, class TLine>
int fillLineBuffer(typename TImage::ConstPointer input,
		   const typename TImage::IndexType StartIndex,
		   const TLine line,  // unit vector
		   const float tol,
		   const typename TBres::OffsetArray LineOffsets,
		   const typename TImage::RegionType AllImage, 
		   typename TImage::PixelType * inbuffer,
		   unsigned &start,
		   unsigned &end)
{

//   if (AllImage.IsInside(StartIndex))
//     {
//     start = 0;
//     }
//   else
  
#if 0
  // we need to figure out where to start
  // this is an inefficient way we'll use for testing
  for (start=0;start < LineOffsets.size();++start)
    {
    if (AllImage.IsInside(StartIndex + LineOffsets[start])) break;
    }
#else
  int status = computeStartEnd<TImage, TBres, TLine>(StartIndex, line, tol, LineOffsets, AllImage,
						     start, end);
  if (!status) return(status);
#endif
#if 1
  unsigned size = end - start + 1;
  for (unsigned i = 0; i < size;i++)
    {
    inbuffer[i] = input->GetPixel(StartIndex + LineOffsets[start + i]);
    }
#else
  typedef ImageRegionConstIteratorWithIndex<TImage> ItType;
  ItType it(input, AllImage);
  it.SetIndex(StartIndex);
  for (unsigned i = 0; i < lastPos;i++)
    {
    inbuffer[i]= it.Get();
    typename TImage::IndexType I = StartIndex + LineOffsets[i];
    typename TImage::OffsetType Off = I - it.GetIndex();
    it += Off;
    }
#endif
  return(1);
}

#else

template <class TImage, class TBres, class TLine, class TFunction>
int fillLineBuffer(typename TImage::ConstPointer input,
		   const typename TImage::IndexType StartIndex,
		   const TLine line,  // unit vector
		   const float tol,
		   const typename TBres::OffsetArray LineOffsets,
		   const typename TImage::RegionType AllImage,
		   const unsigned int KernLen,
		   typename TImage::PixelType * pixbuffer,
		   typename TImage::PixelType * fExtBuffer,
		   unsigned &start,
		   unsigned &end)
{

 int status = computeStartEnd<TImage, TBres, TLine>(StartIndex, line, tol, LineOffsets, AllImage,
						     start, end);
  if (!status) return(status);

  unsigned size = end - start + 1;
  unsigned blocks = size/KernLen;
  unsigned i = 0;
  TFunction m_TF;
//  std::cout << "Line length = " << size << " KernSize = " << KernLen << " Blocks = " << blocks << std::endl;
  for (unsigned j = 0; j<blocks;j++)
    {
//    std::cout << "f1: i = " << i << std::endl;
    typename TImage::PixelType Ext = input->GetPixel(StartIndex + LineOffsets[start + i]);
    pixbuffer[i] = Ext;
    fExtBuffer[i]=Ext;
    ++i;
    for (unsigned k = 1; k < KernLen; k++)
      {
//      std::cout << "f2: i = " << i << std::endl;
      typename TImage::PixelType V = input->GetPixel(StartIndex + LineOffsets[start + i]);
      pixbuffer[i] = V;
      if (m_TF(V, fExtBuffer[i-1]))
	{
	fExtBuffer[i] = V;
	}
      else
	{
	fExtBuffer[i] = fExtBuffer[i-1];
	}
      ++i;
      }
    }
  // finish the rest
  if (i != size - 1)
    {
    typename TImage::PixelType V = input->GetPixel(StartIndex + LineOffsets[start + i]);
    pixbuffer[i] = V;
    fExtBuffer[i] = V;
    i++;
    }
  while (i < size)
    {
    typename TImage::PixelType V = input->GetPixel(StartIndex + LineOffsets[start + i]);
    pixbuffer[i] = V;
    if (m_TF(V, fExtBuffer[i-1]))
      {
      fExtBuffer[i] = V;
      }
    else
      {
      fExtBuffer[i] = fExtBuffer[i-1];
      }
    ++i;
    }
//   for (unsigned h = 0; h < size; h++)
//     {
//     std::cout << (int)pixbuffer[h] << " " << (int)fExtBuffer[h] << std::endl;
//     }
  return(1);
}


#endif

template <class TImage, class TBres>
void copyLineToImage(const typename TImage::Pointer output,
		     const typename TImage::IndexType StartIndex,
		     const typename TBres::OffsetArray LineOffsets,
		     const typename TImage::PixelType * outbuffer,
		     const unsigned start,
		     const unsigned end)
{
  unsigned size = end - start + 1;
  for (unsigned i = 0; i <size; i++)
    {
#if 1
    output->SetPixel(StartIndex + LineOffsets[start + i], outbuffer[i]);
#else
    output->SetPixel(StartIndex + LineOffsets[start + i], 1+outbuffer[i]);
#endif
    }
//  std::cout << "Copy out " << StartIndex << StartIndex + LineOffsets[len-1] << std::endl;
}

#ifdef ANCHOR_ALGORITHM

template <class TImage, class TBres, class TAnchor, class TLine>
void doFace(typename TImage::ConstPointer input,
	    typename TImage::Pointer output,
	    TLine line,
	    TAnchor &AnchorLine,
	    const typename TBres::OffsetArray LineOffsets,
	    typename TImage::PixelType * inbuffer,
	    typename TImage::PixelType * outbuffer,	      
	    const typename TImage::RegionType AllImage, 
	    const typename TImage::RegionType face)
{
  // iterate over the face
  typedef ImageRegionConstIteratorWithIndex<TImage> ItType;
  ItType it(input, face);
  it.GoToBegin();
  TLine NormLine = line;
  NormLine.Normalize();
  // set a generous tolerance
  float tol = 1.0/LineOffsets.size();
  while (!it.IsAtEnd()) 
    {
    typename TImage::IndexType Ind = it.GetIndex();
    unsigned start, end, len;
    if (fillLineBuffer<TImage, TBres, TLine>(input, Ind, NormLine, tol, LineOffsets, 
					     AllImage, inbuffer, start, end))
      {
      len = end - start + 1;
#if 1
      AnchorLine.doLine(outbuffer, inbuffer, len);
      copyLineToImage<TImage, TBres>(output, Ind, LineOffsets, outbuffer, start, end);
#else
      // test the decomposition
      copyLineToImage<TImage, TBres>(output, Ind, LineOffsets, inbuffer, start, end);
#endif
      }
    ++it;
    }

}

#else
template <class TImage, class TBres, class TFunction, class TLine>
void doFace(typename TImage::ConstPointer input,
	    typename TImage::Pointer output,
	    TLine line,
	    const typename TBres::OffsetArray LineOffsets,
	    const unsigned int KernLen,
	    typename TImage::PixelType * pixbuffer,
	    typename TImage::PixelType * fExtBuffer,	      
	    typename TImage::PixelType * rExtBuffer,	      
	    const typename TImage::RegionType AllImage, 
	    const typename TImage::RegionType face)
{
  // iterate over the face
  typedef ImageRegionConstIteratorWithIndex<TImage> ItType;
  ItType it(input, face);
  it.GoToBegin();
  TLine NormLine = line;
  NormLine.Normalize();
  // set a generous tolerance
  float tol = 1.0/LineOffsets.size();

  while (!it.IsAtEnd()) 
    {
    typename TImage::IndexType Ind = it.GetIndex();
    unsigned start, end, len;
    if (fillLineBuffer<TImage, TBres, TLine, TFunction>(input, Ind, NormLine, tol, LineOffsets, 
							AllImage, KernLen, pixbuffer, fExtBuffer, 
							start, end))
      {
      len = end - start + 1;
      // compute the reverse running extreme -- not that we aren't
      // using unsigned types because it is messy using them in
      // reverse loops
      int size = (int)(end - start + 1);
      int blocks = size/(int)KernLen;
      int i = size - 1;
      TFunction m_TF;

      if ((i > (blocks * (int)KernLen - 1)))
	{
	rExtBuffer[i] = pixbuffer[i];
//	std::cout << "r1: i = " << i << " " << (int) rExtBuffer[i] << std::endl;
	--i;
	while (i >= (int)(blocks * KernLen))
	  {
//	  std::cout << "r2: i = " << i << std::endl;
	  
	  typename TImage::PixelType V = pixbuffer[i];
	  if (m_TF(V, rExtBuffer[i+1]))
	    {
	    rExtBuffer[i] = V;
	    }
	  else
	    {
	    rExtBuffer[i] = rExtBuffer[i+1];
	    }
//	  std::cout << "r2: i = " << i << " " << (int) rExtBuffer[i] << std::endl;
	  --i;
	  }	
	}
      for (unsigned j = 0; j<blocks;j++)
	{
//	std::cout << "r3: i = " << i << std::endl;
	typename TImage::PixelType Ext = pixbuffer[i];
	rExtBuffer[i]=Ext;
	--i;
	for (unsigned k = 1; k < KernLen; k++)
	  {
//	  std::cout << "r4: i = " << i << std::endl;
	  typename TImage::PixelType V = pixbuffer[i];
	  if (m_TF(V, rExtBuffer[i+1]))
	    {
	    rExtBuffer[i] = V;
	    }
	  else
	    {
	    rExtBuffer[i] = rExtBuffer[i+1];
	    }
	  --i;
	  }
	}
      // now compute result
      if (size <= KernLen/2)
	{
	for (unsigned j = 0;j < size;j++)
	  {
	  pixbuffer[j] = fExtBuffer[size-1];
	  }
	}
      else if (size <= KernLen)
	{
	for (unsigned j = 0;j < size - KernLen/2;j++)
	  {
	  pixbuffer[j] = fExtBuffer[j + KernLen/2];
	  }
	for (unsigned j =  size - KernLen/2; j <= KernLen/2; j++)
	  {
	  pixbuffer[j] = fExtBuffer[size-1];
	  }
	for (unsigned j =  KernLen/2 + 1; j < size; j++)
	  {
	  pixbuffer[j] = rExtBuffer[j - KernLen/2];
	  }
	}
      else
	{
	// line beginning
	for (unsigned j = 0;j < KernLen/2;j++)
	  {
	  pixbuffer[j] = fExtBuffer[j + KernLen/2];
	  }
	for (unsigned j = KernLen/2, k=KernLen/2 + KernLen/2, l = KernLen/2 - KernLen/2;
	     j < size - KernLen/2; j++, k++, l++)
	  {
	  typename TImage::PixelType V1 = fExtBuffer[k];
	  typename TImage::PixelType V2 = rExtBuffer[l];
	  if (m_TF(V1, V2))
	    {
	    pixbuffer[j] = V1;
	    }
	  else
	    {
	    pixbuffer[j] = V2;
	    }
	  }
	// line end -- involves reseting the end of the reverse
	// extreme array
	for (unsigned j = size - 2; (j > 0) && (j >= (size - KernLen - 1)); j--)
	  {
//	  std::cout << "j1 = " << j << " ";
	  if (m_TF(rExtBuffer[j+1], rExtBuffer[j]))
	    {
	    rExtBuffer[j] = rExtBuffer[j+1];
	    }
//	    std::cout << (int) rExtBuffer[j] << " ";
	  }
//	std::cout << std::endl;
	for (unsigned j = size - KernLen/2; j < size;j++)
	  {
//	  std::cout << "j2 = " << j - KernLen/2 << std::endl;
	  pixbuffer[j]=rExtBuffer[j-KernLen/2];
	  }
	}
      copyLineToImage<TImage, TBres>(output, Ind, LineOffsets, pixbuffer, start, end);
      }
    ++it;
    }
}
#endif
template <class TInputImage, class TLine>
typename TInputImage::RegionType
mkEnlargedFace(const typename TInputImage::ConstPointer input,
	       const typename TInputImage::RegionType AllImage,
	       const TLine line)
{
  // use the face calculator to produce a face list
  typedef itk::NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<TInputImage>
FaceCalculatorType;
  FaceCalculatorType faceCalculator;
  typename TInputImage::SizeType radius;
  radius.Fill(1);
  typename FaceCalculatorType::FaceListType faceList;
  faceList = faceCalculator(input, AllImage, radius);
  typename FaceCalculatorType::FaceListType::iterator fit;
  fit = faceList.begin();
  typename TInputImage::RegionType RelevantRegion;
  bool foundFace = false;
  float MaxComp = NumericTraits<float>::NonpositiveMin();
  unsigned DomDir;
  ++fit;
//  std::cout << "------------" << std::endl;
  // figure out the dominant direction of the line
  for (unsigned i = 0;i< TInputImage::RegionType::ImageDimension;i++) 
    {
    if (fabs(line[i]) > MaxComp)
      {
      MaxComp = fabs(line[i]);
      DomDir = i;
      }
    }
  for (;fit != faceList.end();++fit) 
    {
    // check whether this face is suitable for parallel sweeping - i.e
    // whether the line is within 45 degrees of the perpendicular
    // Figure out the perpendicular using the region size
    unsigned FaceDir;
    for (unsigned i = 0;i< TInputImage::RegionType::ImageDimension;i++) 
      {
      if (fit->GetSize()[i] == 1) FaceDir = i;
      }
    if (FaceDir == DomDir) // within 1 degree 
      {
      // now check whether the line goes inside the image from this face
      if ( needToDoFace<typename TInputImage::RegionType, TLine>(AllImage, *fit, line) ) 
	{
//	std::cout << "Using face: " << *fit << line << std::endl;
	RelevantRegion = *fit;
	foundFace = true;
	break;
	}
      }
    }
  if (foundFace) 
    {
    // enlarge the region so that sweeping the line across it will
    // cause all pixels to be visited.
    // find the dimension not within the face
    unsigned NonFaceDim;
    
    for (unsigned i = 0; i < TInputImage::RegionType::ImageDimension;i++) 
      {
      if (RelevantRegion.GetSize()[i] == 1)
	{
	NonFaceDim=i;
	break;
	}
      }

    // figure out how much extra each other dimension needs to be extended
    typename TInputImage::SizeType NewSize = RelevantRegion.GetSize();
    typename TInputImage::IndexType NewStart = RelevantRegion.GetIndex();
    unsigned NonFaceLen = AllImage.GetSize()[NonFaceDim];
    for (unsigned i = 0; i < TInputImage::RegionType::ImageDimension;i++) 
      {
      if (i != NonFaceDim)
	{
	int Pad = (int)ceil((float)(NonFaceLen) * line[i]/fabs(line[NonFaceDim])); 
//	int Pad = (int)((float)(NonFaceLen) * line[i]/fabs(line[NonFaceDim])); 
	if (Pad < 0)
	  {
	  // just increase the size - no need to change the start
	  NewSize[i] += abs(Pad) + 1;
	  }
	else 
	  {
	  // change the size and index
	  NewSize[i] += Pad + 1;
	  NewStart[i] -= Pad + 1;
	  }
	}
      }
    RelevantRegion.SetSize(NewSize);
    RelevantRegion.SetIndex(NewStart);
    }
  else 
    {
    std::cout << "Line " << line << " doesnt correspond to a face" << std::endl;
    }
//  std::cout << "Result region = " << RelevantRegion << std::endl;

//  std::cout << "+++++++++++++++++" << std::endl;
  return RelevantRegion;
}

template <class TRegion, class TLine>
std::list<TRegion> mkFaceList(const TRegion AllImage,
			      const TLine line)
{
  std::list<TRegion> resultList;
  TRegion AllImageC = AllImage;
  typename TRegion::IndexType ImStart = AllImageC.GetIndex();
  typename TRegion::SizeType ImSize = AllImageC.GetSize();
  typename TRegion::IndexType FStart;
  typename TRegion::SizeType FSize;

  for (unsigned i = 0; i < TRegion::ImageDimension; ++i)
    {
    // check the low coordinate end
    for (unsigned j = 0; j < TRegion::ImageDimension; ++j) 
      {
      FStart[j] = ImStart[j];
      if ( i == j )
	{
	FSize[j]=1;
	}
      else
	{
	FSize[j]=ImSize[j];
	}
      }
    TRegion NewReg;
    NewReg.SetSize(FSize);
    NewReg.SetIndex(FStart);
    if ( needToDoFace<TRegion, TLine>(AllImageC, NewReg, line) ) 
      {
      resultList.push_back(NewReg);
      // modify the input image region
      ImStart[i] += 1;
      ImSize[i] -= 1;
      AllImageC.SetIndex(ImStart);
      AllImageC.SetSize(ImSize);
      }
    // check the high coordinate end
    for (unsigned j = 0; j < TRegion::ImageDimension; ++j) 
      {
      if ( i == j )
	{
	FStart[j] = ImStart[j] + ImSize[j] - 1;
	FSize[j]=1;
	}
      else
	{
	FStart[j] = ImStart[j];
	FSize[j]=ImSize[j];
	}
      }
    NewReg.SetSize(FSize);
    NewReg.SetIndex(FStart);
    if ( needToDoFace<TRegion, TLine>(AllImageC, NewReg, line) ) 
      {
      resultList.push_back(NewReg);
      // modify the input image region
      ImSize[i] -= 1;
      AllImageC.SetIndex(ImStart);
      AllImageC.SetSize(ImSize);
      }
    }
  return (resultList);
}

template <class TLine>
unsigned int getLinePixels(const TLine line)
{
  float N = line.GetNorm();
  float correction = 0.0;
  
  for (unsigned int i = 0; i < TLine::Dimension; i++)
    {
    float tt = fabs(line[i]/N);
    if (tt > correction) correction=tt;
    }
  
  N *= correction;
  return (int)(N + 0.5);
}

} // namespace itk

#endif
