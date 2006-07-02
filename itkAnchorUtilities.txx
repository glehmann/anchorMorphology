#ifndef __itkAnchorUtilities_txx
#define __itkAnchorUtilities_txx

#include "itkAnchorUtilities.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionConstIterator.h"
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

template <class TImage, class TBres>
void fillLineBuffer(typename TImage::ConstPointer input,
		    const typename TImage::IndexType StartIndex,
		    const typename TBres::OffsetArray LineOffsets,
		    const typename TImage::RegionType AllImage, 
		    typename TImage::PixelType * inbuffer,
		    unsigned &len)
{
  // first figure out how long this line will be - use a binary search
  // It may be better to predict the intersection between the line and
  // the box
  unsigned int maxPos= LineOffsets.size() - 1, minPos = 0;
  unsigned int lastPos;
  //unsigned loops = 0;
  while (true)
    {
    unsigned int pos = (int)floor((double)(maxPos - minPos))/2 + minPos;
    typename TImage::IndexType Ind = StartIndex + LineOffsets[pos];
    typename TImage::IndexType IndN = StartIndex + LineOffsets[pos+1];
    bool I1 = AllImage.IsInside(Ind);
    bool I2;
    //++loops;
    if (I1) 
      {
      I2 = !AllImage.IsInside(IndN);
      }
    else
      {
      maxPos = pos - 1;
      continue;
      }

    if (I1 && I2) 
      {
      // found the border
      lastPos = pos;
      break;
      }
    else if (I1)
      {
      // still inside
      minPos = pos + 1;
      }
    else
      {
      maxPos = pos - 1;
      }
    }
  //std::cout << "Border found in " << loops << std::endl;
  // now copy into the buffer
  len = lastPos + 1;
#if 1
  for (unsigned i = 0; i <= lastPos;i++)
    {
    inbuffer[i] = input->GetPixel(StartIndex + LineOffsets[i]);
    }
#else
  typedef ImageRegionConstIterator<TImage> ItType;
  ItType it(input, AllImage);
  it.SetIndex(StartIndex);
  for (unsigned i = 0; i < lastPos;i++)
    {
    inbuffer[i]= it.Get();
    typename TImage::IndexType I = StartIndex + LineOffsets[i];
    typename TImage::OffsetType O = I - it.GetIndex();
    //it += O;
    }
#endif
}

template <class TImage, class TBres>
void copyLineToImage(const typename TImage::Pointer output,
		     const typename TImage::IndexType StartIndex,
		     const typename TBres::OffsetArray LineOffsets,
		     const typename TImage::PixelType * outbuffer,
		     const unsigned len)
{
  for (unsigned i = 0; i < len;i++)
    {
    output->SetPixel(StartIndex + LineOffsets[i], outbuffer[i]);
    }

}

template <class TImage, class TBres, class TAnchor>
void doFace(typename TImage::ConstPointer input,
	    typename TImage::Pointer output,
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
  while (!it.IsAtEnd()) 
    {
    typename TImage::IndexType Ind = it.GetIndex();
    unsigned int len;
    fillLineBuffer<TImage, TBres>(input, Ind, LineOffsets, AllImage, inbuffer, len);
    AnchorLine.doLine(outbuffer, inbuffer, len);
    copyLineToImage<TImage, TBres>(output, Ind, LineOffsets, outbuffer, len);
    ++it;
    }

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
    float SinTheta = fabs(line[i]/N);
    float CosTheta = sqrt(1 - SinTheta*SinTheta);
    float tt = std::max(SinTheta, CosTheta);
    if (tt > correction) correction=tt;
    }
  
  N *= correction;
  return (int)(N + 0.5);
}

} // namespace itk

#endif
