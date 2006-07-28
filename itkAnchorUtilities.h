#ifndef __itkAnchorUtilities_h
#define __itkAnchorUtilities_h

#include <list>

#define ANCHOR_ALGORITHM

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
		  const TLine line);

#ifdef ANCHOR_ALGORITHM
template <class TImage, class TBres, class TLine>
int fillLineBuffer(typename TImage::ConstPointer input,
		   const typename TImage::IndexType StartIndex,
		   const TLine line,
		   const float tol,
		   const typename TBres::OffsetArray LineOffsets,
		   const typename TImage::RegionType AllImage, 
		   typename TImage::PixelType * inbuffer,
		   unsigned &start,
		   unsigned &end);
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
		   unsigned &end);
#endif
template <class TImage, class TBres, class TLine>
int computeStartEnd(const typename TImage::IndexType StartIndex,
		    const TLine line,
		    const float tol,
		    const typename TBres::OffsetArray LineOffsets,
		    const typename TImage::RegionType AllImage, 
		    unsigned &start,
		    unsigned &end);
#if 0
template <class TImage, class TBres>
void fillLineBuffer(typename TImage::ConstPointer input,
		    const typename TImage::IndexType StartIndex,
		    const typename TBres::OffsetArray LineOffsets,
		    const typename TImage::RegionType AllImage, 
		    typename TImage::PixelType * inbuffer,
		    unsigned &len);
#endif

template <class TImage, class TBres>
void copyLineToImage(const typename TImage::Pointer output,
		     const typename TImage::IndexType StartIndex,
		     const typename TBres::OffsetArray LineOffsets,
		     const typename TImage::PixelType * outbuffer,
		     const unsigned start,
		     const unsigned end);

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
	    const typename TImage::RegionType face);
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
	    const typename TImage::RegionType face);

#endif
// This creates a list of non overlapping faces that need to be
// processed for this particular line orientation. We are doing this
// instead of using the Face Calculator to avoid repeated operations
// starting from corners.
template <class TRegion, class TLine>
std::list<TRegion> mkFaceList(const TRegion AllImage,
			      const TLine line);

// This returns a face with a normal between +/- 45 degrees of the
// line. The face is enlarged so that AllImage is entirely filled by
// lines starting from every pixel in the face. This means that some
// of the region will not touch the image. This approach is necessary
// because we want to be able to sweep the lines in a fashion that
// does not have overlap between them.
template <class TInputImage, class TLine>
typename TInputImage::RegionType
mkEnlargedFace(const typename TInputImage::ConstPointer input,
	       const typename TInputImage::RegionType AllImage,
	       const TLine line);

// figure out the correction factor for length->pixel count based on
// line angle
template <class TLine>
unsigned int getLinePixels(const TLine line);

} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkAnchorUtilities.txx"
#endif

#endif
