#ifndef __itkAnchorUtilities_h
#define __itkAnchorUtilities_h

#include <list>

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

template <class TImage, class TBres>
void fillLineBuffer(typename TImage::ConstPointer input,
		    const typename TImage::IndexType StartIndex,
		    const typename TBres::OffsetArray LineOffsets,
		    const typename TImage::RegionType AllImage, 
		    typename TImage::PixelType * inbuffer,
		    unsigned &len);

template <class TImage, class TBres>
void copyLineToImage(const typename TImage::Pointer output,
		     const typename TImage::IndexType StartIndex,
		     const typename TBres::OffsetArray LineOffsets,
		     const typename TImage::PixelType * outbuffer,
		     const unsigned len);

template <class TImage, class TBres, class TAnchor>
void doFace(typename TImage::ConstPointer input,
	    typename TImage::Pointer output,
	    TAnchor &AnchorLine,
	    const typename TBres::OffsetArray LineOffsets,
	    typename TImage::PixelType * inbuffer,
	    typename TImage::PixelType * outbuffer,	      
	    const typename TImage::RegionType AllImage, 
	    const typename TImage::RegionType face);

// This creates a list of non overlapping faces that need to be
// processed for this particular line orientation. We are doing this
// instead of using the Face Calculator to avoid repeated operations
// starting from corners.
template <class TRegion, class TLine>
std::list<TRegion> mkFaceList(const TRegion AllImage,
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
