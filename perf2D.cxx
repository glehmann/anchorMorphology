#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkCommand.h"
#include "itkSimpleFilterWatcher.h"
#include "itkFlatStructuringElement.h"
#include <iomanip>
#include "itkTimeProbe.h"


#include "itkAnchorDilateImageFilter.h"

// test the decomposition by dilating a single pixel

template <class TImage>
typename TImage::Pointer createIm2D(unsigned nx, unsigned ny)
{
  typename TImage::Pointer Result = TImage::New();
  typename TImage::IndexType start;
  typename TImage::SizeType size;
  typename TImage::RegionType region;
  start.Fill(0);
  size[0]=nx;
  size[1]=ny;
  region.SetSize(size);
  region.SetIndex(start);
  Result->SetRegions(region);
  Result->Allocate();
  return(Result);
  
}

int main(int, char * argv[])
{
  const int dim = 2;
  
  typedef unsigned char PType;
  typedef itk::Image< PType, dim > IType;

  // create an image
  IType::Pointer inimage = createIm2D<IType>(201, 201);
  IType::IndexType where;
  where[0]=where[1]=100;
  inimage->FillBuffer(0);
  inimage->SetPixel(where, 255);

  typedef itk::FlatStructuringElement<dim> SEType;

  typedef itk::AnchorDilateImageFilter<IType, SEType > FilterType;

  FilterType::Pointer filter = FilterType::New();
  filter->SetInput( inimage );
  SEType::RadiusType Rad;
  int lines = atoi(argv[1]);
  Rad.Fill(atoi(argv[2]));
  SEType K = SEType::Poly(Rad, lines);

  filter->SetKernel(K);
  itk::TimeProbe timer;
  for (unsigned i = 0;i < 100;i++)
    {
    timer.Start();
    filter->Update();
    timer.Stop();
    filter->Modified();
    }

  std::cout << std::setprecision(3) << argv[1] << " lines, length " << argv[2] << " time = " << timer.GetMeanTime() << std::endl;
  return 0;
}

