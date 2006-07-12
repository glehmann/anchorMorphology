#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkCommand.h"
#include "itkSimpleFilterWatcher.h"
#include "itkFlatStructuringElement.h"

#include "itkAnchorDilateImageFilter.h"

// test the decomposition by dilating a single pixel

template <class TImage>
typename TImage::Pointer createIm3D(unsigned nx, unsigned ny, unsigned nz)
{
  typename TImage::Pointer Result = TImage::New();
  typename TImage::IndexType start;
  typename TImage::SizeType size;
  typename TImage::RegionType region;
  start.Fill(0);
  size[0]=nx;
  size[1]=ny;
  size[2]=nz;
  region.SetSize(size);
  region.SetIndex(start);
  Result->SetRegions(region);
  Result->Allocate();
  return(Result);
  
}

int main(int, char * argv[])
{
  const int dim = 3;
  
  typedef unsigned char PType;
  typedef itk::Image< PType, dim > IType;

  // create an image
  IType::Pointer inimage = createIm3D<IType>(201, 201, 201);
  IType::IndexType where;
  where[0]=100;
  where[1]=100;
  where[2]=100;
  inimage->FillBuffer(0);
  inimage->SetPixel(where, 255);

  typedef itk::FlatStructuringElement<dim> SEType;

  typedef itk::AnchorDilateImageFilter<IType, SEType > FilterType;

  FilterType::Pointer filter = FilterType::New();
  filter->SetInput( inimage );
  SEType::RadiusType Rad;
  int lines = atoi(argv[1]);
  Rad[0]=atoi(argv[2]);
  Rad[1]=atoi(argv[3]);
  Rad[2]=atoi(argv[4]);
  SEType K = SEType::Poly(Rad, lines);
  //SEType K = SEType::Box(Rad);

  filter->SetKernel(K);
  itk::SimpleFilterWatcher watcher(filter, "filter");

  typedef itk::ImageFileWriter< IType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( filter->GetOutput() );
  writer->SetFileName( argv[5] );
  writer->Update();

  return 0;
}

