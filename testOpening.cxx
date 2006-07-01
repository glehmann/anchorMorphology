#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkCommand.h"
#include "itkSimpleFilterWatcher.h"
#include "itkFlatStructuringElement.h"

#include "itkAnchorOpenImageFilter.h"
#include "itkIndent.h"

int main(int, char * argv[])
{
  const int dim = 2;
  
  typedef unsigned char PType;
  typedef itk::Image< PType, dim > IType;

  typedef itk::ImageFileReader< IType > ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );

  typedef itk::FlatStructuringElement<dim> SEType;

  typedef itk::AnchorOpenImageFilter< IType, SEType > FilterType;

  FilterType::Pointer filter = FilterType::New();
  filter->SetInput( reader->GetOutput() );
  SEType::RadiusType Rad;
  Rad.Fill(11);
  SEType K = SEType::Poly(Rad, 4);
  K.PrintSelf(std::cout, itk::Indent(0));
  filter->SetKernel(K);


  itk::SimpleFilterWatcher watcher(filter, "filter");

  typedef itk::ImageFileWriter< IType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( filter->GetOutput() );
  writer->SetFileName( argv[2] );
  writer->Update();

  return 0;
}

