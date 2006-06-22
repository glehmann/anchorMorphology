#include "itkBresenhamLine.h"

int main(int, char * argv[])
{
  const unsigned int dim = 2;
  
  typedef itk::BresenhamLine<dim> BresType;
  BresType::OffsetArray A;
  BresType B;
  BresType::LType dir;
  int len = 100;
  dir[0] = 1;
  dir[1] = 0.25;
  dir.Normalize();
  A = B.buildLine(dir, len);

  std::cout << "Direction = " << dir << std::endl;
  std::cout << "Result length = " << A.size() << std::endl;
  
  for (unsigned i = 0;i<A.size();i++)
    {
    std::cout << A[i] << std::endl;
    }
  return 0;
}

