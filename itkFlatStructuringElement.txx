
#ifndef __itkFlatStructuringElement_txx
#define __itkFlatStructuringElement_txx

#include "itkFlatStructuringElement.h"
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


namespace itk
{
template<unsigned int VDimension>
FlatStructuringElement<VDimension> FlatStructuringElement<VDimension>
::Poly(RadiusType radius, unsigned lines)
{
    FlatStructuringElement res = FlatStructuringElement();
    res = res.PolySub(Dispatch<VDimension>(), radius, lines);
#if 0
    float theta, phi, step;
    theta = phi = 0;
    step = M_PI/(lines - 1);

    while (theta < M_PI) 
      {

      std::cout << "theta= " << theta << " phi = " << phi << std::endl;
      LType O = res.mkOffset(phi, theta);
      std::cout << O << std::endl;
      if (res.checkParallel(O, res.m_Lines))
	{
	std::cout << "Already have this line" << std::endl;
	}
      else
	{
	res.m_Lines.push_back(O);
	}
      theta += step;
      phi += step;
      }
    std::cout << "---------------" << std::endl;
#endif
    return(res);

}

template<unsigned int VDimension>
FlatStructuringElement<VDimension> FlatStructuringElement<VDimension>
::PolySub(const Dispatch<2> &, RadiusType radius, unsigned lines) const
{
  // radial decomposition method from "Radial Decomposition of Discs
  // and Spheres" - CVGIP: Graphical Models and Image Processing
  //std::cout << "2 dimensions" << std::endl;
  FlatStructuringElement res = FlatStructuringElement();
  res.m_Decomposition = true;
  
  unsigned int rr = 0;
  for (unsigned i=0;i<VDimension;i++)
    {
    if (radius[i] > rr) rr = radius[i];
    }
  if (lines == 0)
    {
    // select some default line values
    if (rr <= 3) lines=2;
    else if (rr <= 8) lines=4;
    else lines=6;
    }
  // start with a circle - figure out the length of the structuring
  // element we need -- This method results in a polygon with 2*lines
  // sides, each side with length k, where k is the structuring
  // element length. Therefore the value of k we need to produce the
  // radius we want is: (M_PI * rr * 2)/(2*lines*2)
  float k = (M_PI * (float)rr)/((float)lines);
  std::cout << "k= " << k << std::endl;
  float theta, step;
  theta = 0;
  step = M_PI/lines;
  while (theta < M_PI) 
    {
//    std::cout << "theta= " << theta << std::endl;
    LType O;
    O[0]= k * cos(theta);
    O[1]= k * sin(theta);
    if (res.checkParallel(O, res.m_Lines))
      {
      std::cout << "Already have this line" << std::endl;
      }
    else
      {
      std::cout << O << std::endl;
      res.m_Lines.push_back(O);
      }
    theta += step;
    }
    
  return(res);
}

template<unsigned int VDimension>
FlatStructuringElement<VDimension> FlatStructuringElement<VDimension>
::PolySub(const Dispatch<3> &, RadiusType radius, unsigned lines) const
{
  std::cout << "3 dimensions" << std::endl;

}

template<unsigned int VDimension>
FlatStructuringElement<VDimension> FlatStructuringElement<VDimension>
::PolySub(const DispatchBase &, RadiusType radius, unsigned lines) const
{
  //itkWarningMacro("Don't know how to deal with this many dimensions");
}

template<unsigned int VDimension>
FlatStructuringElement<VDimension> FlatStructuringElement<VDimension>
::Box(RadiusType radius)
{
  // this should work for any number of dimensions
  FlatStructuringElement res = FlatStructuringElement();
  res.m_Decomposition = true;
  for (unsigned i = 0;i<VDimension;i++)
    {
    if (radius[i] != 0)
      {
      LType L;
      L.Fill(0);
      L[i] = radius[i];
      res.m_Lines.push_back(L);
      }
    }
  return(res);
}

template<unsigned int VDimension>
typename FlatStructuringElement<VDimension>::LType FlatStructuringElement<VDimension>
::mkOffset(float phi, float theta)
{
  LType res;
  res[0] = cos(phi)*cos(theta);
  res[1] = cos(phi)*sin(theta);
  res[2] = sin(theta);
  return res;
}

template<unsigned int VDimension>
bool
FlatStructuringElement<VDimension>::
checkParallel(LType NewVec, DecompType Lines)
{
  LType NN = NewVec;
  NN.Normalize();
  for (unsigned i = 0; i < Lines.size(); i++)
    {
    LType LL = Lines[i];
    LL.Normalize();
    float L = NN*LL;
    if ((1.0 - fabs(L)) < 0.001) return(true);
    }
  return(false);
}

template<unsigned int VDimension>
void FlatStructuringElement< VDimension >
::PrintSelf(std::ostream &os, Indent indent) const
{
  //Superclass::PrintSelf(os, indent);
  if (m_Decomposition)
    {
    os << indent << "SE decomposition:" << std::endl;
    for (unsigned i = 0;i < m_Lines.size(); i++)
      {
      os << indent << m_Lines[i] << std::endl;
      }
    }
}

}

#endif
