
#ifndef __itkFlatStructuringElement_txx
#define __itkFlatStructuringElement_txx

#include "itkFlatStructuringElement.h"
#include <math.h>
#include <vector>

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
  res.m_Decomposable = true;
  
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
  // radius we want is: (M_PI * rr * 2)/(2*lines)
  float k1 = (M_PI * (float)radius[0])/((float)lines);
  float k2 = (M_PI * (float)radius[1])/((float)lines);
  //std::cout << "k= " << k << std::endl;
  float theta, step;
  step = M_PI/lines;
  theta = 0;
  // just to ensure that we get the last one
  while (theta <= M_PI/2.0 + 0.0001)
    {
    LType O;
    O[0] = k1 * cos(theta);
    O[1] = k2 * sin(theta);
    if (!res.checkParallel(O, res.m_Lines))
      {
      //std::cout << O << std::endl;
      res.m_Lines.push_back(O);
      }
     O[0] = k1 * cos(-theta);
     O[1] = k2 * sin(-theta);
     if (!res.checkParallel(O, res.m_Lines))
       {
       //std::cout << O << std::endl;
       res.m_Lines.push_back(O);
       }
     theta += step;
     //std::cout << "theta1 = " << theta << " " << M_PI/2.0 << std::endl;
    }
    
  return(res);
}

//    O[0] = k1 * cos(phi) * cos(theta);
//    O[1] = k2 * cos(phi) * sin(theta);
//    O[2] = k3 * sin(theta);

template<unsigned int VDimension>
FlatStructuringElement<VDimension> FlatStructuringElement<VDimension>
::PolySub(const Dispatch<3> &, RadiusType radius, unsigned lines) const
{
  FlatStructuringElement res = FlatStructuringElement();
  res.m_Decomposable = true;
  std::cout << "3 dimensions" << std::endl;
  unsigned int rr = 0;
  int iterations = 1;
  int faces = lines * 2;
  for (unsigned i=0;i<VDimension;i++)
    {
    if (radius[i] > rr) rr = radius[i];
    }
  switch (faces)
    {
    case 12:
    {
    // dodecahedron
    float phi = (1.0 + sqrt(5.0)) / 2.0;
    float b = 1.0/phi;
    float c = 2.0 - phi;
    unsigned facets = 12;
    typedef std::vector<FacetType> FacetArrayType;
    FacetArrayType FacetArray = FacetArrayType(facets);
    // set up vectors normal to the faces - only put in 3 points for
    // each face:
    // face 1
    LType PP;
    FacetType Fc;
    b /= 2.0;
    c /= 2.0;

    PP[0]=c;PP[1]=0;PP[2]=0.5;
    Fc.P1 = PP;
    PP[0]=-c;PP[1]=0;PP[2]=0.5;
    Fc.P2 = PP;
    PP[0]=-b;PP[1]=b;PP[2]=b;
    Fc.P3 = PP;
    FacetArray[0] = Fc;

    PP[0]=-c;PP[1]=0;PP[2]=0.5;
    Fc.P1 = PP;
    PP[0]=c;PP[1]=0;PP[2]=0.5;
    Fc.P2 = PP;
    PP[0]=b;PP[1]=-b;PP[2]=b;
    Fc.P3 = PP;
    FacetArray[1] = Fc;
    
    PP[0]=c;PP[1]=0;PP[2]=-0.5;
    Fc.P1 = PP;
    PP[0]=-c;PP[1]=0;PP[2]=-0.5;
    Fc.P2 = PP;
    PP[0]=-b;PP[1]=-b;PP[2]=-b;
    Fc.P3 = PP;
    FacetArray[2] = Fc;

    PP[0]=-c;PP[1]=0;PP[2]=-0.5;
    Fc.P1 = PP;
    PP[0]=c;PP[1]=0;PP[2]=-0.5;
    Fc.P2 = PP;
    PP[0]=b;PP[1]=b;PP[2]=-b;
    Fc.P3 = PP;
    FacetArray[3] = Fc;

    PP[0]=0;PP[1]=0.5;PP[2]=-c;
    Fc.P1 = PP;
    PP[0]=0;PP[1]=0.5;PP[2]=c;
    Fc.P2 = PP;
    PP[0]=b;PP[1]=b;PP[2]=b;
    Fc.P3 = PP;
    FacetArray[4] = Fc;

    PP[0]=0;PP[1]=0.5;PP[2]=c;
    Fc.P1 = PP;
    PP[0]=0;PP[1]=0.5;PP[2]=-c;
    Fc.P2 = PP;
    PP[0]=-b;PP[1]=b;PP[2]=-b;
    Fc.P3 = PP;
    FacetArray[5] = Fc;

    PP[0]=0;PP[1]=-0.5;PP[2]=-c;
    Fc.P1 = PP;
    PP[0]=0;PP[1]=-0.5;PP[2]=c;
    Fc.P2 = PP;
    PP[0]=-b;PP[1]=-b;PP[2]=b;
    Fc.P3 = PP;
    FacetArray[6] = Fc;

    PP[0]=0;PP[1]=-0.5;PP[2]=c;
    Fc.P1 = PP;
    PP[0]=0;PP[1]=-0.5;PP[2]=-c;
    Fc.P2 = PP;
    PP[0]=b;PP[1]=-b;PP[2]=-b;
    Fc.P3 = PP;
    FacetArray[7] = Fc;

    PP[0]=0.5;PP[1]=c;PP[2]=0;
    Fc.P1 = PP;
    PP[0]=0.5;PP[1]=-c;PP[2]=0;
    Fc.P2 = PP;
    PP[0]=b;PP[1]=-b;PP[2]=b;
    Fc.P3 = PP;
    FacetArray[8] = Fc;

    PP[0]=0.5;PP[1]=-c;PP[2]=0;
    Fc.P1 = PP;
    PP[0]=0.5;PP[1]=c;PP[2]=0;
    Fc.P2 = PP;
    PP[0]=b;PP[1]=b;PP[2]=-b;
    Fc.P3 = PP;
    FacetArray[9] = Fc;

    PP[0]=-0.5;PP[1]=c;PP[2]=0;
    Fc.P1 = PP;
    PP[0]=-0.5;PP[1]=-c;PP[2]=0;
    Fc.P2 = PP;
    PP[0]=-b;PP[1]=-b;PP[2]=-b;
    Fc.P3 = PP;
    FacetArray[10] = Fc;

    PP[0]=-0.5;PP[1]=-c;PP[2]=0;
    Fc.P1 = PP;
    PP[0]=-0.5;PP[1]=c;PP[2]=0;
    Fc.P2 = PP;
    PP[0]=-b;PP[1]=b;PP[2]=b;
    Fc.P3 = PP;
    FacetArray[11] = Fc;
    for (unsigned j = 0; j < facets; j++)
      {
      // Find a line perpendicular to each face
      LType L, A, B;
      A = FacetArray[j].P2 - FacetArray[j].P1;
      B = FacetArray[j].P3 - FacetArray[j].P1;
      L[0] = A[1]*B[2] - B[1]*A[2];
      L[1] = B[0]*A[2] - A[0]*B[2];
      L[2] = A[0]*B[1] - B[0]*A[1];
      
      L.Normalize();
      // Scale to required length
      L*= rr;
      if (!res.checkParallel(L, res.m_Lines))
	{
	res.m_Lines.push_back(L);
	}
      }
    return(res);

    }
    break;
    case 14:
    {
    // cube with the corners cut off
    LType A;
    // The axes
    A[0]=1;A[1]=0;A[2]=0;
    A *= rr;
    res.m_Lines.push_back(A);
    A[0]=0;A[1]=1;A[2]=0;
    A *= rr;
    res.m_Lines.push_back(A);
    A[0]=0;A[1]=0;A[2]=1;
    A *= rr;
    res.m_Lines.push_back(A);
    // Diagonals
    A[0]=1;A[1]=1;A[2]=1;
    A.Normalize();
    A *= rr;
    res.m_Lines.push_back(A);

    A[0]=-1;A[1]=1;A[2]=1;
    A.Normalize();
    A *= rr;
    res.m_Lines.push_back(A);

    A[0]=1;A[1]=-1;A[2]=1;
    A.Normalize();
    A *= rr;
    res.m_Lines.push_back(A);

    A[0]=-1;A[1]=-1;A[2]=1;
    A.Normalize();
    A *= rr;
    res.m_Lines.push_back(A);
    return(res);
    }
    break;
    case 20:
    {
    // Icosahedron
    float phi = (1.0 + sqrt(5.0)) / 2.0;
    float a = 0.5;
    float b = 1.0/(2.0*phi);
    unsigned facets = 20;
    typedef std::vector<FacetType> FacetArrayType;
    FacetArrayType FacetArray = FacetArrayType(facets);
    // set up vectors normal to the faces - only put in 3 points for
    // each face:
    // face 1
    LType PP;
    FacetType Fc;

    PP[0]=0;PP[1]=b;PP[2]=-a;
    Fc.P1 = PP;
    PP[0]=b;PP[1]=a;PP[2]=0;
    Fc.P2 = PP;
    PP[0]=-b;PP[1]=a;PP[2]=0;
    Fc.P3 = PP;
    FacetArray[0] = Fc;
    
    PP[0]=0;PP[1]=b;PP[2]=a;
    Fc.P1 = PP;
    PP[0]=-b;PP[1]=a;PP[2]=0;
    Fc.P2 = PP;
    PP[0]=b;PP[1]=a;PP[2]=0;
    Fc.P3 = PP;
    FacetArray[1] = Fc;
    
    PP[0]=0;PP[1]=b;PP[2]=a;
    Fc.P1 = PP;
    PP[0]=0;PP[1]=-b;PP[2]=a;
    Fc.P2 = PP;
    PP[0]=-a;PP[1]=0;PP[2]=b;
    Fc.P3 = PP;
    FacetArray[2] = Fc;
    
    PP[0]=0;PP[1]=b;PP[2]=a;
    Fc.P1 = PP;
    PP[0]=a;PP[1]=0;PP[2]=b;
    Fc.P2 = PP;
    PP[0]=0;PP[1]=-b;PP[2]=a;
    Fc.P3 = PP;
    FacetArray[3] = Fc;
    
    PP[0]=0;PP[1]=b;PP[2]=-a;
    Fc.P1 = PP;
    PP[0]=0;PP[1]=-b;PP[2]=-a;
    Fc.P2 = PP;
    PP[0]=a;PP[1]=0;PP[2]=-b;
    Fc.P3 = PP;
    FacetArray[4] = Fc;
    
    PP[0]=0;PP[1]=b;PP[2]=-a;
    Fc.P1 = PP;
    PP[0]=-a;PP[1]=0;PP[2]=-b;
    Fc.P2 = PP;
    PP[0]=0;PP[1]=-b;PP[2]=-a;
    Fc.P3 = PP;
    FacetArray[5] = Fc;
    
    PP[0]=0;PP[1]=-b;PP[2]=a;
    Fc.P1 = PP;
    PP[0]=b;PP[1]=-a;PP[2]=0;
    Fc.P2 = PP;
    PP[0]=-b;PP[1]=-a;PP[2]=0;
    Fc.P3 = PP;
    FacetArray[6] = Fc;
    
    PP[0]=0;PP[1]=-b;PP[2]=-a;
    Fc.P1 = PP;
    PP[0]=-b;PP[1]=-a;PP[2]=0;
    Fc.P2 = PP;
    PP[0]=b;PP[1]=-a;PP[2]=0;
    Fc.P3 = PP;
    FacetArray[7] = Fc;
    
    PP[0]=-b;PP[1]=a;PP[2]=0;
    Fc.P1 = PP;
    PP[0]=-a;PP[1]=0;PP[2]=b;
    Fc.P2 = PP;
    PP[0]=-a;PP[1]=0;PP[2]=-b;
    Fc.P3 = PP;
    FacetArray[8] = Fc;
    
    PP[0]=-b;PP[1]=-a;PP[2]=0;
    Fc.P1 = PP;
    PP[0]=-a;PP[1]=0;PP[2]=-b;
    Fc.P2 = PP;
    PP[0]=-a;PP[1]=0;PP[2]=b;
    Fc.P3 = PP;
    FacetArray[9] = Fc;
     
    PP[0]=b;PP[1]=a;PP[2]=0;
    Fc.P1 = PP;
    PP[0]=a;PP[1]=0;PP[2]=-b;
    Fc.P2 = PP;
    PP[0]=a;PP[1]=0;PP[2]=b;
    Fc.P3 = PP;
    FacetArray[10] = Fc;
     
    PP[0]=b;PP[1]=-a;PP[2]=0;
    Fc.P1 = PP;
    PP[0]=a;PP[1]=0;PP[2]=b;
    Fc.P2 = PP;
    PP[0]=a;PP[1]=0;PP[2]=-b;
    Fc.P3 = PP;
    FacetArray[11] = Fc;
     
    PP[0]=0;PP[1]=b;PP[2]=a;
    Fc.P1 = PP;
    PP[0]=-a;PP[1]=0;PP[2]=b;
    Fc.P2 = PP;
    PP[0]=-b;PP[1]=a;PP[2]=0;
    Fc.P3 = PP;
    FacetArray[12] = Fc;
     
    PP[0]=0;PP[1]=b;PP[2]=a;
    Fc.P1 = PP;
    PP[0]=b;PP[1]=a;PP[2]=0;
    Fc.P2 = PP;
    PP[0]=a;PP[1]=0;PP[2]=b;
    Fc.P3 = PP;
    FacetArray[13] = Fc;
     
    PP[0]=0;PP[1]=b;PP[2]=-a;
    Fc.P1 = PP;
    PP[0]=-b;PP[1]=a;PP[2]=0;
    Fc.P2 = PP;
    PP[0]=-a;PP[1]=0;PP[2]=-b;
    Fc.P3 = PP;
    FacetArray[14] = Fc;
     
    PP[0]=0;PP[1]=b;PP[2]=-a;
    Fc.P1 = PP;
    PP[0]=a;PP[1]=0;PP[2]=-b;
    Fc.P2 = PP;
    PP[0]=b;PP[1]=a;PP[2]=0;
    Fc.P3 = PP;
    FacetArray[15] = Fc;
     
    PP[0]=0;PP[1]=-b;PP[2]=-a;
    Fc.P1 = PP;
    PP[0]=-a;PP[1]=0;PP[2]=-b;
    Fc.P2 = PP;
    PP[0]=-b;PP[1]=-a;PP[2]=0;
    Fc.P3 = PP;
    FacetArray[16] = Fc;
      
    PP[0]=0;PP[1]=-b;PP[2]=-a;
    Fc.P1 = PP;
    PP[0]=b;PP[1]=-a;PP[2]=0;
    Fc.P2 = PP;
    PP[0]=a;PP[1]=0;PP[2]=-b;
    Fc.P3 = PP;
    FacetArray[17] = Fc;
       
    PP[0]=0;PP[1]=-b;PP[2]=a;
    Fc.P1 = PP;
    PP[0]=-b;PP[1]=-a;PP[2]=0;
    Fc.P2 = PP;
    PP[0]=-a;PP[1]=0;PP[2]=b;
    Fc.P3 = PP;
    FacetArray[18] = Fc;
       
    PP[0]=0;PP[1]=-b;PP[2]=a;
    Fc.P1 = PP;
    PP[0]=a;PP[1]=0;PP[2]=b;
    Fc.P2 = PP;
    PP[0]=b;PP[1]=-a;PP[2]=0;
    Fc.P3 = PP;
    FacetArray[19] = Fc;
 
    for (unsigned j = 0; j < facets; j++)
      {
      // Find a line perpendicular to each face
      LType L, A, B;
      A = FacetArray[j].P2 - FacetArray[j].P1;
      B = FacetArray[j].P3 - FacetArray[j].P1;
      L[0] = A[1]*B[2] - B[1]*A[2];
      L[1] = B[0]*A[2] - A[0]*B[2];
      L[2] = A[0]*B[1] - B[0]*A[1];
      
      L.Normalize();
      // Scale to required length
      L*= rr;
      if (!res.checkParallel(L, res.m_Lines))
	{
	res.m_Lines.push_back(L);
	}
      }
    return(res);
    }
    break;
    case 32:
    {
    iterations = 1;
    // 2 iterations leads to 128 faces, which is too many
    // subdivision of octahedron
    // create triangular facet approximation to a sphere - begin with
    // unit sphere
    // total number of facets is 8 * (4^iterations)
    int facets = 8 * (int)pow(4, iterations);
    float sqrt2 = sqrt(2.0);
    std::cout << facets << " facets" << std::endl;
    typedef std::vector<FacetType> FacetArrayType;
    FacetArrayType FacetArray = FacetArrayType(facets);
    
    // original corners of octahedron
    LType P0, P1, P2, P3, P4, P5; 
    P0[0]=0;         P0[1]=0;       P0[2]=1;
    P1[0]=0;         P1[1]=0;       P1[2]=-1;
    P2[0]=-1.0/sqrt2;P2[1]=-1/sqrt2;P2[2]=0;
    P3[0]=1/sqrt2;   P3[1]=-1/sqrt2;P3[2]=0;
    P4[0]=1/sqrt2;   P4[1]=1/sqrt2; P4[2]=0;
    P5[0]=-1/sqrt2;  P5[1]=1/sqrt2; P5[2]=0;
    
    FacetType F0, F1, F2, F3, F4, F5, F6, F7;
    F0.P1 = P0; F0.P2 = P3; F0.P3 = P4;
    F1.P1 = P0; F1.P2 = P4; F1.P3 = P5;
    F2.P1 = P0; F2.P2 = P5; F2.P3 = P2;
    F3.P1 = P0; F3.P2 = P2; F3.P3 = P3;
    F4.P1 = P1; F4.P2 = P4; F4.P3 = P3;
    F5.P1 = P1; F5.P2 = P5; F5.P3 = P4;
    F6.P1 = P1; F6.P2 = P2; F6.P3 = P5;
    F7.P1 = P1; F7.P2 = P3; F7.P3 = P2;
    
    FacetArray[0] = F0;
    FacetArray[1] = F1;
    FacetArray[2] = F2;
    FacetArray[3] = F3;
    FacetArray[4] = F4;
    FacetArray[5] = F5;
    FacetArray[6] = F6;
    FacetArray[7] = F7;
    int pos = 8;
    // now subdivide the octahedron
    for (unsigned it = 0; it<iterations; it++)
      {
      // Bisect edges and move to sphere
      unsigned ntold = pos;
      for (unsigned i = 0; i < ntold; i++)
	{
	LType Pa, Pb, Pc;
	for (unsigned d = 0; d<  VDimension;d++)
	  {
	  Pa[d] = (FacetArray[i].P1[d] + FacetArray[i].P2[d])/2;
	  Pb[d] = (FacetArray[i].P2[d] + FacetArray[i].P3[d])/2;
	  Pc[d] = (FacetArray[i].P3[d] + FacetArray[i].P1[d])/2;
	  }
	Pa.Normalize();
	Pb.Normalize();
	Pc.Normalize();
	FacetArray[pos].P1 = FacetArray[i].P1; 
	FacetArray[pos].P2 = Pa; 
	FacetArray[pos].P3 = Pc;
	pos++;
	FacetArray[pos].P1 = Pa;
	FacetArray[pos].P2 = FacetArray[i].P2;  
	FacetArray[pos].P3 = Pb;
	pos++;
	FacetArray[pos].P1 = Pb;
	FacetArray[pos].P2 = FacetArray[i].P3;  
	FacetArray[pos].P3 = Pc;
	pos++;
	FacetArray[i].P1 = Pa;
	FacetArray[i].P2 = Pb;
	FacetArray[i].P3 = Pc;
	}
      }
    
      
    for (unsigned j = 0; j < facets; j++)
      {
      // Find a line perpendicular to each face
      LType L, A, B;
      A = FacetArray[j].P2 - FacetArray[j].P1;
      B = FacetArray[j].P3 - FacetArray[j].P1;
      L[0] = A[1]*B[2] - B[1]*A[2];
      L[1] = B[0]*A[2] - A[0]*B[2];
      L[2] = A[0]*B[1] - B[0]*A[1];
      
      L.Normalize();
      // Scale to required length
      L*= rr;
      if (!res.checkParallel(L, res.m_Lines))
	{
	res.m_Lines.push_back(L);
	}
      }
    return(res);
    }
    break;
    default:
      std::cout << "Unsupported number of lines" << std::endl;
      return(res);
    }
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
  res.m_Decomposable = true;
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
    if ((1.0 - fabs(L)) < 0.000001) return(true);
    }
  return(false);
}

template<unsigned int VDimension>
void FlatStructuringElement< VDimension >
::PrintSelf(std::ostream &os, Indent indent) const
{
  //Superclass::PrintSelf(os, indent);
  if (m_Decomposable)
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
