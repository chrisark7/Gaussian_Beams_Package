(* ::Package:: *)

BeginPackage["GaussianBeams`"]


Beam\[Omega]::usage="Beam\[Omega][z]; Gives the beam size as a function of z.  The Beam_ functions use the waist size, \[Omega]0, the Rayleigh range, zr, and the wavelength, \[Lambda], as if they are already defined or will be replaced later.";
BeamR::usage="BeamR[z]; Gives the radius of curvature as a function of z.  The Beam_ functions use the waist size, \[Omega]0, the Rayleigh range, zr, and the wavelength, \[Lambda], as if they are already defined or will be replaced later.";
Beam\[Eta]::usage="Beam\[Eta][z]; Gives the Gouy phase as a function of z relative to z=0.  The Beam_ functions use the waist size, \[Omega]0, the Rayleigh range, zr, and the wavelength, \[Lambda], as if they are already defined or will be replaced later.";
BeamU::usage="BeamU[n,x,z]; Gives the 1D cross section of a laser field of mode order n.  The Beam_ functions use the waist size, \[Omega]0, the Rayleigh range, zr, and the wavelength, \[Lambda], as if they are already defined or will be replaced later.";
BeamHG::usage="BeamHG[m,n,x,y,z]; Gives the full 2D laser field with unit power of mode order m,n.  The Beam_ functions use the waist size, \[Omega]0, the Rayleigh range, zr, and the wavelength, \[Lambda], as if they are already defined or will be replaced later.";
BeamHGProp::usage="BeamHGProp[m,n,z2,z1]; The laser field propogator from location z1 to location z2.  Includes the propogation phase and the Gouy phase.  The Beam_ functions use the waist size, \[Omega]0, the Rayleigh range, zr, and the wavelength, \[Lambda], as if they are already defined or will be replaced later.";
BeamHGNum::usage="BeamHGNum[m,n]; A renumbering of the HG modes ({m,n}->j) for use with scattering matrix calculations."
BeamHGNumInv::usage="BeamHGNumInv[j]; The inverse of the numbering defined by BeamHGNum (j->{m,n})"


RayT::usage="RayT[d]; Gives the ray matrix for translation by a distance d.  It is assumed that the distance is already scaled by the index of refraction.  For propogation of a q parameter in a medium remember that it is the reduced distance, \!\(\*FractionBox[\(d\), \(n\)]\), which matters.";
RayM::usage="RayM[r]; Gives the ray matrix for reflection off of a curved mirror of radius of curvature r.";
RayL::usage="RayL[f]; Gives the ray matrix for a beam passing through a thin lens of focal length f.";
RayIF::usage="RayIF[ni,nf]; Gives the ray matrix for propogation through a flat interface between two media.  ni is the index of refraction of the initial medium and nf is the index of refraction of the final medium.";
RayIC::usage="RayIC[ni,nf,r]; Gives the ray matrix for propogation through a curved interface between two media.  ni is the index of refraction of the initial medium, nf is the index of refraction of the final medium, and r is the radius of curvature of the interface.";


Q\[Omega]::usage="Q\[Omega][q]; Gives the beam size of a Gaussian beam described by the q parameter.";
Q\[Omega]M::usage="Q\[Omega]M[q,m]; Gives the beam size of a Gaussian beam described by the q parameter with an \!\(\*SuperscriptBox[\(M\), \(2\)]\) value of \!\(\*SuperscriptBox[\(m\), \(2\)]\).";
QR::usage="QR[q]; Gives the radius of curvature of a Gaussian beam described by the q parameter.";
Q\[Omega]0::usage="Q\[Omega]0[q]; Gives the waist size of a Guassian beam described by the q paremeter.";
Q\[Eta]::usage="Q\[Eta][q]; Gives the Gouy phase relative to the waist of the Gaussian beam described the q paremeter.";
Q\[Theta]::usage="Q\[Theta][q]; Gives the divergence angle of the Gaussian beam described by the q parameter.";
QOvlpField::usage="QOvlpField[q1,q2]; Gives the complex value of the overlap integral between two Gaussian beams of unit power described by beam parameters q1 and q2.";
QOvlpPow::usage="QOvlpPow[q1,q2]; Gives the real value of the square of the overlap integral between two Gaussian beams fo unit power described by the beam parameters q1 and q2.";
QProp::usage="QProp[q,rayMat]; Propogates the q parameter through the optical system described by the ABCD matrix rayMat.";
QCav::usage="QCav[rayMat]; Calculates the stable q parameter in a cavity described by the round trip ray matrix rayMat.";
QReverse::usage="QReverse[q]; Reverses the direction of the q parameter by switching the sign of the real part.";


Ray4M::usage="Ray4M[r,\[Epsilon],\[Epsilon]p];  The 4x4 misaligned ray matrix for a curved mirror of racius r and misalignments \[Epsilon] and \[Epsilon]p (position and angle).  Follows the prescription given by Wang Shaomin in 'Matrix Methods in Treating Decenterd Optical Systems.' Opt. Quant. Electron. 17;1;1-14.";
Ray4T::usage="Ray4T[d]; The 4.4 misaligned ray matrix for translation through free space over a distance d.  .  Follows the prescription given by Wang Shaomin in 'Matrix Methods in Treating Decenterd Optical Systems.' Opt. Quant. Electron. 17;1;1-14.";
Ray4IC::usage="Ray4IC[ni,nf,r,\[Epsilon],\[Epsilon]p]; The 4x4 misaligned ray matrix for translation through a curved interface of radius r, initial index of refraction ni, final index of refraction nf, and misalignments \[Epsilon] and \[Epsilon]p (position and angle).  Follows the prescription given by Wang Shaomin in 'Matrix Methods in Treating Decenterd Optical Systems.' Opt. Quant. Electron. 17;1;1-14.";
Ray4IF::usage="Ray4IF[ni,nf,\[Epsilon],\[Epsilon]p]; The 4x4 misaligned ray matrix for translation through a flat interface with initial index of refraction ni, final index of refraction nf, and misalignments \[Epsilon] and \[Epsilon]p (position and angle).  Follows the prescription given by Wang Shaomin in 'Matrix Methods in Treating Decenterd Optical Systems.' Opt. Quant. Electron. 17;1;1-14.";
Ray4L::usage="Ray4L[f,\[Epsilon],\[Epsilon]p];  The 4x4 misaligned ray matrix for translation through a thin lens of focal length f, and misalignments \[Epsilon] and \[Epsilon]p (position and angle).  Follows the prescription given by Wang Shaomin in 'Matrix Methods in Treating Decenterd Optical Systems.' Opt. Quant. Electron. 17;1;1-14.";


Begin["`Private`"]


Beam\[Omega][z_]:=Global`\[Omega]0 Sqrt[1+(z/Global`zr)^2];
BeamR[z_]:=z(1+(Global`zr/z)^2);
Beam\[Eta][z_]:=ArcTan[z/Global`zr];
BeamU[n_,x_,z_]:=(2/Pi)^(1/4) (1/(2^n Factorial[n]Beam\[Omega][z]))^(1/2) HermiteH[n,Sqrt[2]x/Beam\[Omega][z]]Exp[-x^2(1/Beam\[Omega][z]^2+I(Pi/(Global`\[Lambda] BeamR[z])))];
BeamHG[m_,n_,x_,y_,z_]:=BeamU[m,x,z]BeamU[n,y,z];
BeamHGProp[m_,n_,z2_,z1_]:=Exp[-I 2Pi/Global`\[Lambda](z2-z1)]Exp[I(m+n+1)(Beam\[Eta][z2]-Beam\[Eta][z1])];
BeamHGNum[m_,n_]:=(n+m+1)/2(n+m)+m;
BeamHGNumInv[j_]:={m,n}/.Solve[{(n+m+1)/2(n+m)+m==j,m>=0,n>=0},{m,n},Integers][[1]];


RayT[d_]:={{1,d},{0,1}};
RayM[r_]:={{1,0},{-2/r,1}};
RayL[f_]:={{1,0},{-1/f,1}};
RayIF[ni_,nf_]:={{1,0},{0,ni/nf}}; 
RayIC[ni_,nf_,r_]:={{1,0},{(ni-nf)/(r*nf),ni/nf}}; 


Q\[Omega][q_]:=Sqrt[Global`\[Lambda]/Pi Im[q](1+(Re[q]/Im[q])^2)];
Q\[Omega]M[q_,m_]:=Sqrt[Global`\[Lambda]/Pi Im[q](1+m^4(Re[q]/Im[q])^2)];
QR[q_]:=Im[q](Re[q]/Im[q]+Im[q]/Re[q]);
Q\[Omega]0[q_]:=Sqrt[Global`\[Lambda]/Pi Im[q]];
Q\[Eta][q_]:=ArcTan[Re[q]/Im[q]];
Q\[Theta][q_]:=Global`\[Lambda]/(Q\[Omega]0[q] Pi);
QOvlpField[q1_,q2_]:=(4 Im[q1]Im[q2])/Abs[Conjugate[q1]-q2]^2;
QOvlpPow[q1_,q2_]:=Abs[QOvlpField[q1,q2]]^2;
QProp[q_,{{a_,b_},{c_,d_}}]:=(a q+b)/(c q+d);
QCav[{{a_,b_},{c_,d_}}]:=If[c<0,1/(2c)(a-d-Sqrt[(a-d)^2+4b c]),1/(2c)(a-d+Sqrt[(a-d)^2+4b c])]
QReverse[q_]:=q-2*Re[q];


Ray4M[r_,\[Epsilon]_,\[Epsilon]p_]={{1,0,0,0},{-2/r,1,-2/r \[Epsilon],-2\[Epsilon]p},{0,0,1,0},{0,0,0,1}};
Ray4T[d_]={{1,d,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
Ray4IC[ni_,nf_,r_,\[Epsilon]_,\[Epsilon]p_]={{1,0,0,0},{-(1-ni/nf)/r,ni/nf,(1-ni/nf)/r \[Epsilon],(1-ni/nf)\[Epsilon]p},{0,0,1,0},{0,0,0,1}};
Ray4IF[ni_,nf_,\[Epsilon]_,\[Epsilon]p_]={{1,0,0,0},{0,ni/nf,0,(1-ni/nf)\[Epsilon]p},{0,0,1,0},{0,0,0,1}};
Ray4L[f_,\[Epsilon]_,\[Epsilon]p_]={{1,0,0,0},{-1/f,1,1/f \[Epsilon],0},{0,0,1,0},{0,0,0,1}};


End[]
EndPackage[]
