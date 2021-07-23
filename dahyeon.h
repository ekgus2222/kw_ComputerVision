#ifndef DAHYEON_H
#define DAHYEON_H

#include "kfc.h"
#include <iostream>
#include <vector>



void Canny_Edge(KImageGray& imgSrc,double sigma);
KImageDouble Gaussian_Filter(KImageGray& icMain, double sigma);
KImageGray Show(std::vector<std::vector<KImageDouble>>& Octave);
std::vector<std::vector<real_key>> Find_Key_Point(std::vector<std::vector<KImageDouble>>& Dog);
bool IsPeak(std::vector<KImageDouble>& Dog_Scale,int i,int ii, int jj);
void Orientation(std::vector<std::vector<real_key>>& KeyPoint, std::vector<KImageDouble> Octave_0,double sigma);
void Mark_KeyPoint(KImageGray &igScale, double mag, double ori_deg, int pX, int pY);
void KeyPoint_descriptor(std::vector<std::vector<real_key>>& KeyPoint, std::vector<KImageDouble> Octave_0,double sigma);

void Optical_Flow(KImageGray Image0, KImageGray Image1, UV** mat_uv);
void Draw(KImageGray& Img,UV ** mat_uv);
#endif // DAHYEON_H
