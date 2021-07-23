#include "dahyeon.h"
#include <iostream>
#include <QFileDialog>
#include <QPainter>
#include <math.h>
#include <QDebug>
#include <algorithm>
#include <stdlib.h>
//#include <Windows.h>
#include <algorithm>
#include <vector>
#include <stack>
#include <fstream>
#include "eigen-3.3.8/Eigen/Dense"


#include "mainframe.h"
#include "ui_mainframe.h"
#include "imageform.h"

KImageDouble Gaussian_Filter(KImageGray& icMain, double sigma)
{

    KImageDouble Return(icMain);

    int row = icMain.Row();
    int col = icMain.Col();

    int Gau_filter_size = 8*sigma + 1; //filter size

    int size = std::sqrt(Gau_filter_size);

    //qDebug()<<size;

    double mask;



    for(int ii = 0; ii<row;ii++)
    {
        for(int jj = 0; jj<col;jj++)
        {

            if(ii<size || ii>=row-size || jj<size || jj>=col-size)
            {
                Return[ii][jj] = 0;
                continue;
            }
            mask = 0;


            for(int i = -size; i<size+1; i++)
            {
                for(int j = -size; j<size+1;j++)
                {
                    mask += std::exp(-0.5*((i*i+j*j)/(sigma*sigma)))*icMain[ii-i][jj-j];

                }
            }

            Return[ii][jj] = (1/(2*M_PI*sigma*sigma))*mask;
        }

    }

    //qDebug()<<"good";
    return Return;
}


KImageGray Show(std::vector<std::vector<KImageDouble>>& Octave){

    int octave_que_size = Octave.size();
    int merged_row = Octave.front().front().Row() * 3;
    int merged_col = Octave.front().front().Col() * Octave.front().size();

    if (merged_col >= 1920) {
        merged_col = 1920;
    }
    KImageGray Base(merged_row, merged_col);


    std::vector<KImageDouble> img_vec;
    KImageGray now;

    int sum_of_prev_row = 0;


    for(int octave = 0; octave<Octave.size();octave++)
    {

        img_vec = Octave[octave];

        int img_vec_size = img_vec.size(); //6
        int each_row = img_vec.front().Row(); //359
        int each_col = img_vec.front().Col(); //640
        int poor_column = 0;


        for(int img = 0; img<img_vec_size;img++)
        {

            //qDebug()<<"ok";
            now = img_vec[img].ToGray();

            for(int i=0;i<each_row;i++)
            {
                for(int j =0; j<each_col;j++)
                {
                    Base[i+sum_of_prev_row][j+(img-poor_column)*each_col] = now[i][j];
                }
            }

            // 다음에 이어 붙일 이미지가 존재하지만 모니터의 가로 길이가 부족할 때 한 칸 아래로 내림
            if (each_col - 1 + ((img - poor_column) + 1) * each_col > merged_col && img + 1 < img_vec_size)
            {
                  sum_of_prev_row += each_row;
                  poor_column = img + 1;
            }
            //qDebug()<<"end";

        }
        sum_of_prev_row += each_row;

    }


    return Base;
}

bool IsPeak(std::vector<KImageDouble>& Dog_Scale,int i,int ii, int jj)
{
    if(Dog_Scale[i][ii][jj]<0.03) return false;

    if(Dog_Scale[i][ii][jj]>0) //극대
    {
         if(Dog_Scale[i][ii][jj]<Dog_Scale[i-1][ii-1][jj-1]||Dog_Scale[i][ii][jj]<Dog_Scale[i-1][ii-1][jj]||Dog_Scale[i][ii][jj]<Dog_Scale[i-1][ii-1][jj+1]||
            Dog_Scale[i][ii][jj]<Dog_Scale[i-1][ii][jj-1]||Dog_Scale[i][ii][jj]<Dog_Scale[i-1][ii][jj]||Dog_Scale[i][ii][jj]<Dog_Scale[i-1][ii][jj+1]||
            Dog_Scale[i][ii][jj]<Dog_Scale[i-1][ii+1][jj-1]||Dog_Scale[i][ii][jj]<Dog_Scale[i-1][ii+1][jj]||Dog_Scale[i][ii][jj]<Dog_Scale[i-1][ii+1][jj+1]||


            Dog_Scale[i][ii][jj]<Dog_Scale[i][ii-1][jj-1]||Dog_Scale[i][ii][jj]<Dog_Scale[i][ii-1][jj]||Dog_Scale[i][ii][jj]<Dog_Scale[i][ii-1][jj+1]||
            Dog_Scale[i][ii][jj]<Dog_Scale[i][ii][jj-1]||                                               Dog_Scale[i][ii][jj]<Dog_Scale[i][ii][jj+1]||
            Dog_Scale[i][ii][jj]<Dog_Scale[i][ii+1][jj-1]||Dog_Scale[i][ii][jj]<Dog_Scale[i][ii+1][jj]||Dog_Scale[i][ii][jj]<Dog_Scale[i][ii+1][jj+1]||

            Dog_Scale[i][ii][jj]<Dog_Scale[i+1][ii-1][jj-1]||Dog_Scale[i][ii][jj]<Dog_Scale[i+1][ii-1][jj]||Dog_Scale[i][ii][jj]<Dog_Scale[i+1][ii-1][jj+1]||
            Dog_Scale[i][ii][jj]<Dog_Scale[i+1][ii][jj-1]||Dog_Scale[i][ii][jj]<Dog_Scale[i+1][ii][jj]||Dog_Scale[i][ii][jj]<Dog_Scale[i][ii+1][jj+1]||
            Dog_Scale[i][ii][jj]<Dog_Scale[i+1][ii+1][jj-1]||Dog_Scale[i][ii][jj]<Dog_Scale[i+1][ii+1][jj]||Dog_Scale[i][ii][jj]<Dog_Scale[i+1][ii+1][jj+1]) return false;

     }
     else
     {
         if(Dog_Scale[i][ii][jj]>Dog_Scale[i-1][ii-1][jj-1]||Dog_Scale[i][ii][jj]>Dog_Scale[i-1][ii-1][jj]||Dog_Scale[i][ii][jj]>Dog_Scale[i-1][ii-1][jj+1]||
            Dog_Scale[i][ii][jj]>Dog_Scale[i-1][ii][jj-1]||Dog_Scale[i][ii][jj]>Dog_Scale[i-1][ii][jj]||Dog_Scale[i][ii][jj]>Dog_Scale[i-1][ii][jj+1]||
            Dog_Scale[i][ii][jj]>Dog_Scale[i-1][ii+1][jj-1]||Dog_Scale[i][ii][jj]>Dog_Scale[i-1][ii+1][jj]||Dog_Scale[i][ii][jj]>Dog_Scale[i-1][ii+1][jj+1]||


            Dog_Scale[i][ii][jj]>Dog_Scale[i][ii-1][jj-1]||Dog_Scale[i][ii][jj]>Dog_Scale[i][ii-1][jj]||Dog_Scale[i][ii][jj]>Dog_Scale[i][ii-1][jj+1]||
            Dog_Scale[i][ii][jj]>Dog_Scale[i][ii][jj-1]||                                               Dog_Scale[i][ii][jj]>Dog_Scale[i][ii][jj+1]||
            Dog_Scale[i][ii][jj]>Dog_Scale[i][ii+1][jj-1]||Dog_Scale[i][ii][jj]>Dog_Scale[i][ii+1][jj]||Dog_Scale[i][ii][jj]>Dog_Scale[i][ii+1][jj+1]||

            Dog_Scale[i][ii][jj]>Dog_Scale[i+1][ii-1][jj-1]||Dog_Scale[i][ii][jj]>Dog_Scale[i+1][ii-1][jj]||Dog_Scale[i][ii][jj]>Dog_Scale[i+1][ii-1][jj+1]||
            Dog_Scale[i][ii][jj]>Dog_Scale[i+1][ii][jj-1]||Dog_Scale[i][ii][jj]>Dog_Scale[i+1][ii][jj]||Dog_Scale[i][ii][jj]>Dog_Scale[i][ii+1][jj+1]||
            Dog_Scale[i][ii][jj]>Dog_Scale[i+1][ii+1][jj-1]||Dog_Scale[i][ii][jj]>Dog_Scale[i+1][ii+1][jj]||Dog_Scale[i][ii][jj]>Dog_Scale[i+1][ii+1][jj+1]) return false;
      }

    return true;

}

std::vector<std::vector<real_key>> Find_Key_Point(std::vector<std::vector<KImageDouble>>& Dog)
{

    double thres = 1.5;
    double r = thres;
    double b = pow(r+1,2)/r;
    double Dxx,Dyy,Dxy;
    double TrH, DetH;
    double a;

    real_key key;
    std::vector<std::vector<real_key>> key_vec(2);

    //scale = 0;
    for(int dog = 0; dog<3; dog++){ //octave 별 dog 선택
        for(int i=1;i<4;i++) //scale 변화
        {
            for(int ii=1;ii<Dog[dog][i].Row()-1;ii++)
            {
                for(int jj=1;jj<Dog[dog][i].Col()-1;jj++)
                {
                    //qDebug() << dog <<" "<<ii<<jj;


                    if(IsPeak(Dog[dog],i,ii,jj)) //6개 이미지 넘겨줌
                    {


                        Dxx = Dog[dog][i][ii][jj+1]+Dog[dog][i][ii][jj-1]-2*Dog[dog][i][ii][jj];
                        Dyy = Dog[dog][i][ii+1][jj]+Dog[dog][i][ii-1][jj]-2*Dog[dog][i][ii][jj];
                        Dxy = ((Dog[dog][i][ii+1][jj+1]-Dog[dog][i][ii+1][jj-1]) - (Dog[dog][i][ii-1][jj+1] - Dog[dog][i][ii-1][jj-1]))/4.0;
                        TrH = Dxx + Dyy;
                        DetH = Dxx*Dyy - Dxy*Dxy;

                        a = TrH*TrH/DetH;

                        //qDebug() << Dxx<<" "<<Dyy<<" "<<" " <<Dxy<<" " << DetH;
                        if(DetH<=0) continue;
                        //if(a < b) KeyPoint[ii][jj] = 255;
                        if(a >= b) continue;

                        //KeyPoint
                        key.x = ii*pow(2,dog);
                        key.y = jj*pow(2,dog);
                        key.scale = i;


                        key_vec[i-1].push_back(key);

                    }

                }

            }

        }
    }

    return key_vec;
}

//double get_mag(KImageGray& icMain, )

void Orientation(std::vector<std::vector<real_key>>& KeyPoint, std::vector<KImageDouble> Octave_0,double sigma)
{

    KImageDouble Scale_img;
    int x,y;
    double mag;
    double phase;
    double gau_kenel;

    double bucket[36] = {0,};

    int x_t,y_t;

    int index;

    double max = 0.;
    int dir;

    //int cnt = 1;

    int KeyPoint_scale_size = 0;

    for(int scale = 0; scale<KeyPoint.size(); scale++) //root 2 sigma, 2 sigma
    {
        Scale_img = Octave_0[scale+1]; //origin roo2, origin 2
        sigma = sigma*pow(1.4,scale+1);
        KeyPoint_scale_size = KeyPoint[scale].size();
        for(int i=0;i<KeyPoint_scale_size;i++)
        {
            x = KeyPoint[scale][i].x;
            y = KeyPoint[scale][i].y;

            //16개 좌표 x+w_x, y+w_y
            for(int w_x = -5; w_x < 6; w_x++)
            {
                for(int w_y = -5; w_y < 6; w_y++)
                {
                    x_t = x+w_x;
                    y_t = y+w_y;

                    if(x_t>Scale_img.Row()-2||x_t<1||y_t>Scale_img.Col()-2||y_t<1) continue;
                    gau_kenel = 1./(2*M_PI*sigma*sigma) * exp(-0.5*(w_x*w_x + w_y*w_y)/(sigma*sigma));

                    mag = gau_kenel * sqrt(pow(Scale_img[x_t+1][y_t]-Scale_img[x_t-1][y_t],2)
                                       +pow(Scale_img[x_t][y_t+1]-Scale_img[x_t][y_t-1],2));

                    phase = atan2(Scale_img[x_t][y_t+1]-Scale_img[x_t][y_t-1],
                                  Scale_img[x_t+1][y_t]-Scale_img[x_t-1][y_t])*180/M_PI+180.0;

                    //if(phase>360) phase -= 360;

                    //if(phase<0) phase += 180;

                    index = phase/10;

                    bucket[index] += mag;

                    //qDebug() << scale << " " << i << " " << x << " " << y << " " << x_t << " " << y_t << " " << index << mag;
                    //qDebug() << index;
                }
            }

            for(int b = 0; b<36; b++)
            {
                if(bucket[b]>max)
                {
                    max = bucket[b];
                    dir = b*10;
                }

                qDebug()<< b << " " << bucket[b];
            }

            //qDebug()<< dir << " " << max;

            for(int b = 0; b<36; b++)
            {
                if(bucket[b]>max*0.8&& b != dir)
                {
                    KeyPoint[scale].push_back({KeyPoint[scale][i].x,KeyPoint[scale][i].y,(double)b*10,bucket[b]});
                   //qDebug()<< b << " " << bucket[b];
                }

                //qDebug()<< b << " " << bucket[b];
            }


            KeyPoint[scale][i].direction = dir;
            KeyPoint[scale][i].magnitude = max;
            //qDebug()<< cnt << " " << dir << bucket[8] << bucket[9];
            //cnt++;


        }

    }


}


void KeyPoint_descriptor(std::vector<std::vector<real_key>>& KeyPoint, std::vector<KImageDouble> Octave_0,double sigma)
{

    KImageDouble Scale_img;
    int x,y;
    double mag[81];
    double phase[81];
    double gau_kenel;

    double bucket[8] = {0,};

    int x_t,y_t;




    //int cnt = 1;

    int KeyPoint_scale_size = 0;

    for(int scale = 0; scale<KeyPoint.size(); scale++) //root 2 sigma, 2 sigma
    {
        Scale_img = Octave_0[scale+1]; //origin roo2, origin 2
        sigma = sigma*pow(1.4,scale+1);
        KeyPoint_scale_size = KeyPoint[scale].size();
        for(int i=0;i<KeyPoint_scale_size;i++)
        {


            x = KeyPoint[scale][i].x;
            y = KeyPoint[scale][i].y;

            //11개 좌표 x+w_x, y+w_y
            for(int w_x = -4; w_x < 5; w_x++)
            {
                for(int w_y = -4; w_y < 5; w_y++)
                {
                    x_t = x+w_x;
                    y_t = y+w_y;

                    if(x_t>Scale_img.Row()-2||x_t<1||y_t>Scale_img.Col()-2||y_t<1) continue;
                    gau_kenel = 1./(2*M_PI*sigma*sigma) * exp(-0.5*(w_x*w_x + w_y*w_y)/(sigma*sigma));

                    mag[(w_x+4)*8+w_y+4] = gau_kenel * sqrt(pow(Scale_img[x_t+1][y_t]-Scale_img[x_t-1][y_t],2)
                                                            +pow(Scale_img[x_t][y_t+1]-Scale_img[x_t][y_t-1],2));

                    phase[(w_x+4)*8+w_y+4] = atan2(Scale_img[x_t][y_t+1]-Scale_img[x_t][y_t-1],
                                                    Scale_img[x_t+1][y_t]-Scale_img[x_t-1][y_t])*180/M_PI+180.0;

                }
            }

            //index = phase/45
            for(int p = 0; p < 2; p++)
            {
                std::vector<double> dou_vec;

                bucket[(int)phase[0+5*p]/45] += mag[0+5*p];
                bucket[(int)phase[1+5*p]/45] += mag[1+5*p];
                bucket[(int)phase[2+5*p]/45] += mag[2+5*p];
                bucket[(int)phase[3+5*p]/45] += mag[3+5*p];

                bucket[(int)phase[9+5*p]/45] += mag[9+5*p];
                bucket[(int)phase[10+5*p]/45] += mag[10+5*p];
                bucket[(int)phase[11+5*p]/45] += mag[11+5*p];
                bucket[(int)phase[12+5*p]/45] += mag[12+5*p];

                bucket[(int)phase[18+5*p]/45] += mag[18+5*p];
                bucket[(int)phase[19+5*p]/45] += mag[19+5*p];
                bucket[(int)phase[20+5*p]/45] += mag[20+5*p];
                bucket[(int)phase[21+5*p]/45] += mag[21+5*p];

                bucket[(int)phase[27+5*p]/45] += mag[27+5*p];
                bucket[(int)phase[28+5*p]/45] += mag[28+5*p];
                bucket[(int)phase[29+5*p]/45] += mag[29+5*p];
                bucket[(int)phase[30+5*p]/45] += mag[30+5*p];


                for(int b = 0;b<8;b++)
                {
                    dou_vec.push_back(bucket[b]);
                }

                KeyPoint[scale][i].feature.push_back(dou_vec);
            }


            for(int p = 0; p < 2; p++)
            {
                std::vector<double> dou_vec;

                bucket[(int)phase[45+5*p]/45] += mag[45+5*p];
                bucket[(int)phase[46+5*p]/45] += mag[46+5*p];
                bucket[(int)phase[47+5*p]/45] += mag[47+5*p];
                bucket[(int)phase[48+5*p]/45] += mag[48+5*p];

                bucket[(int)phase[54+5*p]/45] += mag[54+5*p];
                bucket[(int)phase[55+5*p]/45] += mag[55+5*p];
                bucket[(int)phase[56+5*p]/45] += mag[56+5*p];
                bucket[(int)phase[57+5*p]/45] += mag[57+5*p];

                bucket[(int)phase[63+5*p]/45] += mag[63+5*p];
                bucket[(int)phase[64+5*p]/45] += mag[64+5*p];
                bucket[(int)phase[65+5*p]/45] += mag[65+5*p];
                bucket[(int)phase[66+5*p]/45] += mag[66+5*p];

                bucket[(int)phase[72+5*p]/45] += mag[72+5*p];
                bucket[(int)phase[73+5*p]/45] += mag[73+5*p];
                bucket[(int)phase[74+5*p]/45] += mag[74+5*p];
                bucket[(int)phase[75+5*p]/45] += mag[75+5*p];


                for(int b = 0;b<8;b++)
                {
                    dou_vec.push_back(bucket[b]);
                }

                KeyPoint[scale][i].feature.push_back(dou_vec);
            }


        }

    }


}


void Mark_KeyPoint(KImageGray &igScale, double mag, double ori_deg, int pX, int pY){
    double dX_circle , dY_circle;
    int iX_circle, iY_circle;

    double radius = 0;

    int iX_dir, iY_dir;
    double dX_dir, dY_dir;

    double custom = 4;

    int tmp = (int)(mag / custom);
    if(tmp < 2 ){
        radius = 5;
    }
    else if(tmp < 5){
        radius = 7;
    }
    else if(tmp < 10){
        radius = 8;
    }
    else{
        radius = 10;
    }

    double theta_rad;
    //Draw circle
    for(int angle = 0; angle <= 360; angle+= 1){

        theta_rad = (double)angle * 0.01745329; // 0.017453..=>1 / 180 * 3.14592; degree to rad

        dX_circle = (double)pX - (double)(radius * cos(theta_rad));
        dY_circle = (double)pY - (double)(radius * sin(theta_rad));

        iX_circle = (int)dX_circle; iY_circle = (int)dY_circle;

        if(iX_circle > 0 && iY_circle > 0 && iX_circle < (int)igScale.Row() && iY_circle < (int)igScale.Col()){
           igScale[iX_circle][iY_circle] = 255;
        }
    }

    //Draw Direction
    for(int range = 0; range <= radius; range++){
        theta_rad = ori_deg * 0.01745329; // 0.017453..=>1 / 180 * 3.14592; degree to rad

        dX_dir = (double)pX - (double)(range * cos(theta_rad));
        dY_dir = (double)pY - (double)(range * sin(theta_rad));

        iX_dir = (int)dX_dir; iY_dir = (int)dY_dir;

        if(iX_dir > 0 && iY_dir > 0 && iX_dir < (int)igScale.Row() && iY_dir < (int)igScale.Col()){
           igScale[iX_dir][iY_dir] = 255;
        }

    }

}

















void Canny_Edge(KImageGray& imgSrc,double sigma)
{
    KImageGray icMain = imgSrc;


    int row = icMain.Row();
    int col = icMain.Col();

    //Gaussian Filter

    int Gau_filter_size = 8*sigma + 1; //filter size

    int size = std::sqrt(Gau_filter_size);

    //qDebug()<<size;

    double mask;


    for(int ii = size; ii<row-size;ii++){
        for(int jj = size; jj<col-size;jj++){

            mask = 0;


            for(int i = -size; i<size+1; i++){
                for(int j = -size; j<size+1;j++){
                    mask += std::exp(-0.5*((i*i+j*j)/(sigma*sigma)))*icMain[ii-i][jj-j];

                }
            }

            icMain[ii][jj] = (1/(2*M_PI*sigma*sigma))*mask;
        }
    }

    //Sobel Mask

    int dLow = 80;
    int dHigh = 120;
    double mag_x;
    double mag_y;
    double phase;

    double** magnitude = new double*[row]{0,};
    int** direction = new int*[row] {0,};
    for(int i = 0; i < row; i++)
    {
        magnitude[i] = new double[col]{0,};  // magnitude
        direction[i] = new int[col]{0,};
    }

    double Mask_w[3][3] = {{-1., 0., 1.},
                           {-2., 0., 2.},
                           {-1., 0. ,1.}};
    double Mask_h[3][3] = {{1., 2., 1.},
                           {0., 0., 0.},
                           {-1., -2. ,-1.}};

    //qDebug()<<row<< " "<<col;


    //qDebug()<<(unsigned char)((((int)(113/22.5)+1)>>1) & 0x00000003);

    for(int ii = 1; ii<row-1;ii++){
        for(int jj = 1; jj<col-1;jj++){

            mag_x = 0;
            mag_y = 0;

            for(int i = -1; i<2; i++){
                for(int j = -1; j<2;j++){
                    mag_x += (Mask_w[i+1][j+1]*icMain[ii+i][jj+j]);
                    mag_y += (Mask_h[i+1][j+1]*icMain[ii+i][jj+j]);
                }
            }

            magnitude[ii][jj] = fabs(mag_x) + fabs(mag_y);

            if(magnitude[ii][jj] > dLow)
            {
                //icMain2[ii][jj] = 255;
                phase = atan2(mag_x,mag_y)*180/M_PI;
                if(phase>360) phase -= 360;

                if(phase<0) phase += 180;

                direction[ii][jj] = (unsigned char)((((int)(phase/22.5)+1)>>1) & 0x00000003);


                //qDebug()<<ii<<" "<<jj<<" "<<phase<<" "<<direction[ii][jj];


            }
            else {
                magnitude[ii][jj] = 0;
                //icMain2[ii][jj] = 0;
            }

        }
    }


    // Non-Maxima Suppression
    //int nDx[4] = {0, 1, 1, 1};
    //int nDy[4] = {1, 1, 0, -1};
    int nDx[4] = {0, 1, 1, -1};
    int nDy[4] = {1, 1, 0, 1};

    double** buffer = new double*[row] {0,};
    int** realedge = new int*[row]{0,};
    for(int i = 0; i < row; i++)
    {
        buffer[i] = new double[col]{0,};
        realedge[i] = new int[col]{0,};
    }

    class HighEdge{
        public:
            int x;
            int y;
    };

    HighEdge h;
    std::stack<HighEdge> highedge;

    for(int ii = 1; ii<row-1;ii++){
        for(int jj = 1; jj<col-1;jj++){
            if(magnitude[ii][jj] == 0) continue;

            //if(magnitude[ii][jj] > magnitude[ii+nDx[direction[ii][jj]]][jj+nDy[direction[ii][jj]]] && magnitude[ii][jj] > magnitude[ii-nDx[direction[ii][jj]]][jj-nDy[direction[ii][jj]]])
            if(magnitude[ii][jj] > magnitude[ii+nDy[direction[ii][jj]]][jj+nDx[direction[ii][jj]]] && magnitude[ii][jj] > magnitude[ii-nDy[direction[ii][jj]]][jj-nDx[direction[ii][jj]]])
            {
                if(magnitude[ii][jj] > dHigh)
                {
                    h.x = ii;
                    h.y = jj;

                    highedge.push(h);
                    realedge[ii][jj] = 255;

                }

                //realedge[ii][jj] = 255;
                buffer[ii][jj] = magnitude[ii][jj];
            }
        }
    }


    for(int ii = 1; ii<row-1;ii++){
        for(int jj = 1; jj<col-1;jj++){
            //icMain2[ii][jj] = realedge[ii][jj];
        }
    }


    //Thresholding
    int x,y;
    while(!highedge.empty())
    {
        x = highedge.top().x;
        y = highedge.top().y;
        for(int i = -1; i < 2; i++)
        {
            for(int j = -1; j < 2; j++)
            {
                if(buffer[x+i][y+j] && buffer[x+i][y+j]<=dHigh)
                {
                    h.x = x+i;

                    h.y = y+j;
                    highedge.push(h);
                    realedge[x+i][y+j] = 255;
                    buffer[x+i][y+j] = 0;
                }
            }
        }
        highedge.pop();
    }

    for(int ii = 1; ii<row-1;ii++){
        for(int jj = 1; jj<col-1;jj++){
            imgSrc[ii][jj] = realedge[ii][jj];
        }
    }

    delete[] direction;
    delete[] magnitude;
}












UV finduv(Eigen::Matrix<double,25,2> mat_A, Eigen::Matrix<double,25,1> mat_b)
{

    Eigen::Matrix<double,2,1> mat_d;
    Eigen::Matrix<double,2,25> mat_AT;
    Eigen::Matrix<double,25,25> mat_AA_inverse;
    Eigen::Matrix<double,2,2> A_TA;
    Eigen::Matrix<double,2,2> A_TA_inverse;
    //A_T
    mat_AT = mat_A.transpose();

    //A_T*A
    A_TA = mat_AT*mat_A;

    A_TA_inverse=A_TA.inverse();

    mat_d = A_TA_inverse*mat_AT*mat_b; //2x2 2x25 25x1

    //qDebug()<<mat_d(1,0);

    UV a;

    a.u = mat_d(0,0);
    a.v = mat_d(1,0);

    //qDebug()<<a.u<<a.v;

    return a;

}

void Optical_Flow(KImageGray Image0, KImageGray Image1, UV** mat_uv)
{
    int row = Image0.Row();
    int col = Image0.Col();

    double** mat_It = new double*[row]{0,};
    double** mat_Ix = new double*[row]{0,};
    double** mat_Iy = new double*[row]{0,};

    //UV** mat_uv = new UV*[row]{0,};

    //KImageGray temp(row,col);

    //int dx[4] = {-1,1,-1,1};
    //int dy[4]= {1,1,-1,-1};

    int dx[4] = {0,1,-1,0};
    int dy[4]= {1,0,0,-1};



    //double mat_b[1][25];
    Eigen::Matrix<double, 25, 1> mat_b;
    //double mat_A[25][2];
    Eigen::Matrix<double, 25, 2> mat_A;




    for(int i = 0; i < row; i++)
    {
        mat_It[i] = new double[col]{0,};
        mat_Ix[i] = new double[col]{0,};
        mat_Iy[i] = new double[col]{0,};
        //mat_uv[i] = new UV[col]{{0,0},};
    }


    for(int i=0; i < row - 1; i++)
    {
        for(int j=0; j<col - 1; j++)
        {
            mat_It[i][j] = Image1[i][j]-Image0[i][j];
            mat_Ix[i][j] = Image1[i+1][j+1]*dx[0]+Image1[i+1][j+1]*dx[1]+Image1[i][j]*dx[2]+Image1[i+1][j]*dx[3];
            mat_Iy[i][j] = Image1[i+1][j+1]*dy[0]+Image1[i+1][j+1]*dy[1]+Image1[i][j]*dy[2]+Image1[i+1][j]*dy[3];
        }
    }

    //qDebug()<<"mat_A.size()";


    for(int i=2;i<row-2;i++)
    {
        for(int j=2;j<col-2;j++)
        {
           for(int l=-2; l<3; l++)
           {
               for(int m=-2; m<3; m++)
               {
                   mat_b(5*(l+2) + m+2,0) = -mat_It[i+l][j+m];

                   mat_A(5*(l+2) + m+2,0) = mat_Ix[i+l][j+m];
                   mat_A(5*(l+2) + m+2,1) = mat_Iy[i+l][j+m];
               }
           }

           mat_uv[i][j] = finduv(mat_A,mat_b);
        }
    }


    //return temp;

}


void Mark_KeyPoint_Line(KImageGray &igScale, double mag, double ori_deg, int pX, int pY){
    double dX_circle , dY_circle;
    int iX_circle, iY_circle;

    double radius = 0;

    int iX_dir, iY_dir;
    double dX_dir, dY_dir;

    double custom = 4;

    int tmp = (int)(mag / custom);
    if(tmp < 2 ){
        radius = 5;
    }
    else if(tmp < 5){
        radius = 7;
    }
    else if(tmp < 10){
        radius = 8;
    }
    else{
        radius = 10;
    }

    double theta_rad;
    //Draw circle
    for(int angle = 0; angle <= 360; angle+= 1){

        theta_rad = (double)angle * 0.01745329; // 0.017453..=>1 / 180 * 3.14592; degree to rad

        dX_circle = (double)pX - (double)(radius * cos(theta_rad));
        dY_circle = (double)pY - (double)(radius * sin(theta_rad));

        iX_circle = (int)dX_circle; iY_circle = (int)dY_circle;

        if(iX_circle > 0 && iY_circle > 0 && iX_circle < (int)igScale.Row() && iY_circle < (int)igScale.Col()){
           //igScale[iX_circle][iY_circle] = 255;
        }
    }

    // !!!Main!!! -> Draw Line


    //Draw Direction
    for(int range = 0; range <= mag*2; range++){
        theta_rad = ori_deg * 0.01745329; // 0.017453..=>1 / 180 * 3.14592; degree to rad

        dX_dir = (double)pX - (double)(range * cos(theta_rad));
        dY_dir = (double)pY - (double)(range * sin(theta_rad));

        iX_dir = (int)dX_dir; iY_dir = (int)dY_dir;

        if(iX_dir > 0 && iY_dir > 0 && iX_dir < (int)igScale.Row() && iY_dir < (int)igScale.Col()){
           igScale[iX_dir][iY_dir] = 255;
        }

    }

    // !!!Main!!!

}


void Draw(KImageGray& Img, UV** mat_uv)
{

    double mag, phase;
    double u,v;

    int row = Img.Row();
    int col = Img.Col();


    for(int i = 0 ; i < row; i+=10 )
    {
        for(int j = 0; j < col; j+=5)
        {

            u = mat_uv[i][j].u;
            v = mat_uv[i][j].v;

            mag = sqrt(u*u+v*v);
            phase = atan2(u,v)*180/M_PI;

            if(phase>360) phase -= 360;

            if(phase<0) phase += 180;



            Mark_KeyPoint_Line(Img,mag, phase, i, j);

        }
    }

}





