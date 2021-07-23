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
#include <queue>
#include <fstream>
#include "eigen-3.3.8/Eigen/Dense"


#include "mainframe.h"
#include "ui_mainframe.h"
#include "imageform.h"
#include "dahyeon.h"


MainFrame::MainFrame(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::MainFrame)
{
    ui->setupUi(this);

    _plpImageForm       = new KPtrList<ImageForm*>(10,false,false);
    _q_pFormFocused     = 0;

    //객체 맴버의 초기화

    //get a current directory
    char st[100];
    GetCurrentDirectoryA(100,st);

    //리스트 출력창을 안보이게
    ui->listWidget->setVisible(false);
    this->adjustSize();

    //UI 활성화 갱신
    UpdateUI();
}

MainFrame::~MainFrame()
{
    delete ui;
    delete _plpImageForm;

}

void MainFrame::CloseImageForm(ImageForm *pForm)
{
    //ImageForm 포인터 삭제
    _plpImageForm->Remove(pForm);

    //활성화 ImageForm 초기화
    _q_pFormFocused     = 0;

    //관련 객체 삭제

    //UI 활성화 갱신
    UpdateUI();
}

void MainFrame::UpdateUI()
{

    if(ui->tabWidget->currentIndex() == 0)
    {
        ui->buttonSepiaTone->setEnabled( _q_pFormFocused &&  _q_pFormFocused->ID() == "OPEN" );

    }
    else if(ui->tabWidget->currentIndex() == 1)
    {

    }
    else if(ui->tabWidget->currentIndex() == 2)
    {

    }

}

void MainFrame::OnMousePos(const int &nX, const int &nY, ImageForm* q_pForm)
{

    UpdateUI();
}

void MainFrame::closeEvent(QCloseEvent* event)
{
    //생성된 ImageForm을 닫는다.
    for(int i=_plpImageForm->Count()-1; i>=0; i--)
        _plpImageForm->Item(i)->close();

    //리스트에서 삭제한다.
    _plpImageForm->RemoveAll();
}


void MainFrame::on_buttonOpen_clicked()
{
    //이미지 파일 선택
    QFileDialog::Options    q_Options   =  QFileDialog::DontResolveSymlinks  | QFileDialog::DontUseNativeDialog; // | QFileDialog::ShowDirsOnly
    QString                 q_stFile    =  QFileDialog::getOpenFileName(this, tr("Select a Image File"),  "./data", "Image Files(*.bmp *.ppm *.pgm *.png)",0, q_Options);

    if(q_stFile.length() == 0)
        return;

    //이미지 출력을 위한 ImageForm 생성
    ImageForm*              q_pForm   = new ImageForm(q_stFile, "OPEN", this);

    _plpImageForm->Add(q_pForm);
    q_pForm->show();

    //UI 활성화 갱신
    UpdateUI();
}

void MainFrame::on_buttonDeleteContents_clicked()
{

    //생성된 ImageForm을 닫는다.
    for(int i=_plpImageForm->Count()-1; i>=0; i--)
        _plpImageForm->Item(i)->close();

    //리스트에서 삭제한다.
    _plpImageForm->RemoveAll();
}

void MainFrame::on_buttonSepiaTone_clicked()
{
    KImageColor   icMain;


    //포커스 된 ImageForm으로부터 영상을 가져옴

    if(_q_pFormFocused != 0 && _q_pFormFocused->ImageColor().Address() &&  _q_pFormFocused->ID() == "OPEN")
        icMain = _q_pFormFocused->ImageColor();
    else
        return;

    KImageGray img_Hue(icMain.Row(),icMain.Col()),
               img_Sat(icMain.Row(),icMain.Col()),
               img_Val(icMain.Row(),icMain.Col());

    //hue, saturation 값 가져오기
    double dHue2 = ui->spinHue->value();//text().toDouble();
    double dSat2 = ui->spinSaturation->value(); // 위와 같은 식으로

    //
    //icMain 변환
    //:

    double dMin, dMax, dDiff;
    double dVal;
    double C,X,m;
    double R,G,B;
    double dHue,dSat;

    for(int i=icMain.Row(),ii=0; i; i--,ii++)
        for(int j=icMain.Col(),jj=0; j; j--,jj++)
        {
            dMin  = std::min(std::min(icMain[ii][jj].r,icMain[ii][jj].g),icMain[ii][jj].b);
            dMax  = std::max(std::max(icMain[ii][jj].r,icMain[ii][jj].g),icMain[ii][jj].b);
            dDiff = dMax - dMin;

            //value
            dVal = dMax/255.0;
            img_Val[ii][jj]= dVal*255;

            //saturation
            dSat = dDiff/dMax;

            //hue
            if(dMax == (double)(icMain[ii][jj].r)) //maximun = red
                dHue = 60.0 * (double)(icMain[ii][jj].g-icMain[ii][jj].b)/dDiff;
            else if(dMax == (double)(icMain[ii][jj].g))
                dHue = 60.0 * (double)(icMain[ii][jj].b-icMain[ii][jj].r)/dDiff + 120.0;
            else
                dHue = 60.0 * (double)(icMain[ii][jj].r-icMain[ii][jj].g)/dDiff + 240.0;

            if(dHue == 360.0)
                dHue = 0.0;

            img_Hue[ii][jj]=dHue/360.0 * 255;
            img_Sat[ii][jj]=dSat * 255.0;

            dHue = dHue2;
            dSat = dSat2;

            C = dVal*dSat;
            X = C*(1-std::fabs(fmod(dHue/60,2)-1));
            m = dVal - C;

            if(dHue>=0&&dHue<60)
            {
                R = C;
                G = X;
                B = 0;
            }
            else if (dHue>=60&&dHue<120)
            {
                R = X;
                G = C;
                B = 0;
            }
            else if (dHue>=120&&dHue<180)
            {
                R = 0;
                G = C;
                B = X;
            }
            else if (dHue>=180&&dHue<240)
            {
                R = 0;
                G = X;
                B = C;
            }
            else if (dHue>=240&&dHue<300)
            {
                R = X;
                G = 0;
                B = C;
            }
            else if (dHue>=300&&dHue<360)
            {
                R = C;
                G = 0;
                B = X;
            }


            icMain[ii][jj].r=(R+m)*255.0;
            icMain[ii][jj].g=(G+m)*255.0;
            icMain[ii][jj].b=(B+m)*255.0;

        }

    //출력을 위한 ImageForm 생성
    //sepia
    ImageForm*  q_pForm = new ImageForm(icMain, "Sepia Tone", this);
    _plpImageForm->Add(q_pForm);
    q_pForm->show();


    //HSV
    ImageForm*  q_pHue = new ImageForm(img_Hue,"Hue",this);
    _plpImageForm->Add(q_pHue);
    q_pHue->show();

    ImageForm*  q_pSat = new ImageForm(img_Sat,"Sat",this);
    _plpImageForm->Add(q_pSat);
    q_pSat->show();

    ImageForm*  q_pVal = new ImageForm(img_Val,"Val",this);
    _plpImageForm->Add(q_pVal);
    q_pVal->show();


    //UI 활성화 갱신
    UpdateUI();
}


void MainFrame::on_buttonShowList_clicked()
{
    static int nWidthOld = ui->tabWidget->width();

    if(ui->listWidget->isVisible())
    {
        nWidthOld = ui->listWidget->width();
        ui->listWidget->hide();
        this->adjustSize();
    }
    else
    {
        ui->listWidget->show();
        QRect q_rcWin = this->geometry();

        this->setGeometry(q_rcWin.left(), q_rcWin.top(), q_rcWin.width()+nWidthOld, q_rcWin.height());
    }
}

void MainFrame::on_spinHue_valueChanged(const QString &arg1)
{

}

void MainFrame::on_spinSaturation_valueChanged(const QString &arg1)
{


}

void MainFrame::on_spinHue_textChanged(const QString &arg1)
{
    ImageForm* q_pForm = 0;

        for(int i=0; i<_plpImageForm->Count();i++)
            if((*_plpImageForm)[i]->ID()=="Sepia Tone")
            {
                q_pForm = (*_plpImageForm)[i];
                break;
            }



        KImageColor icMain = q_pForm->ImageColor();

        KImageGray img_Hue(icMain.Row(),icMain.Col()),
                   img_Sat(icMain.Row(),icMain.Col()),
                   img_Val(icMain.Row(),icMain.Col());

        //hue, saturation 값 가져오기
        double dHue2 = ui->spinHue->value();//text().toDouble();
        double dSat2 = ui->spinSaturation->value(); // 위와 같은 식으로

        //
        //icMain 변환
        //:

        //double dValue = max(icMain[i][j].r,icMain.[i][j].g,icMain.[i][j].b)
        //double dValue <= icMain.[i][j].r
        //                 icMain.[i][j].g
        //                  icMain.[i][j].b

        //icMain[i][j].r,g,b <= dHue, dSat, dValue;

        double dMin, dMax, dDiff;
        double dVal;
        double C,X,m;
        double R,G,B;
        double dHue,dSat;

        for(int i=icMain.Row(),ii=0; i; i--,ii++)
            for(int j=icMain.Col(),jj=0; j; j--,jj++)
            {
                dMin  = std::min(std::min(icMain[ii][jj].r,icMain[ii][jj].g),icMain[ii][jj].b);
                dMax  = std::max(std::max(icMain[ii][jj].r,icMain[ii][jj].g),icMain[ii][jj].b);
                dDiff = dMax - dMin;

                //value
                dVal = dMax/255.0;
                img_Val[ii][jj]= dVal*255;

                //saturation
                dSat = dDiff/dMax;

                //hue
                if(dMax == (double)(icMain[ii][jj].r)) //maximun = red
                    dHue = 60.0 * (double)(icMain[ii][jj].g-icMain[ii][jj].b)/dDiff;
                else if(dMax == (double)(icMain[ii][jj].g))
                    dHue = 60.0 * (double)(icMain[ii][jj].b-icMain[ii][jj].r)/dDiff + 120.0;
                else
                    dHue = 60.0 * (double)(icMain[ii][jj].r-icMain[ii][jj].g)/dDiff + 240.0;

                if(dHue == 360.0)
                    dHue = 0.0;

                img_Hue[ii][jj]=dHue/360.0 * 255;
                img_Sat[ii][jj]=dSat * 255.0;

                dHue = dHue2;
                dSat = dSat2;

                C = dVal*dSat;
                X = C*(1-std::fabs(fmod(dHue/60,2)-1));
                m = dVal - C;

                if(dHue>=0&&dHue<60)
                {
                    R = C;
                    G = X;
                    B = 0;
                }
                else if (dHue>=60&&dHue<120)
                {
                    R = X;
                    G = C;
                    B = 0;
                }
                else if (dHue>=120&&dHue<180)
                {
                    R = 0;
                    G = C;
                    B = X;
                }
                else if (dHue>=180&&dHue<240)
                {
                    R = 0;
                    G = X;
                    B = C;
                }
                else if (dHue>=240&&dHue<300)
                {
                    R = X;
                    G = 0;
                    B = C;
                }
                else if (dHue>=300&&dHue<360)
                {
                    R = C;
                    G = 0;
                    B = X;
                }


                icMain[ii][jj].r=(R+m)*255.0;
                icMain[ii][jj].g=(G+m)*255.0;
                icMain[ii][jj].b=(B+m)*255.0;

            }

        if(q_pForm)
            q_pForm->Update(icMain);
        else
        {
            q_pForm = new ImageForm(icMain, "Sepia Tone", this);
            _plpImageForm->Add(q_pForm);
            q_pForm->show();
        }




}

void MainFrame::on_spinSaturation_textChanged(const QString &arg1)
{
    ImageForm* q_pForm = 0;

        for(int i=0; i<_plpImageForm->Count();i++)
            if((*_plpImageForm)[i]->ID()=="Sepia Tone")
            {
                q_pForm = (*_plpImageForm)[i];
                break;
            }

        KImageColor icMain=q_pForm->ImageColor();

        KImageGray img_Hue(icMain.Row(),icMain.Col()),
                   img_Sat(icMain.Row(),icMain.Col()),
                   img_Val(icMain.Row(),icMain.Col());

        //hue, saturation 값 가져오기
        double dHue2 = ui->spinHue->value();//text().toDouble();
        double dSat2 = ui->spinSaturation->value(); // 위와 같은 식으로

        //
        //icMain 변환
        //:

        //double dValue = max(icMain[i][j].r,icMain.[i][j].g,icMain.[i][j].b)
        //double dValue <= icMain.[i][j].r
        //                 icMain.[i][j].g
        //                  icMain.[i][j].b

        //icMain[i][j].r,g,b <= dHue, dSat, dValue;

        double dMin, dMax, dDiff;
        double dVal;
        double C,X,m;
        double R,G,B;
        double dHue,dSat;

        for(int i=icMain.Row(),ii=0; i; i--,ii++)
            for(int j=icMain.Col(),jj=0; j; j--,jj++)
            {
                dMin  = std::min(std::min(icMain[ii][jj].r,icMain[ii][jj].g),icMain[ii][jj].b);
                dMax  = std::max(std::max(icMain[ii][jj].r,icMain[ii][jj].g),icMain[ii][jj].b);
                dDiff = dMax - dMin;

                //value
                dVal = dMax/255.0;
                img_Val[ii][jj]= dVal*255;

                //saturation
                dSat = dDiff/dMax;

                //hue
                if(dMax == (double)(icMain[ii][jj].r)) //maximun = red
                    dHue = 60.0 * (double)(icMain[ii][jj].g-icMain[ii][jj].b)/dDiff;
                else if(dMax == (double)(icMain[ii][jj].g))
                    dHue = 60.0 * (double)(icMain[ii][jj].b-icMain[ii][jj].r)/dDiff + 120.0;
                else
                    dHue = 60.0 * (double)(icMain[ii][jj].r-icMain[ii][jj].g)/dDiff + 240.0;

                if(dHue == 360.0)
                    dHue = 0.0;

                img_Hue[ii][jj]=dHue/360.0 * 255;
                img_Sat[ii][jj]=dSat * 255.0;

                dHue = dHue2;
                dSat = dSat2;

                C = dVal*dSat;
                X = C*(1-std::fabs(fmod(dHue/60,2)-1));
                m = dVal - C;

                if(dHue>=0&&dHue<60)
                {
                    R = C;
                    G = X;
                    B = 0;
                }
                else if (dHue>=60&&dHue<120)
                {
                    R = X;
                    G = C;
                    B = 0;
                }
                else if (dHue>=120&&dHue<180)
                {
                    R = 0;
                    G = C;
                    B = X;
                }
                else if (dHue>=180&&dHue<240)
                {
                    R = 0;
                    G = X;
                    B = C;
                }
                else if (dHue>=240&&dHue<300)
                {
                    R = X;
                    G = 0;
                    B = C;
                }
                else if (dHue>=300&&dHue<360)
                {
                    R = C;
                    G = 0;
                    B = X;
                }


                icMain[ii][jj].r=(R+m)*255.0;
                icMain[ii][jj].g=(G+m)*255.0;
                icMain[ii][jj].b=(B+m)*255.0;

            }

        if(q_pForm)
            q_pForm->Update(icMain);
        else
        {
            q_pForm = new ImageForm(icMain, "Sepia Tone", this);
            _plpImageForm->Add(q_pForm);
            q_pForm->show();
        }

}

void MainFrame::on_buttonOtsu_clicked()
{

    KImageGray igImg;

    //포커스 된 ImageForm으로부터 영상을 가져옴
    if(_q_pFormFocused != 0 && _q_pFormFocused->ImageGray().Address() &&  _q_pFormFocused->ID() == "OPEN")
    {
        igImg = _q_pFormFocused->ImageGray();

    }
    else
        return;
/*
    if(_q_pFormFocused != 0 && _q_pFormFocused->ImageColor().Address() &&  _q_pFormFocused->ID() == "OPEN")
    {
        imglabel = _q_pFormFocused->ImageColor();
    }
    else
        return;
*/


    //KImageGray labeling;

    //To get its histogramming
    int histo[256] = {0, };
    double sigmab[256] = {0,};
    double N = igImg.Size() ,p,q1,q2,u,u1=0,u2=0,temp_q1,T;
    double max = -1;


    //histogram
    for(unsigned int i=0; i<igImg.Row(); i++){
        for(unsigned int j=0; j<igImg.Col(); j++)
        {
            histo[igImg[i][j]] += 1;

        }
    }


    for(int t = 0; t<255;t++){
       p = (double)histo[t] / (double) N;
       u += ((double)t)*p;
    }

    //Find T
    // t=0
    p = (double)histo[0]/(double)N;

    //qDebug() << p;

    q1 = p;
    q2 = 1.0 - q1;
    sigmab[0]=q1*q2*(u1-u2)*(u1-u2);

    //qDebug() << sigmab[0];


    for(int t=1;t<256;t++)
    {
        p = (double)histo[t]/(double)N;
        temp_q1 = q1;
        q1 += p;
        q2 = 1.0 - q1;
        if(q1 ==0)
        {
            u1=0;
            u2 = (u - q1 * u1) / (1.0 - q1);
        }
        else if(q2 == 0)
        {
            u2=0;
            u1 = ((temp_q1 * u1) + ((double)(t)*p)) / (q1);
        }
        else
        {
            u1 = ((temp_q1 * u1) + ((double)(t)*p)) / (q1);
            u2 = (u - q1 * u1) / (1.0 - q1);
        }

        sigmab[t] = q1*q2*(u1-u2)*(u1-u2);

        //qDebug() << "p : " << p <<" , "<< t << " : " << sigmab[t];
    }

    for(int i = 0; i < 256; i++)
    {
        if(max < sigmab[i])
        {
            max = sigmab[i];
            T = i;
            //qDebug() <<" T : " <<  T;
        }
    }



    //Thresholding
    for(int i=0; i<igImg.Row(); i++)
    {
        for(int j=0; j<igImg.Col(); j++)
        {
            if(igImg[i][j] > T) igImg[i][j] = 255; //foreground
            else igImg[i][j] = 0; //background
        }
    }


    //Image Labeling

    KImageGray igImg2 = igImg;
    int label[igImg2.Row()][igImg2.Col()];
    std::fill(&label[0][0], &label[igImg2.Row()-1][igImg2.Col()-1], 0);

    int label_num = 1;

    for(int i = 1; i<igImg2.Row()-1;i++)
    {
        for(int j = 1; j<igImg2.Col()-1;j++)
        {
            if(igImg2[i][j] == 255) //일단 내가 전경
            {
                if( (igImg2[i-1][j-1] == 0 && igImg2[i-1][j] == 0 && igImg2[i][j-1] == 0) || (label[i-1][j-1] == 0 && label[i-1][j] == 0 && label[i][j-1] == 0)) // 대각선 위 오른쪽이 전부 배경
                {
                    label[i][j] = label_num;
                    label_num++;
                }
                else if(igImg2[i-1][j-1] == 0 && igImg2[i-1][j] == 255 && igImg2[i][j-1] == 0) //위쪽만 전경
                {
                    label[i][j] = label[i-1][j];
                }
                else if(igImg2[i-1][j-1] == 0 && igImg2[i-1][j] == 0 && igImg2[i][j-1] == 255) //왼쪽만 전경
                {
                    label[i][j] = label[i][j-1];
                }
                else if(igImg2[i-1][j-1] == 255 && igImg2[i-1][j] == 0 && igImg2[i][j-1] == 0) //대각선만 전경
                {
                    label[i][j] = label[i-1][j-1];
                }
                else if(igImg2[i-1][j-1] == 255 && igImg2[i-1][j] == 0 && igImg2[i][j-1] == 255) //대각선이랑 왼쪽만 전경
                {
                    label[i][j] = label[i][j-1];
                    label[i-1][j] = label[i][j-1];
                }
                else if(igImg2[i-1][j-1] == 255 && igImg2[i-1][j] == 255 && igImg2[i][j-1] == 0) //대각선이랑 위쪽만 전경
                {
                    label[i][j] = label[i-1][j];
                    label[i][j-1] = label[i-1][j];
                }
                else if(igImg2[i-1][j-1] == 0 && igImg2[i-1][j] == 255 && igImg2[i][j-1] == 255) //왼쪽이랑 위쪽만 전경
                {
                    label[i][j] = label[i-1][j];
                    label[i][j-1] = label[i-1][j];
                }
                else if(igImg2[i-1][j-1] == 255 && igImg2[i-1][j] == 255 && igImg2[i][j-1] == 255) // 전부 전경
                {
                    label[i][j] = label[i-1][j]; //위쪽 , 왼쪽, 대각선 순 기준
                    label[i-1][j-1] = label[i-1][j];
                    label[i][j-1] = label[i-1][j];
                }
            }
        }
    }
    for(int i = 1; i<igImg2.Row()-1;i++)
    {
        for(int j = igImg2.Col()-2; j > 1;j--)
        {
            if(igImg2[i][j] == 255) //일단 내가 전경
            {
                if(igImg2[i-1][j-1] == 0 && igImg2[i-1][j] == 0 && igImg2[i][j-1] == 255) //왼쪽만 전경
                {
                    label[i][j-1] = label[i][j];
                }
                else if(igImg2[i-1][j-1] == 0 && igImg2[i-1][j] == 255 && igImg2[i][j-1] == 0) //위쪽만 전경
                {
                    label[i][j] = label[i-1][j];
                }
                else if(igImg2[i-1][j-1] == 255 && igImg2[i-1][j] == 255 && igImg2[i][j-1] == 255) // 전부 전경
                {
                    label[i][j] = label[i-1][j]; //위쪽 , 왼쪽, 대각선 순 기준
                    label[i-1][j-1] = label[i-1][j];
                    label[i][j-1] = label[i-1][j];
                }

            }
        }
    }


    KImageColor imglabel = KImageColor(igImg2.Row(),igImg2.Col());

    for(unsigned int ii = 1;ii<igImg2.Row();ii++)
              for(unsigned int jj = 1;jj<igImg2.Col();jj++)
          {
              imglabel[ii][jj].b = label[ii][jj];
              imglabel[ii][jj].r = (5* label[ii][jj]) % 256 ;
              imglabel[ii][jj].g = (10* label[ii][jj]) % 256;
          }



    //show
    ImageForm*  q_pForm = new ImageForm(igImg, "Otsu's Thresholding", this);
    _plpImageForm->Add(q_pForm);
    q_pForm->show();

    ImageForm*  q_pForm2 = new ImageForm(imglabel, "Labeling", this);
    _plpImageForm->Add(q_pForm2);
    q_pForm2->show();



    //UI 활성화 갱신
    UpdateUI();



}

void MainFrame::on_button_dil_n_ero_clicked()
{



}

void MainFrame::on_button3x3_clicked()
{
    ImageForm* q_pForm = 0;

        for(int i=0; i<_plpImageForm->Count();i++)
            if((*_plpImageForm)[i]->ID()=="Otsu's Thresholding")
            {
                q_pForm = (*_plpImageForm)[i];
                break;
            }



        KImageGray img_dil = q_pForm->ImageGray();
        KImageGray img_ero = q_pForm->ImageGray();
        KImageGray img_tmp = q_pForm->ImageGray();

         //qDebug() << "1";

        //Dilation
        for(unsigned int i = 1;i<img_dil.Row()-1;i++)
        {
            for(unsigned int j = 1;j<img_dil.Col()-1;j++)
            {
                if(img_tmp[i-1][j-1] == 255 || img_tmp[i-1][j] == 255 || img_tmp[i-1][j+1] == 255 || img_tmp[i][j-1] == 255 || img_tmp[i][j+1] == 255 || img_tmp[i+1][j-1] == 255 || img_tmp[i+1][j] == 255 || img_tmp[i+1][j+1] == 255 )
                    img_dil[i][j] = 255;

                //qDebug() << "i : " << i << "j : " << j;
            }

        }

        //Eroison
        for(unsigned int i = 1;i<img_ero.Row()-1;i++)
        {
            for(unsigned int j = 1;j<img_ero.Col()-1;j++)
            {
                if(img_tmp[i-1][j-1] == 0 || img_tmp[i-1][j] == 0 || img_tmp[i-1][j+1] == 0 || img_tmp[i][j-1] == 0 || img_tmp[i][j+1] == 0 || img_tmp[i+1][j-1] == 0 || img_tmp[i+1][j] == 0 || img_tmp[i+1][j+1] == 0 )
                    img_ero[i][j] = 0;

                //qDebug() << "i : " << i << "j : " << j;
            }

        }



        ImageForm*  q_dil = new ImageForm(img_dil, "Dilation3x3", this);
        _plpImageForm->Add(q_dil);
        q_dil->show();

        ImageForm*  q_ero = new ImageForm(img_ero, "Eroison3x3", this);
        _plpImageForm->Add(q_ero);
        q_ero->show();
}

void MainFrame::on_button5x5_clicked()
{
    ImageForm* q_pForm = 0;

        for(int i=0; i<_plpImageForm->Count();i++)
            if((*_plpImageForm)[i]->ID()=="Otsu's Thresholding")
            {
                q_pForm = (*_plpImageForm)[i];
                break;
            }



        KImageGray img_dil = q_pForm->ImageGray();
        KImageGray img_ero = q_pForm->ImageGray();
        KImageGray img_tmp = q_pForm->ImageGray();

         //qDebug() << "1";

        //Dilation
        for(unsigned int i = 2;i<img_dil.Row()-2;i++)
        {
            for(unsigned int j = 2;j<img_dil.Col()-2;j++)
            {
                if(img_tmp[i-2][j-2] == 255 || img_tmp[i-2][j-1] == 255|| img_tmp[i-2][j] == 255 || img_tmp[i-2][j+1] == 255 || img_tmp[i-2][j+2] == 255 ||
                   img_tmp[i-1][j-2] == 255 || img_tmp[i-1][j-1] == 255|| img_tmp[i-1][j] == 255 || img_tmp[i-1][j+1] == 255 || img_tmp[i-1][j+2] == 255 ||
                   img_tmp[i][j-2] == 255 || img_tmp[i][j-1] == 255 || img_tmp[i][j+1] == 255 || img_tmp[i][j+2] == 255 ||
                   img_tmp[i+1][j-2] == 255 || img_tmp[i+1][j-1] == 255|| img_tmp[i+1][j] == 255 || img_tmp[i+1][j+1] == 255 || img_tmp[i+1][j+2] == 255 ||
                   img_tmp[i+2][j-2] == 255 || img_tmp[i+2][j-1] == 255|| img_tmp[i+2][j] == 255 || img_tmp[i+2][j+1] == 255 || img_tmp[i+2][j+2] == 255)
                    img_dil[i][j] = 255;

                //qDebug() << "i : " << i << "j : " << j;
            }

        }

        //Eroison
        for(unsigned int i = 2;i<img_ero.Row()-2;i++)
        {
            for(unsigned int j = 2;j<img_ero.Col()-2;j++)
            {
                if(img_tmp[i-2][j-2] == 0 || img_tmp[i-2][j-1] == 0|| img_tmp[i-2][j] == 0 || img_tmp[i-2][j+1] == 0 || img_tmp[i-2][j+2] == 0 ||
                   img_tmp[i-1][j-2] == 0 || img_tmp[i-1][j-1] == 0|| img_tmp[i-1][j] == 0 || img_tmp[i-1][j+1] == 0 || img_tmp[i-1][j+2] == 0 ||
                   img_tmp[i][j-2] == 0 || img_tmp[i][j-1] == 0 || img_tmp[i][j+1] == 0 || img_tmp[i][j+2] == 0 ||
                   img_tmp[i+1][j-2] == 0 || img_tmp[i+1][j-1] == 0|| img_tmp[i+1][j] == 0 || img_tmp[i+1][j+1] == 0 || img_tmp[i+1][j+2] == 0 ||
                   img_tmp[i+2][j-2] == 0 || img_tmp[i+2][j-1] == 0|| img_tmp[i+2][j] == 0 || img_tmp[i+2][j+1] == 0 || img_tmp[i+2][j+2] == 0 )
                    img_ero[i][j] = 0;

                //qDebug() << "i : " << i << "j : " << j;
            }

        }



        ImageForm*  q_dil = new ImageForm(img_dil, "Dilation5x5", this);
        _plpImageForm->Add(q_dil);
        q_dil->show();

        ImageForm*  q_ero = new ImageForm(img_ero, "Eroison5x5", this);
        _plpImageForm->Add(q_ero);
        q_ero->show();
}

void MainFrame::on_buttonHEQ_clicked()
{
    KImageColor icMain;

    //포커스 된 ImageForm으로부터 영상을 가져옴
    if(_q_pFormFocused != 0 && _q_pFormFocused->ImageColor().Address() &&  _q_pFormFocused->ID() == "OPEN")
    {
        icMain = _q_pFormFocused->ImageColor();
    }
    else
        return;


    //To get its histogramming
    int histo_R[256] = {0, };
    int histo_G[256] = {0, };
    int histo_B[256] = {0, };

    double P_R[256] = {0, };
    double P_G[256] = {0, };
    double P_B[256] = {0, };

    double T_R[256] = {0, };
    double T_G[256] = {0, };
    double T_B[256] = {0, };

    double r_R[256] = {0, };
    double r_G[256] = {0, };
    double r_B[256] = {0, };


    //histogram
    for(unsigned int i=0; i<icMain.Row(); i++){
        for(unsigned int j=0; j<icMain.Col(); j++)
        {
            histo_R[icMain[i][j].r] += 1;
            histo_G[icMain[i][j].g] += 1;
            histo_B[icMain[i][j].b] += 1;

        }
    }

    for(unsigned int t=0; t<256; t++){

        P_R[t] = (double)histo_R[t]/(double)icMain.Size();
        P_G[t] = (double)histo_G[t]/(double)icMain.Size();
        P_B[t] = (double)histo_B[t]/(double)icMain.Size();

        //qDebug() << t << " : " << P_R[t];

    }

    T_R[0] = P_R[0];
    T_G[0] = P_G[0];
    T_B[0] = P_B[0];

    // 1/255 * r = T_R[r]

    for(unsigned int r=1; r<256; r++){
        T_R[r] = T_R[r-1] + P_R[r];
        r_R[r] = T_R[r] * 255;
        T_B[r] = T_G[r-1] + P_G[r];
        r_G[r] = T_G[r] * 255;
        T_G[r] = T_B[r-1] + P_B[r];
        r_B[r] = T_B[r] * 255;

        //qDebug() << r << " : " << T_R[r];

    }


    for(unsigned int i=0; i<icMain.Row(); i++){
        for(unsigned int j=0; j<icMain.Col(); j++)
        {
            icMain[i][j].r = r_R[icMain[i][j].r];
            icMain[i][j].g = r_R[icMain[i][j].g];
            icMain[i][j].b = r_R[icMain[i][j].b];
        }
    }

    ImageForm*  q_pForm2 = new ImageForm(icMain, "HEQ", this);
    _plpImageForm->Add(q_pForm2);
    q_pForm2->show();
}



void MainFrame::on_buttonHMA_clicked()
{
    //포커스 된 ImageForm으로부터 영상을 가져옴
    KImageColor Source;

    //포커스 된 ImageForm으로부터 영상을 가져옴
    if(_q_pFormFocused != 0 && _q_pFormFocused->ImageColor().Address() &&  _q_pFormFocused->ID() == "OPEN")
    {
        Source = _q_pFormFocused->ImageColor();
    }
    else
        return;


    ImageForm* q_pForm = 0;

        for(int i=0; i<_plpImageForm->Count();i++)
            if((*_plpImageForm)[i]->ID()=="target")
            {
                q_pForm = (*_plpImageForm)[i];
                break;
            }

    KImageColor Target=q_pForm->ImageColor();



    //get histogram
    int histo_S_R[256] = {0, };
    int histo_S_G[256] = {0, };
    int histo_S_B[256] = {0, };

    int histo_T_R[256] = {0, };
    int histo_T_G[256] = {0, };
    int histo_T_B[256] = {0, };

    for(unsigned int i=0; i<Source.Row(); i++){
        for(unsigned int j=0; j<Source.Col(); j++)
        {
            histo_S_R[Source[i][j].r] += 1;
            histo_S_G[Source[i][j].g] += 1;
            histo_S_B[Source[i][j].b] += 1;
        }
    }

    for(unsigned int i=0; i<Target.Row(); i++){
        for(unsigned int j=0; j<Target.Col(); j++)
        {
            histo_T_R[Target[i][j].r] += 1;
            histo_T_G[Target[i][j].g] += 1;
            histo_T_B[Target[i][j].b] += 1;
        }
    }


    //get P
    double P_S_R[256] = {0, };
    double P_S_G[256] = {0, };
    double P_S_B[256] = {0, };

    double P_T_R[256] = {0, };
    double P_T_G[256] = {0, };
    double P_T_B[256] = {0, };


    for(unsigned int t=0; t<256; t++){

        P_S_R[t] = (double)histo_S_R[t]/(double)Source.Size();
        P_S_G[t] = (double)histo_S_G[t]/(double)Source.Size();
        P_S_B[t] = (double)histo_S_B[t]/(double)Source.Size();

        P_T_R[t] = (double)histo_T_R[t]/(double)Target.Size();
        P_T_G[t] = (double)histo_T_G[t]/(double)Target.Size();
        P_T_B[t] = (double)histo_T_B[t]/(double)Target.Size();

        //qDebug() << t << " : " << P_S_R[t];

    }


    //get y , yp
    double y_R[256] = {0, };
    double y_G[256] = {0, };
    double y_B[256] = {0, };

    double yp_R[256] = {0, };
    double yp_G[256] = {0, };
    double yp_B[256] = {0, };

    y_R[0] = P_S_R[0];
    y_G[0] = P_S_G[0];
    y_B[0] = P_S_B[0];

    yp_R[0] = P_T_R[0];
    yp_G[0] = P_T_G[0];
    yp_B[0] = P_T_B[0];

    for(unsigned int r=1; r<256; r++){
        y_R[r] = y_R[r-1] + P_S_R[r];
        //r_R[r] = T_R[r] * 255;
        y_G[r] = y_G[r-1] + P_S_G[r];
        //r_G[r] = T_G[r] * 255;
        y_B[r] = y_B[r-1] + P_S_B[r];
        //r_B[r] = T_B[r] * 255;

        yp_R[r] = yp_R[r-1] + P_T_R[r];
        //r_R[r] = T_R[r] * 255;
        yp_G[r] = yp_G[r-1] + P_T_G[r];
        //r_G[r] = T_G[r] * 255;
        yp_B[r] = yp_B[r-1] + P_T_B[r];
        //r_B[r] = T_B[r] * 255;

        //qDebug() << r << " : " << yp_R[r];

    }


    int tr_R[256] = {0, };
    int tr_G[256] = {0, };
    int tr_B[256] = {0, };

    for(unsigned int i=0; i<256; i++){
         double min_R = 100000.0, min_G = 100000.0, min_B = 100000.0;
        for(unsigned int j = 0; j<256; j++)
        {
            if(min_R>std::fabs(y_R[i]-yp_R[j]))
            {
                min_R = std::fabs(y_R[i]-yp_R[j]);
                tr_R[i] = j;
            }

            if(min_G>std::fabs(y_G[i]-yp_G[j]))
            {
                min_G = std::fabs(y_G[i]-yp_G[j]);
                tr_G[i] = j;
            }

            if(min_B>std::fabs(y_B[i]-yp_B[j]))
            {
                min_B = std::fabs(y_B[i]-yp_B[j]);
                tr_B[i] = j;
            }
        }
        //qDebug() << i << " : " << min_R;
    }


    for(unsigned int i=0; i<256; i++)
    {
        qDebug() << i << " : " << tr_R[i];
    }


    for(unsigned int i=0; i<Source.Row(); i++){
        for(unsigned int j=0; j<Source.Col(); j++)
        {
            Source[i][j].r=tr_R[Source[i][j].r];
            Source[i][j].g=tr_G[Source[i][j].g];
            Source[i][j].b=tr_B[Source[i][j].b];
        }
    }



    ImageForm*  q_pForm2 = new ImageForm(Source, "Histogram Matching", this);
    _plpImageForm->Add(q_pForm2);
    q_pForm2->show();

    //ImageForm*  q_pForm1 = new ImageForm(Target, "target", this);
    //_plpImageForm->Add(q_pForm1);
    //q_pForm1->show();

}

void MainFrame::on_buttonHistoTarget_clicked()
{
    KImageColor Target;

    //포커스 된 ImageForm으로부터 영상을 가져옴
    if(_q_pFormFocused != 0 && _q_pFormFocused->ImageColor().Address() &&  _q_pFormFocused->ID() == "OPEN")
    {
        Target = _q_pFormFocused->ImageColor();
    }
    else
        return;

    ImageForm*  q_pForm2 = new ImageForm(Target, "target", this);
    _plpImageForm->Add(q_pForm2);
    q_pForm2->show();
}

void MainFrame::on_button_BoxFilter_clicked()
{
    //포커스 된 ImageForm으로부터 영상을 가져옴
    KImageColor Gaussian_Noise;
    KImageColor Salt_Noise;

    //포커스 된 ImageForm으로부터 영상을 가져옴
    if(_q_pFormFocused != 0 && _q_pFormFocused->ImageColor().Address() &&  _q_pFormFocused->ID() == "OPEN")
    {
        Gaussian_Noise = _q_pFormFocused->ImageColor();
        Salt_Noise = _q_pFormFocused->ImageColor();
    }
    else
        return;




    //Gaussian Noise 생성
    double GauP = ui->spin_GaussianParameter->value();

    KGaussian GauN;
    GauN.Create(0, GauP*GauP);
    GauN.OnRandom(Gaussian_Noise.Size());

    double dNoise = 0;
    double Kp = 8;
    for(int i = 0; i<Gaussian_Noise.Row();i++){
        for(int j = 0; j<Gaussian_Noise.Col();j++){
            dNoise = GauN.Generate();

            if(Gaussian_Noise[i][j].r + dNoise*Kp > 255)
                Gaussian_Noise[i][j].r = 255;
            else if(Gaussian_Noise[i][j].r +dNoise*Kp < 0)
                Gaussian_Noise[i][j].r = 0;
            else
                Gaussian_Noise[i][j].r += dNoise*Kp;

            if(Gaussian_Noise[i][j].g + dNoise*Kp > 255)
                Gaussian_Noise[i][j].g = 255;
            else if(Gaussian_Noise[i][j].g +dNoise*Kp < 0)
                Gaussian_Noise[i][j].g = 0;
            else
                Gaussian_Noise[i][j].g += dNoise*Kp;

            if(Gaussian_Noise[i][j].b + dNoise*Kp > 255)
                Gaussian_Noise[i][j].b = 255;
            else if(Gaussian_Noise[i][j].b +dNoise*Kp < 0)
                Gaussian_Noise[i][j].b = 0;
            else
                Gaussian_Noise[i][j].b += dNoise*Kp;
        }
    }


    //Salt Noise 생성
    srand((int)time(NULL));


    double Thres = 0.005;
    double minus_Thres = (double)(1-Thres);
    double random;


    for(int i=0; i<Salt_Noise.Row(); i++){
        for(int j=0;j<Salt_Noise.Col();j++){

            random = (double) rand()/RAND_MAX;

            if(random<Thres)
            {
                Salt_Noise[i][j].r = 0;
                Salt_Noise[i][j].g = 0;
                Salt_Noise[i][j].b = 0;
            }
            else if(random>minus_Thres){
                Salt_Noise[i][j].r = 255;
                Salt_Noise[i][j].g = 255;
                Salt_Noise[i][j].b = 255;
            }


        }
    }



    //show noise
    ImageForm*  q_pForm1 = new ImageForm(Salt_Noise, "Salt_and_Pepper Noise", this);
    _plpImageForm->Add(q_pForm1);
    q_pForm1->show();

    ImageForm*  q_pForm2 = new ImageForm(Gaussian_Noise, "Gaussian Noise", this);
    _plpImageForm->Add(q_pForm2);
    q_pForm2->show();





    //필터링
    int filtersize = ui->spin_FilterSize->value();

    //Box Filtering
        //S_n_P noise 이미지 불러오기
    ImageForm* q_pForm_s = 0;

        for(int i=0; i<_plpImageForm->Count();i++)
            if((*_plpImageForm)[i]->ID()=="Salt_and_Pepper Noise")
            {
                q_pForm_s = (*_plpImageForm)[i];
                break;
            }

    KImageColor Salt_Noise_2 = q_pForm_s->ImageColor();

    //Gaussian noise 이미지 불러오기

    ImageForm* q_pForm_g = 0;

        for(int i=0; i<_plpImageForm->Count();i++)
            if((*_plpImageForm)[i]->ID()=="Gaussian Noise")
            {
                q_pForm_g = (*_plpImageForm)[i];
                break;
            }
    KImageColor Gaussian_Noise_2 = q_pForm_g->ImageColor();

    int size = 0.5*filtersize-0.5;
    int box = filtersize*filtersize;
    //qDebug()<<size;
        //Salt Noise
    double sum_r_s,sum_g_s,sum_b_s;
    double sum_r_g,sum_g_g,sum_b_g;

    for(int ii = size; ii<Salt_Noise_2.Row()-size;ii++){
        for(int jj = size; jj<Salt_Noise_2.Col()-size;jj++){
            sum_r_s = 0;
            sum_g_s = 0;
            sum_b_s = 0;

            sum_r_g = 0;
            sum_g_g = 0;
            sum_b_g = 0;
            for(int i = -size; i<size+1; i++){
                for(int j = -size; j<size+1;j++){
                    sum_r_s+=Salt_Noise[ii-i][jj-j].r;
                    sum_g_s+=Salt_Noise[ii-i][jj-j].g;
                    sum_b_s+=Salt_Noise[ii-i][jj-j].b;

                    sum_r_g+=Gaussian_Noise[ii-i][jj-j].r;
                    sum_g_g+=Gaussian_Noise[ii-i][jj-j].g;
                    sum_b_g+=Gaussian_Noise[ii-i][jj-j].b;
                }
            }

            //qDebug()<<sum_r/(filtersize*filtersize);

            Salt_Noise_2[ii][jj].r=sum_r_s/box;
            Salt_Noise_2[ii][jj].g=sum_g_s/box;
            Salt_Noise_2[ii][jj].b=sum_b_s/box;

            Gaussian_Noise_2[ii][jj].r=sum_r_g/box;
            Gaussian_Noise_2[ii][jj].g=sum_g_g/box;
            Gaussian_Noise_2[ii][jj].b=sum_b_g/box;
        }
    }


    //show noise 제거
    ImageForm*  q_pForm3 = new ImageForm(Salt_Noise_2, "Box_Salt_n_Pepper", this);
    _plpImageForm->Add(q_pForm3);
    q_pForm3->show();

    ImageForm*  q_pForm4 = new ImageForm(Gaussian_Noise_2, "Box_Gaussian", this);
    _plpImageForm->Add(q_pForm4);
    q_pForm4->show();

}

void MainFrame::on_button_GaussianFilter_clicked()
{
    //포커스 된 ImageForm으로부터 영상을 가져옴
    KImageColor Gaussian_Noise;
    KImageColor Salt_Noise;

    //포커스 된 ImageForm으로부터 영상을 가져옴
    if(_q_pFormFocused != 0 && _q_pFormFocused->ImageColor().Address() &&  _q_pFormFocused->ID() == "OPEN")
    {
        Gaussian_Noise = _q_pFormFocused->ImageColor();
        Salt_Noise = _q_pFormFocused->ImageColor();
    }
    else
        return;




    //Gaussian Noise 생성
    double sigma = ui->spin_GaussianParameter->value();

    KGaussian GauN;
    GauN.Create(0, sigma*sigma);
    GauN.OnRandom(Gaussian_Noise.Size());

    double dNoise = 0;
    double Kp = 8;
    for(int i = 0; i<Gaussian_Noise.Row();i++){
        for(int j = 0; j<Gaussian_Noise.Col();j++){
            dNoise = GauN.Generate();

            if(Gaussian_Noise[i][j].r + dNoise*Kp > 255)
                Gaussian_Noise[i][j].r = 255;
            else if(Gaussian_Noise[i][j].r +dNoise*Kp < 0)
                Gaussian_Noise[i][j].r = 0;
            else
                Gaussian_Noise[i][j].r += dNoise*Kp;

            if(Gaussian_Noise[i][j].g + dNoise*Kp > 255)
                Gaussian_Noise[i][j].g = 255;
            else if(Gaussian_Noise[i][j].g +dNoise*Kp < 0)
                Gaussian_Noise[i][j].g = 0;
            else
                Gaussian_Noise[i][j].g += dNoise*Kp;

            if(Gaussian_Noise[i][j].b + dNoise*Kp > 255)
                Gaussian_Noise[i][j].b = 255;
            else if(Gaussian_Noise[i][j].b +dNoise*Kp < 0)
                Gaussian_Noise[i][j].b = 0;
            else
                Gaussian_Noise[i][j].b += dNoise*Kp;
        }
    }


    //Salt Noise 생성
    srand((int)time(NULL));


    double Thres = 0.005;
    double minus_Thres = (double)(1-Thres);
    double random;


    for(int i=0; i<Salt_Noise.Row(); i++){
        for(int j=0;j<Salt_Noise.Col();j++){

            random = (double) rand()/RAND_MAX;

            if(random<Thres)
            {
                Salt_Noise[i][j].r = 0;
                Salt_Noise[i][j].g = 0;
                Salt_Noise[i][j].b = 0;
            }
            else if(random>minus_Thres){
                Salt_Noise[i][j].r = 255;
                Salt_Noise[i][j].g = 255;
                Salt_Noise[i][j].b = 255;
            }


        }
    }



    //show noise
    ImageForm*  q_pForm1 = new ImageForm(Salt_Noise, "Salt_and_Pepper Noise", this);
    _plpImageForm->Add(q_pForm1);
    q_pForm1->show();

    ImageForm*  q_pForm2 = new ImageForm(Gaussian_Noise, "Gaussian Noise", this);
    _plpImageForm->Add(q_pForm2);
    q_pForm2->show();





    //Gaussian 필터링

    //S_n_P noise 이미지 불러오기
    ImageForm* q_pForm_s = 0;

    for(int i=0; i<_plpImageForm->Count();i++)
        if((*_plpImageForm)[i]->ID()=="Salt_and_Pepper Noise")
        {
            q_pForm_s = (*_plpImageForm)[i];
            break;
        }

    KImageColor Salt_Noise_2 = q_pForm_s->ImageColor();

    //Gaussian noise 이미지 불러오기

    ImageForm* q_pForm_g = 0;

    for(int i=0; i<_plpImageForm->Count();i++)
        if((*_plpImageForm)[i]->ID()=="Gaussian Noise")
        {
            q_pForm_g = (*_plpImageForm)[i];
            break;
        }
    KImageColor Gaussian_Noise_2 = q_pForm_g->ImageColor();

    int Gau_filter_size = 8*sigma + 1; //filter size




    int size = std::sqrt(Gau_filter_size);


    //qDebug()<<size;
        //Salt Noise

    double mask_r_s, mask_g_s, mask_b_s;
    double mask_r_g, mask_g_g, mask_b_g;


    for(int ii = size; ii<Salt_Noise_2.Row()-size;ii++){
        for(int jj = size; jj<Salt_Noise_2.Col()-size;jj++){

            mask_r_s = 0;
            mask_g_s = 0;
            mask_b_s = 0;

            mask_r_g = 0;
            mask_g_g = 0;
            mask_b_g = 0;


            for(int i = -size; i<size+1; i++){
                for(int j = -size; j<size+1;j++){
                    mask_r_s += std::exp(-0.5*((i*i+j*j)/(sigma*sigma)))*Salt_Noise[ii-i][jj-j].r;
                    mask_g_s += std::exp(-0.5*((i*i+j*j)/(sigma*sigma)))*Salt_Noise[ii-i][jj-j].g;
                    mask_b_s += std::exp(-0.5*((i*i+j*j)/(sigma*sigma)))*Salt_Noise[ii-i][jj-j].b;

                    mask_r_g += std::exp(-0.5*((i*i+j*j)/(sigma*sigma)))*Gaussian_Noise[ii-i][jj-j].r;
                    mask_g_g += std::exp(-0.5*((i*i+j*j)/(sigma*sigma)))*Gaussian_Noise[ii-i][jj-j].g;
                    mask_b_g += std::exp(-0.5*((i*i+j*j)/(sigma*sigma)))*Gaussian_Noise[ii-i][jj-j].b;
                }
            }

            Salt_Noise_2[ii][jj].r = (1/(2*M_PI*sigma*sigma))*mask_r_s;
            Salt_Noise_2[ii][jj].g = (1/(2*M_PI*sigma*sigma))*mask_g_s;
            Salt_Noise_2[ii][jj].b = (1/(2*M_PI*sigma*sigma))*mask_b_s;

            Gaussian_Noise_2[ii][jj].r = (1/(2*M_PI*sigma*sigma))*mask_r_g;
            Gaussian_Noise_2[ii][jj].g = (1/(2*M_PI*sigma*sigma))*mask_g_g;
            Gaussian_Noise_2[ii][jj].b = (1/(2*M_PI*sigma*sigma))*mask_b_g;



        }
    }


    //show noise 제거
    ImageForm*  q_pForm3 = new ImageForm(Salt_Noise_2, "Gaussian_Salt_n_Pepper", this);
    _plpImageForm->Add(q_pForm3);
    q_pForm3->show();

    ImageForm*  q_pForm4 = new ImageForm(Gaussian_Noise_2, "Gaussian_Gaussian", this);
    _plpImageForm->Add(q_pForm4);
    q_pForm4->show();
}

void MainFrame::on_button_MedianFilter_clicked()
{
    //포커스 된 ImageForm으로부터 영상을 가져옴
    KImageColor Gaussian_Noise;
    KImageColor Salt_Noise;

    //포커스 된 ImageForm으로부터 영상을 가져옴
    if(_q_pFormFocused != 0 && _q_pFormFocused->ImageColor().Address() &&  _q_pFormFocused->ID() == "OPEN")
    {
        Gaussian_Noise = _q_pFormFocused->ImageColor();
        Salt_Noise = _q_pFormFocused->ImageColor();
    }
    else
        return;




    //Gaussian Noise 생성
    double GauP = ui->spin_GaussianParameter->value();

    KGaussian GauN;
    GauN.Create(0, GauP*GauP);
    GauN.OnRandom(Gaussian_Noise.Size());

    double dNoise = 0;
    double Kp = 8;
    for(int i = 0; i<Gaussian_Noise.Row();i++){
        for(int j = 0; j<Gaussian_Noise.Col();j++){
            dNoise = GauN.Generate();

            if(Gaussian_Noise[i][j].r + dNoise*Kp > 255)
                Gaussian_Noise[i][j].r = 255;
            else if(Gaussian_Noise[i][j].r +dNoise*Kp < 0)
                Gaussian_Noise[i][j].r = 0;
            else
                Gaussian_Noise[i][j].r += dNoise*Kp;

            if(Gaussian_Noise[i][j].g + dNoise*Kp > 255)
                Gaussian_Noise[i][j].g = 255;
            else if(Gaussian_Noise[i][j].g +dNoise*Kp < 0)
                Gaussian_Noise[i][j].g = 0;
            else
                Gaussian_Noise[i][j].g += dNoise*Kp;

            if(Gaussian_Noise[i][j].b + dNoise*Kp > 255)
                Gaussian_Noise[i][j].b = 255;
            else if(Gaussian_Noise[i][j].b +dNoise*Kp < 0)
                Gaussian_Noise[i][j].b = 0;
            else
                Gaussian_Noise[i][j].b += dNoise*Kp;
        }
    }


    //Salt Noise 생성
    srand((int)time(NULL));


    double Thres = 0.005;
    double minus_Thres = (double)(1-Thres);
    double random;


    for(int i=0; i<Salt_Noise.Row(); i++){
        for(int j=0;j<Salt_Noise.Col();j++){

            random = (double) rand()/RAND_MAX;

            if(random<Thres)
            {
                Salt_Noise[i][j].r = 0;
                Salt_Noise[i][j].g = 0;
                Salt_Noise[i][j].b = 0;
            }
            else if(random>minus_Thres){
                Salt_Noise[i][j].r = 255;
                Salt_Noise[i][j].g = 255;
                Salt_Noise[i][j].b = 255;
            }


        }
    }



    //show noise
    ImageForm*  q_pForm1 = new ImageForm(Salt_Noise, "Salt_and_Pepper Noise", this);
    _plpImageForm->Add(q_pForm1);
    q_pForm1->show();

    ImageForm*  q_pForm2 = new ImageForm(Gaussian_Noise, "Gaussian Noise", this);
    _plpImageForm->Add(q_pForm2);
    q_pForm2->show();





    //필터링
    int filtersize = ui->spin_FilterSize->value();

    //Median Filtering
        //S_n_P noise 이미지 불러오기
    ImageForm* q_pForm_s = 0;

        for(int i=0; i<_plpImageForm->Count();i++)
            if((*_plpImageForm)[i]->ID()=="Salt_and_Pepper Noise")
            {
                q_pForm_s = (*_plpImageForm)[i];
                break;
            }

    KImageColor Salt_Noise_2 = q_pForm_s->ImageColor();

    //Gaussian noise 이미지 불러오기

    ImageForm* q_pForm_g = 0;

        for(int i=0; i<_plpImageForm->Count();i++)
            if((*_plpImageForm)[i]->ID()=="Gaussian Noise")
            {
                q_pForm_g = (*_plpImageForm)[i];
                break;
            }
    KImageColor Gaussian_Noise_2 = q_pForm_g->ImageColor();


    int size = 0.5*filtersize-0.5;
    std::vector<double> v_r_s,v_g_s,v_b_s;
    std::vector<double> v_r_g,v_g_g,v_b_g;
    int box = filtersize*filtersize;

    //qDebug()<<size;
        //Salt Noise


    for(int ii = size; ii<Salt_Noise_2.Row()-size;ii++){
        for(int jj = size; jj<Salt_Noise_2.Col()-size;jj++){

            v_r_s.clear();
            v_g_s.clear();
            v_b_s.clear();

            v_r_g.clear();
            v_g_g.clear();
            v_b_g.clear();


            for(int i = -size; i<size+1; i++){
                for(int j = -size; j<size+1;j++){
                    v_r_s.push_back(Salt_Noise[ii-i][jj-j].r);
                    v_g_s.push_back(Salt_Noise[ii-i][jj-j].g);
                    v_b_s.push_back(Salt_Noise[ii-i][jj-j].b);

                    v_r_g.push_back(Gaussian_Noise[ii-i][jj-j].r);
                    v_g_g.push_back(Gaussian_Noise[ii-i][jj-j].g);
                    v_b_g.push_back(Gaussian_Noise[ii-i][jj-j].b);
                }
            }

            //qDebug()<<sum_r/(filtersize*filtersize);
            sort(v_r_s.begin(),v_r_s.end());
            sort(v_g_s.begin(),v_g_s.end());
            sort(v_b_s.begin(),v_b_s.end());

            sort(v_r_g.begin(),v_r_g.end());
            sort(v_g_g.begin(),v_g_g.end());
            sort(v_b_g.begin(),v_b_g.end());

            Salt_Noise_2[ii][jj].r = v_r_s[box/2];
            Salt_Noise_2[ii][jj].g = v_g_s[box/2];
            Salt_Noise_2[ii][jj].b = v_b_s[box/2];

            Gaussian_Noise_2[ii][jj].r = v_r_g[box/2];
            Gaussian_Noise_2[ii][jj].g = v_g_g[box/2];
            Gaussian_Noise_2[ii][jj].b = v_b_g[box/2];


        }
    }


    //show noise 제거
    ImageForm*  q_pForm3 = new ImageForm(Salt_Noise_2, "Median_Salt_n_Pepper", this);
    _plpImageForm->Add(q_pForm3);
    q_pForm3->show();

    ImageForm*  q_pForm4 = new ImageForm(Gaussian_Noise_2, "Median_Gaussian", this);
    _plpImageForm->Add(q_pForm4);
    q_pForm4->show();

}

void MainFrame::on_button_CannyEdge_clicked()
{

    KImageGray icMain;
    KImageGray icMain2;

    //포커스 된 ImageForm으로부터 영상을 가져옴
    if(_q_pFormFocused != 0 && _q_pFormFocused->ImageGray().Address() &&  _q_pFormFocused->ID() == "OPEN")
    {
        icMain = _q_pFormFocused->ImageGray();
        icMain2 = _q_pFormFocused->ImageGray();
    }
    else
        return;


    //Canny edge
    int row = icMain.Row();
    int col = icMain.Col();

    //Gaussian Filter
    double sigma = ui->SpinBox_sigma->value();

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


    qDebug()<<(unsigned char)((((int)(113/22.5)+1)>>1) & 0x00000003);

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


                qDebug()<<ii<<" "<<jj<<" "<<phase<<" "<<direction[ii][jj];


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
            icMain2[ii][jj] = realedge[ii][jj];
        }
    }



    //Show
    ImageForm*  q_pForm1 = new ImageForm(icMain2, "Canny Edge Operator", this);
    _plpImageForm->Add(q_pForm1);
    q_pForm1->show();

}

typedef struct loc{
    int x;
    int y;
}Location;

void MainFrame::on_button_GeneralHough_clicked()
{
    KImageGray icMain;
    //포커스 된 ImageForm으로부터 영상을 가져옴
    if(_q_pFormFocused != 0 && _q_pFormFocused->ImageGray().Address() &&  _q_pFormFocused->ID() == "OPEN")
    {
        icMain = _q_pFormFocused->ImageGray();
    }
    else
        return;

    KImageGray icMain2(icMain.Row(),icMain.Col());

    int f_row[86] = {0,};
    int f_col[86] = {0,};

    //파일 읽기
    std::ifstream readFile;
    readFile.open("C:/Qt/plug.txt");
    int fx,fy,i=0;

    while(!readFile.eof()){
        readFile>>fx>>fy;
        f_col[i] = fy;
        f_row[i] = fx;
        i++;
    }




    //table
    int Xc=0, Yc=0; //중심점
    for(int i=0; i<86;i++)
    {
        Xc += f_row[i];
        Yc += f_col[i];
    }
    Xc /= 86;
    Yc /= 86;


    typedef struct r{
        double Ri,Ai;
    }R;

    std::vector<R> table[4];
    R t;
    int Xi,Yi;

        //Grad 구하기
    double dTmp, dGradX, dGradY;
    int dir;

    for(int j = 1,jj = 86-2; jj; j++, jj--)
    {
        dGradX = (float)(f_row[j+1] - f_row[j-1]) + 1e-8;
        dGradY = (float)(f_col[j+1] - f_col[j-1]) + 1e-8;

        dTmp = (180.0/M_PI)*(atan2(dGradY,dGradX)) + 90;
        dir = ((((int)(dTmp/22.5) + 1) >>1) & 0x00000003);

        //qDebug() << f_row[j] << " " << f_col[j] << " " <<dir;

        t.Ri = sqrt((Xc - f_row[j])*(Xc - f_row[j]) + (Yc - f_col[j])*(Yc - f_col[j]));
        t.Ai = (180.0/M_PI) * atan2((double)(Yc - f_col[j]),(double)(Xc - f_row[j]));

        table[dir].push_back(t);

    }


    //canny edge
    double sigma = ui->SpinBox_sigma->value();

    Canny_Edge(icMain,sigma);

    std::vector<Location> edge;
    Location location;

    for(int i=0;i<icMain.Row();i++)
    {
        for(int j=0;j<icMain.Col();j++)
        {
            if(icMain[i][j] == 255)
            {
                location.x = i;
                location.y = j;
                edge.push_back(location);
                //qDebug() << "x : " << i << "y : " << j;
            }
        }
    }




    //Hough Transform
    int** voting = new int*[icMain.Row()]{0,};
    for(int i = 0; i<icMain.Row(); i++)
    {
        voting[i] = new int[icMain.Col()]{0,};
    }


    /*
    for(int i=0; i < 4; i++)
    {
        qDebug() << i << " " << table[i].size();
    }
    //qDebug() << edge.size();
    */

        //voting
    for(int i = 0; i < edge.size(); i++)
    {
        Xi = edge[i].x;
        Yi = edge[i].y;

            //qDebug() << "ok";
        for(int j=0;j<4;j++)
        {

           for(int k=0;k<table[j].size();k++)
           {

                Xc = Xi - 0.85*table[j][k].Ri*cos(table[j][k].Ai*M_PI/180);
                Yc = Yi - 0.85*table[j][k].Ri*sin(table[j][k].Ai*M_PI/180);

                //qDebug() << "i : " << i << "j : " << j << "k : " << k << "x : " << Xc << "y : " << Yc;


                if(Xc>0 && Xc<icMain.Row() && Yc>0 && Yc<icMain.Col())
                {
                    voting[Xc][Yc]++;
                    if(icMain2[Xc][Yc]+50 < 256)
                        icMain2[Xc][Yc] += 50;
                    else icMain2[Xc][Yc]  = 255;
                }

           }

        }

    }



        //max voting -> 원의 중심점 찾기
    int max = 0;

    for(int ii = 0; ii<icMain.Row();ii++){
        for(int jj = 0; jj<icMain.Col();jj++){
            if(voting[ii][jj]>max){
                Xc = ii;
                Yc = jj;
                max = voting[ii][jj];
            }
        }
    }

    for(int ii = 0; ii<icMain.Row();ii++){
        for(int jj = 0; jj<icMain.Col();jj++){
            icMain[ii][jj] = 0;
        }
    }


        //find
    for(int j = 0; j < 4; j++)
    {
        for(int k = 0; k<table[j].size();k++)
        {
            Xi = Xc + 0.85*table[j][k].Ri*cos(table[j][k].Ai*M_PI/180);
            Yi = Yc + 0.85*table[j][k].Ri*sin(table[j][k].Ai*M_PI/180);

            if(Xi>0 && Xi<icMain.Row() && Yi>0 && Yi<icMain.Col())
            {
                icMain[Xi][Yi] = 255;
            }
        }
    }





    //Show
    ImageForm*  q_pForm1 = new ImageForm(icMain, "General Hough Transform", this);
    _plpImageForm->Add(q_pForm1);
    q_pForm1->show();

    ImageForm*  q_pForm2 = new ImageForm(icMain2, "General Hough Transform_voting", this);
    _plpImageForm->Add(q_pForm2);
    q_pForm2->show();

}



void MainFrame::on_button_CircleHough_clicked()
{
    KImageGray icMain;

    //포커스 된 ImageForm으로부터 영상을 가져옴
    if(_q_pFormFocused != 0 && _q_pFormFocused->ImageGray().Address() &&  _q_pFormFocused->ID() == "OPEN")
    {
        icMain = _q_pFormFocused->ImageGray();
    }
    else
        return;

    KImageGray icMain2(icMain.Row(),icMain.Col());

    //canny edge operator
    double sigma = ui->SpinBox_sigma->value();

    Canny_Edge(icMain,sigma);

    //Hough Transform
    int** voting = new int*[icMain.Row()]{0,};
    for(int i = 0; i<icMain.Row(); i++)
    {
        voting[i] = new int[icMain.Col()]{0,};
    }





    std::vector<Location> edge;
    Location location;



        //Edge 좌표 받아오기
    for(int ii = 0; ii<icMain.Row();ii++){
        for(int jj = 0; jj<icMain.Col();jj++){

            if(icMain[ii][jj] == 255)
            {
                location.x = ii;
                location.y = jj;
                edge.push_back(location);

                 //qDebug() << "x : " << ii << "y : " << jj;
            }

        }
    }


        //voting
    int x,y;
    int a,b;
    for(int i=0;i<edge.size();i++)
    {
        x = edge[i].x;
        y = edge[i].y;

        //qDebug() << "x : " << x << "y : " << y;

        for(int j=0;j<360;j++){
            a = x-51.5*cos(j*M_PI/180);
            b = y-51.5*sin(j*M_PI/180);


            if(a>0 && a<icMain.Row() && b>0 && b<icMain.Col())
            {
                voting[a][b]++;
                if(icMain2[a][b]+30 < 256)
                    icMain2[a][b] += 30;
                else icMain2[a][b]  = 255;
            }


        }


    }



        //max voting -> 원의 중심점 찾기
    int max = 0;

    for(int ii = 0; ii<icMain.Row();ii++){
        for(int jj = 0; jj<icMain.Col();jj++){
            if(voting[ii][jj]>max){
                a = ii;
                b = jj;
                max = voting[ii][jj];
            }
        }
    }

    for(int ii = 0; ii<icMain.Row();ii++){
        for(int jj = 0; jj<icMain.Col();jj++){
            icMain[ii][jj] = 0;
            //icMain2[ii][jj] = 0;
        }
    }


        //circle
    for(int j=0;j<360;j++)
    {
        x = a-51.5*cos(j*M_PI/180);
        y = b-51.5*sin(j*M_PI/180);
        icMain[x][y] = 255;
    }





    //Show
    ImageForm*  q_pForm1 = new ImageForm(icMain, "Circle Hough Transform", this);
    _plpImageForm->Add(q_pForm1);
    q_pForm1->show();

    ImageForm*  q_pForm2 = new ImageForm(icMain2, "Circle Hough Transform_voting", this);
    _plpImageForm->Add(q_pForm2);
    q_pForm2->show();


}

void MainFrame::on_button_GSS_DOG_clicked()
{
    KImageGray icMain,Origin,Half1,Half2;



    //포커스 된 ImageForm으로부터 영상을 가져옴
    if(_q_pFormFocused != 0 && _q_pFormFocused->ImageGray().Address() &&  _q_pFormFocused->ID() == "OPEN")
    {
        icMain = _q_pFormFocused->ImageGray();
    }

    else
        return;

    KImageColor icMain2(icMain.GrayToRGB()), icMain3(icMain.GrayToRGB());

    //Scale space 생성
    double sigma = ui->Spinsigma->value();

    KImageDouble igMain(icMain);
    Origin = igMain.ToGray();
    Half1 = igMain.HalfSize().ToGray();
    Half2 = igMain.HalfSize().HalfSize().ToGray();

    std::vector<std::vector<KImageDouble>> Octave(3); //Octave1,2,3


    for(int j=0;j<5;j++)
    {
        //qDebug() <<"1";
        Octave[0].push_back(Gaussian_Filter(Origin,sigma*pow(1.41,(j-1))));
        //qDebug() <<"2";
        Octave[1].push_back(Gaussian_Filter(Half1,sigma*pow(1.41,(j-1))));
        //qDebug() <<"3";
        Octave[2].push_back(Gaussian_Filter(Half2,sigma*pow(1.41,(j-1))));

    }



    //Dog
    std::vector<std::vector<KImageDouble>> Dog(3);
    for(int i = 0; i < 3;i++)
    {
        KImageDouble Dog_img(Octave[i][0].Row(),Octave[i][0].Col());

        for(int octave=1;octave<5;octave++)
        {

            for(int ii=0;ii<Octave[i][0].Row();ii++)
            {
                for(int jj=0;jj<Octave[i][0].Col();jj++)
                {
                    Dog_img[ii][jj] = Octave[i][octave][ii][jj]-Octave[i][octave-1][ii][jj];
                }
            }

            Dog[i].push_back(Dog_img);
            //qDebug() <<Octave[i][0].Row();
        }
    }



    //KeyPoint & filtering
    std::vector<std::vector<real_key>> KeyPoint = Find_Key_Point(Dog); //루트2 시그마, 2시그마

    qDebug() << KeyPoint[0].size();
    qDebug() << KeyPoint[1].size();

    int x, y;


    Orientation(KeyPoint,Octave[0],sigma); //orientation
    KeyPoint_descriptor(KeyPoint,Octave[0],sigma);


    for(int scale = 0; scale<KeyPoint.size(); scale++) //root 2 sigma, 2 sigma
    {
        for(int i=0;i<KeyPoint[scale].size();i++)
        {
            Mark_KeyPoint(icMain, KeyPoint[scale][i].magnitude, KeyPoint[scale][i].direction, KeyPoint[scale][i].x, KeyPoint[scale][i].y);
        }
    }

    std::ofstream writeFile("sift.txt");


    for(int i = 0; i<KeyPoint.size();i++)
    {
        for(int j=0;j<KeyPoint[i].size();j++)
        {
            //KeyPoint[i][j].x " " KeyPoint[i][j].y " " KeyPoint[i][j].direction " " KeyPoint[i][j].magnitude
            //"\n"
            writeFile << KeyPoint[i][j].x << " " << KeyPoint[i][j].y << " " << KeyPoint[i][j].direction << " " << KeyPoint[i][j].magnitude
                      << "\n";


            for(int f = 0;f<4;f++)
            {
                for(int r = 0; r<KeyPoint[i][j].feature[f].size();r++)
                //KeyPoint[i][j].feature[f][r] " "
                    writeFile << KeyPoint[i][j].feature[f][r] << " ";
                writeFile << "\n";

            }
            //"\n"
            writeFile << "\n";
        }
      }




    //show

    ImageForm*  q_pForm1 = new ImageForm(icMain, "KeyPoint", this);
    _plpImageForm->Add(q_pForm1);
    q_pForm1->show();



    ImageForm*  q_pForm2 = new ImageForm(Show(Octave), "Scale-Space", this);
    _plpImageForm->Add(q_pForm2);
    q_pForm2->show();

    ImageForm*  q_pForm3 = new ImageForm(Show(Dog), "Dog", this);
    _plpImageForm->Add(q_pForm3);
    q_pForm3->show();

    ImageForm*  q_pForm4 = new ImageForm(icMain2, "Origin", this);
    _plpImageForm->Add(q_pForm4);
    q_pForm4->show();


}




void MainFrame::on_button_SiftMatching_clicked()
{

}

void MainFrame::on_button_opticalflow_clicked()
{
    int StartImgnum = ui->spinBox_StartImg->value();
    int EndImgnum = ui->spinBox_EndImg->value();

    std::vector<KImageGray> Imgvec;

    QString q_fileName;

    for(int i = StartImgnum; i <= EndImgnum; i++)
    {
        if(i<10)
        {
            q_fileName = QString::fromStdString("./data/yos.0" + std::to_string(i) + ".pgm");
        }
        else q_fileName = QString::fromStdString("./data/yos." + std::to_string(i) + ".pgm");

        ImageForm*  q_pForm = new ImageForm(q_fileName, "Open", this);
        _plpImageForm->Add(q_pForm);

        KImageGray* PGM = &q_pForm->ImageGray();//.GaussianSmoothed(2);
        Imgvec.emplace_back(*PGM);
    }

    qDebug()<<"1";

    int row = Imgvec[0].Row();
    int col = Imgvec[0].Col();

    UV** mat_uv = new UV*[row]{0,};

    for(int i = 0; i < row; i++)
    {
        mat_uv[i] = new UV[col]{{0,0},};
    }


    //Optical_Flow(Imgvec[0],Imgvec[2],mat_uv);
    //Draw(Imgvec[0],mat_uv);


    for(int i = 0; i < Imgvec.size()-2; i++)
    {
        Optical_Flow(Imgvec[i],Imgvec[i+2],mat_uv);
        Draw(Imgvec[i],mat_uv);

        ImageForm*  q_pForm2 = new ImageForm(Imgvec[i],"Optical Flow", this);
        _plpImageForm->Add(q_pForm2);
        q_pForm2->show();
    }



    /*
    for(int i=0;i<Imgvec.size();i++)
    {
        ImageForm*  q_pForm2 = new ImageForm(Imgvec[i], "Match", this);
        _plpImageForm->Add(q_pForm2);
        q_pForm2->show();
    }
    */











}
