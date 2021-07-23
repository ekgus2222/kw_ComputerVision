#ifndef MAINFRAME_H
#define MAINFRAME_H

#include <QDialog>
#include "kfc.h"

namespace Ui {
class MainFrame;

}

class ImageForm;
class KVoronoiDgm;
class KPotentialField;

class MainFrame : public QDialog
{
    Q_OBJECT

private:
    Ui::MainFrame *ui;

    KPtrList<ImageForm*>*   _plpImageForm;
    ImageForm*              _q_pFormFocused;

public:
    explicit MainFrame(QWidget *parent = 0);
    ~MainFrame();

    void            ImageFormFocused(ImageForm* q_pImageForm)
                    {   _q_pFormFocused  = q_pImageForm;   //활성화된 창의 포인터를 저장함
                        UpdateUI();                        //UI 활성화 갱신
                    }
    void            UpdateUI();
    void            CloseImageForm(ImageForm* pForm);

public:
    void            OnMousePos(const int& nX, const int& nY, ImageForm* q_pForm);

private slots:
    void on_buttonOpen_clicked();
    void on_buttonDeleteContents_clicked();    
    void on_buttonSepiaTone_clicked();
    void on_buttonShowList_clicked();

    void on_spinSaturation_valueChanged(const QString &arg1);

    void on_spinHue_valueChanged(const QString &arg1);

    void on_spinHue_textChanged(const QString &arg1);

    void on_spinSaturation_textChanged(const QString &arg1);

    void on_buttonOtsu_clicked();

    void on_button_dil_n_ero_clicked();

    void on_button3x3_clicked();

    void on_button5x5_clicked();

    void on_buttonHEQ_clicked();

    void on_buttonHMA_clicked();

    void on_buttonHistoTarget_clicked();

    void on_button_BoxFilter_clicked();

    void on_button_GaussianFilter_clicked();

    void on_button_MedianFilter_clicked();

    void on_button_CannyEdge_clicked();

    void on_button_GeneralHough_clicked();

    void on_button_CircleHough_clicked();

    void on_button_GSS_DOG_clicked();

    void on_button_SiftMatching_clicked();

    void on_button_opticalflow_clicked();

protected:
    void closeEvent(QCloseEvent* event);
};

#endif // MAINFRAME_H
