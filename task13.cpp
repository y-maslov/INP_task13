#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <fstream>
#include <iostream>
#include <TMath.h>
#include <TF1.h>
#include <TF1Convolution.h>

Double_t E_min = 3; // GeV
Double_t E_max = 3.2; // GeV
Double_t E_step = 0.5 / 1000.0; // GeV
Double_t Mass = 3.096; // GeV
Double_t Width = 92.9e-6; // GeV

Double_t dopel_gaus(Double_t *x, Double_t *a)
{
    return (TMath::Gaus(x[0], 0, a[0], true) + TMath::Gaus(x[0], a[1], a[2], true));
}

Double_t BW_plus_linear(Double_t *x, Double_t *a)
{
    return a[0] * TMath::BreitWigner(x[0], Mass, Width) + (a[1] + a[2] * x[0]);
}


void task13()
{
    Int_t nbins = (E_max - E_min) / E_step;
    

    auto hist = new TH1D("Hist", "Invariant mass of 3#pi, GeV/c^{2}; m_{3#pi}", nbins, E_min, E_max);

    Double_t input = 0;
    std::ifstream file("m3piJPSI_cut.dat");
    while (file >> input) hist->Fill(input);
  
    TF1 *BW_plus_linear_func = new TF1("BW_plus_linear_func", BW_plus_linear, E_min, E_max, 3);
    TF1 *dopel_gaus_func = new TF1("dopel_gaus_func", dopel_gaus, E_min, E_max, 3);
    
    TF1Convolution *Conv = new TF1Convolution(BW_plus_linear_func, dopel_gaus_func, E_min, E_max, true);
    Conv->SetRange(0, E_max + E_max);
    Conv->SetNofPointsFFT(5000);

    TF1Convolution *Conv1 = new TF1Convolution(BW_plus_linear_func, dopel_gaus_func, E_min, E_max, true);
    Conv1->SetRange(0, E_max + E_max);
    Conv1->SetNofPointsFFT(1000);

    TF1Convolution *Conv2 = new TF1Convolution(BW_plus_linear_func, dopel_gaus_func, E_min, E_max, true);
    Conv2->SetRange(0, E_max + E_max);
    Conv2->SetNofPointsFFT(3000);

    TF1 *fit_func = new TF1("fit", *Conv, E_min, E_max, 6);
    TF1 *fit_func1 = new TF1("fit1", *Conv1, E_min, E_max, 6);
    TF1 *fit_func2 = new TF1("fit2", *Conv2, E_min, E_max, 6);
    fit_func->SetLineColor(kRed);
    fit_func1->SetLineColor(kBlue);
    fit_func2->SetLineColor(kGreen);

    Double_t f_param[6] = {1.83170e+01, 1.00000e+03, -1.00000e+03, -1.60913e+00, 4.34670e-03, 1.17265e-02};
    Double_t f1_param[6] = {1, 1, -1, -1, 4, 1};
    Double_t f2_param[6] = {1.83170e+01, 1.00000e+03, -1.00000e+03, -1.60913e+00, 4.34670e-03, 1.17265e-02};
    fit_func->SetParameters(f_param);
    fit_func1->SetParameters(f1_param);
    fit_func2->SetParameters(f2_param);

    fit_func->SetParLimits(0, -1000, 1000);
    fit_func->SetParLimits(1, -1000, 100000);
    fit_func->SetParLimits(2, -1000, 100000);
    fit_func->SetParLimits(3, -1000, 1000);
    fit_func->SetParLimits(4, -1000, 1000);
    fit_func->SetParLimits(5, -1000, 1000);

    fit_func1->SetParLimits(0, -1000, 1000);
    fit_func1->SetParLimits(1, -1000, 1000);
    fit_func1->SetParLimits(2, -1000, 1000);
    fit_func1->SetParLimits(3, -1000, 1000);
    fit_func1->SetParLimits(4, -1000, 1000);
    fit_func1->SetParLimits(5, -1000, 1000);

    fit_func2->SetParLimits(0, -1000, 1000);
    fit_func2->SetParLimits(1, -1000, 1000);
    fit_func2->SetParLimits(2, -1000, 1000);
    fit_func2->SetParLimits(3, -1000, 1000);
    fit_func2->SetParLimits(4, -1000, 1000);
    fit_func2->SetParLimits(5, -1000, 1000);


    auto canvas = new TCanvas("Canvas", "Canvas", 800, 600);
    auto legend = new TLegend();
    legend->AddEntry(fit_func, "5000 points", "l"); 
    legend->AddEntry(fit_func1, "1000 points", "l"); 
    legend->AddEntry(fit_func2, "3000 points", "l"); 
    auto fit = hist->Fit("fit");
    auto fit1 = hist->Fit("fit1", "+");
    auto fit2= hist->Fit("fit2", "+");
    canvas->SetLogy();
    hist->Draw();
    legend->Draw();


 
}