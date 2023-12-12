#include <TH1D.h>
#include <TCanvas.h>
#include <fstream>
#include <iostream>
#include <TMath.h>
#include <TF1.h>
#include <TF1Convolution.h>

Double_t E_min = 3; // GeV
Double_t E_max = 3.2; // GeV
Double_t E_step = 0.5 / 1000.0; // GeV
Double_t Mass = 1.3; // GeV
Double_t Width = 0.4; // GeV

Double_t dopel_gaus(Double_t *x, Double_t *a)
{
    return (TMath::Gaus(x[0], 0, a[0]) + TMath::Gaus(x[0], a[1], a[2]));
}

Double_t BW_plus_linear(Double_t *x, Double_t *a)
{
    return TMath::BreitWigner(x[0], Mass, Width) + a[1] + a[2] * x[0];
}

void task13()
{
    Int_t nbins = (E_max - E_min) / E_step;
    auto canvas = new TCanvas("Canvas", "Canvas", 800, 600);
    auto hist = new TH1D("Hist", "Hist", nbins, E_min, E_max);

    Double_t input = 0;
    std::ifstream file("m3piJPSI_cut.dat");
    while (file >> input) hist->Fill(input);
  

    TF1 *dopel_gaus_func = new TF1("dopel_gaus_func", dopel_gaus, E_min, E_max, 3);
    TF1 *BW_plus_linear_func = new TF1("BW_plus_linear_func", BW_plus_linear, E_min, E_max, 3);
    TF1Convolution *Conv = new TF1Convolution(BW_plus_linear_func,dopel_gaus_func, E_min, E_max, true);
    Conv->SetNofPointsFFT(10000);

    TF1 *fit_func = new TF1("fit", *Conv, E_min, E_max, 6);
    fit_func->SetParameter(0, 1);
    fit_func->SetParameter(1, 1);
    fit_func->SetParameter(2, 1);
    fit_func->SetParameter(3, 1);
    fit_func->SetParameter(4, 1);
    fit_func->SetParameter(5, 1);



    fit_func->SetParLimits(0, -1000, 10000);
    fit_func->SetParLimits(1, -1000, 10000);
    fit_func->SetParLimits(2, 0.1, 10000);
    fit_func->SetParLimits(3, -1000, 10000);
    fit_func->SetParLimits(4, 0.1, 10000);
    fit_func->SetParLimits(5, 0.1, 10000);
    auto fit = hist->Fit("fit");

    hist->Draw();
}