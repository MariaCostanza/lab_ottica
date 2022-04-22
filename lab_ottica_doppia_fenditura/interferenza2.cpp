#include <iostream>
#include "TMath.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TFile.h"
#include "TMatrixD.h"
#include "TFitResult.h"
#include "TROOT.h"

/*
par0=i0
par1=posizione centro
*/

Double_t myFunction(Double_t *x, Double_t *par){
Double_t xx= TMath::Abs(x[0]- par[1]);
Double_t d = 0.15E-3;
Double_t D = 0.50E-3;
Double_t L = 0.9;
Double_t k = 2*TMath::Pi()/632.8E-9;
Double_t sint = xx/sqrt(xx*xx+L*L);
Double_t f = sin(k*d*sint*0.5)/(k*d*sint*0.5);
Double_t val= par[0]*pow(cos(0.5*k*D*sint),2)*pow(f,2);
return val;
}

void interferenza() {

gStyle->SetOptStat(2210);
gStyle->SetOptFit(1111);
gStyle->SetFitFormat("7.6g");


TGraph *graph = new TGraph("int2_1.txt", "%lg %lg");

TCanvas *c = new TCanvas ("Young","Interferenza doppia fenditura 1");


Double_t x, y;
    graph->GetPoint(0, x, y);
    Double_t max_x = x, max_y = y;
    for(int i = 1; i < graph->GetN(); i++)
    { graph->GetPoint(i, x, y);
        if(y > max_y) {
           max_x = x;
           max_y = y;
        }
    }
    std::cout<<max_x<<"     "<<max_y<<'\n';
    
TF1 *f = new TF1("f", myFunction, 0.06, 0.075, 2);
f->SetParameter(0,max_y);
f->SetParLimits(0,max_y-1,max_y+1);
f->SetParameter(1,max_x);
f->SetParLimits(1,max_x-0.001,max_x+0.001);
f->SetLineColor(kRed);

//-----cosmetics-----
graph->SetMarkerStyle(7);
graph->SetMarkerColor(kBlue);
graph->GetYaxis()->SetTitleOffset(1.2);
graph->GetXaxis()->SetTitleSize(0.04);
graph->GetYaxis()->SetTitleSize(0.04);
graph->GetXaxis()->SetTitle("Posizione (m)");
graph->GetYaxis()->SetTitle("Intensita' luminosa (u arbitraria)");
graph->SetTitle("Doppia fenditura 2");


graph->GetXaxis()->SetRangeUser(0.06, 0.075);
graph->Fit("f","R");
graph->Draw("AP");


Double_t chiSquare=f->GetChisquare();
Double_t nDOF=f->GetNDF();
Double_t Prob=f->GetProb();
std::cout << "ChiSquare = " <<chiSquare << " , nDOF " << nDOF <<endl;
std::cout << "ChiSquare ridotto = " <<chiSquare/nDOF<<endl;
std::cout << "ChiSquare Probability= " <<Prob <<endl; 


c->Print("int2_1.jpg");
c->Print("int2_1.root");

}