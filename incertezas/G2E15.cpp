#include <iostream>
using namespace std;

/*
void histo_y_funcion(){
  
  gROOT->Reset();
  gRandom->GetSeed();
  
  int n_bines = 20;
  int corridas = 1000;
  float rango_min = 5;
  float rango_max = 15;
  float sigma = 1;
  float mu = 3;
  float r = 0;
  
  //crea la ventana, le define el color de fondo y hace que lo que se dibuje en ella tenga grilla
  TCanvas *c1 = new TCanvas("c1","Experiencia (B)", 200, 10, 1000, 600);
  c1->SetFillColor(42);
  c1->SetGrid();

  //crea el objeto TF1 con la funcion a graficar. Pone nombre a los tres parametros y fifa sus valores
  TF1 *funcion_gauss = new TF1("funcion_gauss","[0]/(sqrt(2*3.14159)*[1])*TMath::Exp(-0.5*pow((x-[2])/[1],2))",rango_min,rango_max);
  funcion_gauss->SetParNames("normalizacion","sigma","media");
  funcion_gauss->SetParameter(0,corridas*(rango_max-rango_min)/n_bines);
  funcion_gauss->SetParameter(1,0.75);
  funcion_gauss->SetParameter(2,0);
  
  //dibuja el histograma y la funcion
  funcion_gauss->Draw("same");
}
*/

int binomial(int n, float p){ //Obtiene cantidad de exitos en una experiencia binomial
	gRandom->GetSeed();
	int A = 0;
	float r = 0;
	int i = 0;
	for (i = 0; i < n; i++){
		r = gRandom->Uniform(1);
		if (r < p){
			A += 1;
		}
	}
	return A;
}

void G2_E15(){
	int i = 0;
	int j = 0;
	int corridas = 1000;
	float I = 15;
	int m = 1000;
	int A = 0
	float pEm = (float)I/m;
	float e = 0.75;
	
	int op = 1;
	
	int n_bins = 10;
	float min = 0;
	float max = 30;
	
	//crea la ventana, le define el color de fondo y hace que lo que se dibuje en ella tenga grilla
	TCanvas *c1 = new TCanvas("c1","Experiencia (C)", 200, 10, 1000, 600);
	c1->SetFillColor(42);
	c1->SetGrid();
	
	//crea el histograma y define el color de relleno
	TH1F *histo  = new TH1F("histo","Experimento", n_bins, min, max);
	histo->SetFillColor(45);
	
	//llena el histograma con valores del la experiencia
	for(int j=0; j < corridas; j++){
		switch (op){
			case 1:
				histo->Fill(binomial(n,e));
				break;
			case 2:
				for(int i = 0; i < m; i++) {
					A += binomial(1,pEm);
				}
				histo->Fill(A);
				break;
			case 3:
				for(int i = 0; i < m; i++) {
					A += binomial(1, p) * binomial(1,e);
				}
				histo->Fill(A);
				break;
			default:
				break;
		}
	}	
	
	//Dibuja el histograma y la funcion
	histo->Draw();
}
