#include "lar_utils.C"

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TObjString.h"

#include <CLHEP/Random/RandBinomial.h>
#include <CLHEP/Random/JamesRandom.h>
#include <CLHEP/Vector/ThreeVector.h>

#include <random>
#include <vector>
#include <iostream>
#include <fstream>

void main(){
	
	TString filename = "./data/muon_correlated.root";
	int PDG =  13;
	
	TFile *file = TFile::Open(filename);
	if(!file || file->IsZombie()) {
		std::cerr << "Error al abrir el archivo." << std::endl;
		return;
	}
	
	// Extraer "particula_correlated"
    	TObjArray* parts = filename.Tokenize("/");
    	TString file_with_extension = ((TObjString*)parts->Last())->GetString();
    	TString file_without_extension = file_with_extension.ReplaceAll(".root", "");
	// Extraer "particula"
    	TObjArray* subparts = file_without_extension.Tokenize("_");
    	TString base_name = ((TObjString*)subparts->First())->GetString();
	
	
	TTree *tree = (TTree*)file->Get("lightcalo/LightCaloSBND");
	if(!tree) {
		std::cerr << "Error al obtener el TTree." << std::endl;
        	file->Close();
        	return;
	}
	
	Long64_t nEntries = tree->GetEntries();
	
	std::vector<double> *EnDep_v = new std::vector<double>();
	std::vector<double> *EnDepStep_v = new std::vector<double>();
	std::vector<double> *EnDepPDG_v = new std::vector<double>();
	std::vector<double> *EnDepX_v = new std::vector<double>();
	std::vector<double> *EnDepY_v = new std::vector<double>();
	std::vector<double> *EnDepZ_v = new std::vector<double>();
	std::vector<double> *NScintPhotons_v = new std::vector<double>();
	std::vector<double> *NDriftElectrons_v = new std::vector<double>();
		 
	tree->SetBranchAddress("EnDep_v", &EnDep_v);
	tree->SetBranchAddress("EnDepStep_v", &EnDepStep_v);
	tree->SetBranchAddress("EnDepPDG_v", &EnDepPDG_v);
	tree->SetBranchAddress("EnDepX_v", &EnDepX_v);
	tree->SetBranchAddress("EnDepY_v", &EnDepY_v);
	tree->SetBranchAddress("EnDepZ_v", &EnDepZ_v);
	tree->SetBranchAddress("NScintPhotons_v", &NScintPhotons_v);
	tree->SetBranchAddress("NDriftElectrons_v", &NDriftElectrons_v);
    	
    	CLHEP::HepJamesRandom engine;
	CLHEP::RandBinomial fBinomialGen(engine);
    	
	int nbins = 1000;
	
	TH1F *hist1 = new TH1F("hist1", " ", nbins, 0, 1000);
	TH1F *hist2 = new TH1F("hist2", " ", nbins, 0, 5e4);
	TH1F *hist3 = new TH1F("hist3", " ", nbins, 0, 5e5);
	
	std::vector<double> dEdx_v;
	std::vector<double> num_photons_v;
    	
    	for(Long64_t i=0; i < nEntries; i++){
    		tree->GetEntry(i);
    		Long64_t SIZE = EnDep_v->size();
    		if(i%1000==0)
    		{
    			std::cout << "########################################" << std::endl;
    			std::cout << "Leyendo entrada del TTree: "<< i << std::endl;
    			std::cout << "########################################" << std::endl;
    		}
    		for(Long64_t j = 0; j < SIZE-1; j++){		
    			if(EnDepPDG_v->at(j) == PDG){
		    		double energy_deposit = EnDep_v->at(j); // MeV
    				double energy_step = EnDepStep_v->at(j); // cm
    				
    				double num_ions = 0.;
    				if (energy_deposit >= fWion) num_ions = energy_deposit / fWion;
    				
    				double num_quanta = 0; 	 
    				if (energy_deposit >= fWph) num_quanta = energy_deposit / fWph;
    				
    				double dEdx = (energy_step <= 0.0) ? 0.0 : energy_deposit / energy_step;
    				dEdx = (dEdx < 1.) ? 1. : dEdx;
				
    			    	
    			    	double x0 = EnDepX_v->at(j);
    				double y0 = EnDepY_v->at(j);
    				double z0 = EnDepZ_v->at(j);			
    			    	
    			    	double EfieldStep = EfieldAtStep(x0, y0, z0);
				double R = 0.;
				double num_electrons = 0.;
    			    	
    				if (EfieldStep > 0.){
	    				
	    				double x1 = EnDepX_v->at(j+1);
	    				double y1 = EnDepY_v->at(j+1);
	    				double z1 = EnDepZ_v->at(j+1);
	    				
	    				if(fUseModBoxRecomb){
						R = fModBoxRecomb(EfieldStep, dEdx);
		    			}

		    			else if(fUseEllipsModBoxRecomb){
		 				double phi = AngleToEFieldStep(EfieldStep, x0, y0, z0, x1, y1, z1);
		 				if (std::isnan(phi)) {
		  					R = fModBoxRecomb(EfieldStep, dEdx);
						}
						else {
							double B_ellips =
		    					fEllipsModBoxB * dEdx / (EfieldStep * std::hypot(std::sin(phi), std::cos(phi) / fEllipsModBoxR));
		  					R = std::log(fEllipsModBoxA + B_ellips) / B_ellips;
						}
	    				}
		    			
		    			else{
		    				R = fRecombA / (1. + dEdx * (fRecombk) / EfieldStep);
		    			}
		    		}
    					
    				if(LArQL && PDG != 1000020040){
	    				R += fEscapingEFraction(dEdx) * fFieldCorrection(EfieldStep, dEdx);
	    			}
	    				
	    			if(R < 0.){
					R = 0.;
				} else if(R > 1.){
					R = 1.;
				}
					
				if (num_ions > 0.) {
      					num_electrons = (fUseBinomialFlucts) ? fBinomialGen.fire(num_ions, R) : (num_ions * R);
				}

				double num_photons = (num_quanta - num_electrons);
				
				
				if (PDG == 2212){
					num_electrons = num_electrons*fQProton;
					num_photons = num_photons*fQProton;
				}
				
				if (PDG == 1000020040) {
     					num_electrons = num_electrons * fQAlpha;
      					num_photons = num_photons * fQAlpha;
    				}
				
								
				hist1->Fill(dEdx);
				hist2->Fill(num_photons/energy_deposit);
				hist3->Fill(NScintPhotons_v->at(j)/fScintPreScale/energy_deposit);

				dEdx_v.push_back(dEdx);
				num_photons_v.push_back(num_photons/energy_deposit);
    			}
    		}
    	}
	
	double num_photons_E5_med = calcularMedia(num_photons_v);
	double dEdx_E5_med = calcularMedia(dEdx_v);
	
	TCanvas *c1 = new TCanvas();
	hist1->Draw();
	hist1->GetXaxis()->SetTitle("dE/dx (MeV/cm)");
	hist2->SetFillColor(kBlue-9);
	
	TCanvas *c2 = new TCanvas();
	c2->SetLogy();
	c2->SetLeftMargin(0.15);
	c2->SetRightMargin(0.15);
		
	hist2->GetYaxis()->SetTitle("Eventos");
	hist2->GetXaxis()->SetTitle("L (fotones/MeV)");
	hist2->SetFillColor(kBlue-9); 
    	hist2->SetLineColor(kBlue);
	
    	hist3->SetLineColor(kRed);
	
	hist2->Draw();
	hist2->Rebin(5);
	hist3->Draw("SAME");
	hist3->Rebin(5);
	
	hist2->SetStats(0);
	hist3->SetStats(0);
	
	TLegend *leg1 = new TLegend(0.7, 0.7, 0.9, 0.9);
	leg1->SetBorderSize(0);
	leg1->AddEntry(hist2, "Toy MC", "f");
	leg1->AddEntry(hist3, "LArSoft", "l");
	leg1->Draw();
	
	TString imageName2 = "./images/" + file_without_extension + "_NScintPhotonsR.pdf";
	c2->Print(imageName2);
	
	std::cout << " " << std::endl;
	std::cout << "MEDIAS DE CADA HISTOGRAMA: " << base_name << std::endl;
	std::cout << "ToyModel : L (MeV^{-1}) Media = " << hist2->GetMean() << std::endl;
	std::cout << "LArSoft : L (MeV^{-1}) Media = " << hist3->GetMean() << std::endl;
	
	std::cout << " " << std::endl;
	
	std::cout << "Deposición de energía media 0.5 kV/cm: " << dEdx_E5_med  << " MeV/(g/cm^2)" <<std::endl;
	std::cout << "Yield medio a campo eléctrico 0.5 kV/cm: " << num_photons_E5_med << "(MeV^{-1})" << std::endl;
	std::cout << "Yield relativo (0.5 kV/cm): " << num_photons_E5_med/(1/fWph) << std::endl;

	std::cout << " " << std::endl;
	
}

int main() {
    main();
    return 0;
}
