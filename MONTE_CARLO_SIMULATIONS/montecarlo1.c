#include <math.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>

#include "TCanvas.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TAxis.h"
//____________________________________________________________________________________________________________
// TOY MONTECARLO

// Calcolo funzione inversa Erf (normalizzata tra -1 e 1)

float myErfInv2(float x){
    float tt1, tt2, lnx, sgn, PI = M_PI;
    sgn = (x < 0) ? -1.0f : 1.0f;
    
    x = (1 - x)*(1 + x);        // x = 1 - x*x;
    lnx = logf(x);
    
    tt1 = 2/(PI*0.147) + 0.5f * lnx;
    tt2 = 1/(0.147) * lnx;
    
    return(sgn*sqrtf(-tt1 + sqrtf(tt1*tt1 - tt2)));
}

// Calcolo funzione inversa esponenziale (normalizzata tra 1300 e 1600)

double myExpInv(double x, double c){
    double Na, E;
    
    Na = c/(pow(-M_E,(-c*1600))+pow(M_E,(-c*1300)));
    
    // Na=17.5438 costante di normalizzazione
    
    
    E = (-1.0/c)*log((-c*(x-1))/Na);
    
    
    return E;
}
//____________________________________________________________________________________________________________


#define N 100000    // Servirà per generare l'Istogramma

void montecarlo1 () {
    FILE *fd;       //apertura file
    
    double r, k, med, sigma, pi=M_PI;
    int i;
    srand(time(NULL));
    
    // printf("inserisci media: \n");
    // scanf("%lf",&med);
    // printf("inserisci deviazione: \n");
    // scanf("%lf",&sigma);
    
//__________________________________________________________________________________________________________
    
    // PLOT con ROOT
    TCanvas* canv1 = new TCanvas("canv1","My canvas1", 800, 600);  // 800, 600 rappresentano i pixel
    canv1 -> Divide(1,1);       // Divide lo schermo in una matrice 1x1
    canv1->cd(1);               // Inserisce il plot nella posizione (1,1) della matrice
    
    TCanvas* canv2 = new TCanvas("canv2","My canvas2", 800, 600);
    canv2 -> Divide(1,1);
    canv2->cd(1);
    
    TCanvas* canv3 = new TCanvas("canv3","My canvas3", 800, 600);
    canv3 -> Divide(1,1);
    canv3->cd(1);
    
    TCanvas* canv4 = new TCanvas("canv4","My canvas4", 800, 600);
    canv4 -> Divide(1,1);
    canv4->cd(1);
    
    TCanvas* canv5 = new TCanvas("canv5","My canvas5", 800, 600);
    canv5 -> Divide(1,1);
    canv5->cd(1);
    
    TCanvas* canv6 = new TCanvas("canv6","My canvas6", 800, 600);
    canv6 -> Divide(1,1);
    canv6->cd(1);
    
    TCanvas* canv7 = new TCanvas("canv7","My canvas7", 800, 600);
    canv7 -> Divide(1,1);
    canv7->cd(1);
    
    TCanvas* canv8 = new TCanvas("canv8","My canvas8", 800, 600);
    canv8 -> Divide(1,1);
    canv8->cd(1);
    
    TCanvas* canv9 = new TCanvas("canv9","My canvas9", 800, 600);
    canv9 -> Divide(1,1);
    canv9->cd(1);
    
    TCanvas* canv10 = new TCanvas("canv10","My canvas10", 800, 600);
    canv10 -> Divide(1,1);
    canv10->cd(1);
    
    
    // Titolo del file
    
    canv1->SetTitle("Eta_X");
    canv2->SetTitle("Theta_X");
    canv3->SetTitle("pt_X");
    canv4->SetTitle("pboos_X");
    canv5->SetTitle("ptot_X");
    canv6->SetTitle("Eta_Z");
    canv7->SetTitle("Theta_Z");
    canv8->SetTitle("pt_Z");
    canv9->SetTitle("pboos_Z");
    canv10->SetTitle("ptot_Z");
    
  
    // Particella X
    
    TH1D* histo_eta_X = new TH1D("hist_eta_X","Eta_X distribution; ; # counts", 100, -8, 8);
    TH1D* histo_theta_X = new TH1D("hist_theta_X","Theta_X distribution; [rad]; # counts; ", 100, -0.4, 4.5);
    TH1D* histo_pt_X = new TH1D("hist_pt_X","pt_X distribution; [MeV/c]; # counts", 100, -10, 170);
    TH1D* histo_pboos_X = new TH1D("hist_pboos_X","pboos_X distribution; [MeV/c]; # counts", 100, -1500, 1000);
    TH1D* histo_ptot_X = new TH1D("hist_tot_X","ptot_X distribution; [MeV/c]; # counts", 100, -300, 1450);
    
    
    // Particella Z
    
    TH1D* histo_eta_Z = new TH1D("hist_eta_Z","Eta_Z distribution; ; # counts", 100, -8, 8);
    TH1D* histo_theta_Z = new TH1D("hist_theta_Z","Theta_Z distribution; [rad]; # counts", 100, -0.4, 4.5);
    TH1D* histo_pt_Z = new TH1D("hist_pt_Z","pt_Z distribution; [MeV/c]; # counts", 100, -10, 180);
    TH1D* histo_pboos_Z = new TH1D("hist_pboos_Z","pboos_Z distribution; [MeV/c]; # counts", 100, -300, 300);
    TH1D* histo_ptot_Z = new TH1D("hist_tot_Z","ptot_Z distribution; [MeV/c]; # counts", 100, -50, 300);
    
//_____________________________________________________________________________________________________________
    
    
    med = 1400; // Valore medio dell'energia nel SCM in GeV (√s = 13 TeV)
    sigma = 92; // Sigma dell'energia nel SCM in GeV
    
    // Definisco le variabili del sistema
    
    double M, m1, m2, phi, theta, p[3][N], p_SCM, c;
    // M è l'energia del sistema, m1 la massa della X,
    // m2 la massa della Z, p_SCM l'impulso del SCM, c lo shape dell'exp
    
 
    
    // Richiede quale valore attribuire a c
    
    printf("inserisci c: \n");
    scanf("%lf",&c);
    
    // Massa particella Z (in GeV)
    m1 = 91;
    
    // Massa particella X (in GeV)
    m2 = 1200;
    
//________________________________________________________________________________________________________
    
    // Definisco impulso nel sistema del LAB
    
    double p_lab, med_lab = 180.0, sigma_lab = 450.0;
    
//________________________________________________________________________________________________________
    
    // Creazione effettiva dell'istogramma
    
    for (i = 0; i<N; i++) {
        
        // Genera le energie (M) secondo distribuzione esponenziale
        
        r = (rand() /(RAND_MAX+1.0));
        k = myExpInv(r,c);
        M = k;
        
        // if(M>1291){
        
        // Generazione impulso nel sistema di laboratorio (SLAB)
        
        p_lab = 2*(rand() /(RAND_MAX+1.0))-1;
        p_lab = myErfInv2(p_lab);
        p_lab = med_lab + (sigma_lab*p_lab);
        
        
        double dpz1, dpz2, Xi1, Xi2;
        
        dpz1 = ((2*6500-p_lab)-sqrt(pow((2*6500-p_lab),2)-4*4*pow((M),2)))/2; // Delta pz del protone 1
        dpz2 = pow(M,2)/(4*dpz1);                                             // Delta pz del protone 2
        Xi1 = dpz1/6500;                                                      // dpz1/[impulso protone 1]
        Xi2 = dpz2/6500;                                                      // dpz2/[impulso protone 2]
        
        
        p_SCM = pow(((pow(M,2)-pow(m1,2)+pow(m2,2))/(2*M)),2)-pow(m2,2);  // (Impulso sistema centro massa)^2
        p_SCM = sqrt(p_SCM);                                              // Impulso sistema centro massa (SCM)
        
        theta = (pi)*(rand() /(RAND_MAX+1.0));             // Angolo con asse z
        phi = 2*pi*(rand() /(RAND_MAX+1.0));               // Angolo perpendicolare con asse z
        
        // Impulso del SCM in coordinate cartesiane
        
        p[0][i] = p_SCM*sin(theta)*cos(phi);     // Proiezione impulso lungo asse x
        p[1][i] = p_SCM*sin(theta)*sin(phi);     // Proiezione impulso lungo asse y
        p[2][i] = p_SCM*cos(theta);              // Proiezione impulso lungo asse z
        
        
        // Definisco l'impulso trasverso pt, l'impulso bustato pboos,
        //l'impulso totale p_tot della particella X
        
        double pt_X, pboos_X, ptot_X, thetalab_X, eta_X, gamma, Beta;
        
        gamma = 1.0/(sqrt(1.0-pow((p_lab/med),2)));
        Beta = p_lab/med;
        
        
        // Particella X
        
        pt_X = sqrt(pow(p[0][i],2)+pow(p[1][i],2));
        pboos_X = (gamma)*(p[2][i]-((p_lab/M)*sqrt(pow(m2,2)+pow(p_SCM,2))));
        ptot_X = sqrt(pow(pboos_X,2)+pow(pt_X,2));
        thetalab_X = acos(pboos_X/ptot_X);
        eta_X = -log(tan(thetalab_X/2));
        
        // printf("  pt_X = %.3f\n, pboos_X = %.3f\n, ptot_X = %.3f\n, thetalab_X = %.3f\n, eta_X = %.3f\n", pt_X, pboos_X, ptot_X, thetalab_X, eta_X);
        
        // Definisco l'impulso trasverso pt, l'impulso bustato pboos,
        //l'impulso totale ptot della particella Z
        
        double pt_Z, pboos_Z, ptot_Z, thetalab_Z, eta_Z;
        
        
        // Particella Z
        
        pt_Z = sqrt(pow(-p[0][i],2)+pow(-p[1][i],2));
        pboos_Z = (gamma)*(-p[2][i]-((p_lab/M)*sqrt(pow(m1,2)+pow(p_SCM,2))));
        ptot_Z = sqrt(pow(pboos_Z,2)+pow(pt_Z,2));
        thetalab_Z = acos(pboos_Z/ptot_Z);
        eta_Z = -log(tan(thetalab_Z/2));
        
        // printf("  pt_Z = %.3f\n, pboos_Z\n = %.3f, ptot_Z\n = %.3f, thetalab_Z\n = %.3f, eta_Z = %.3f\n", pt_Z, pboos_Z, ptot_Z, thetalab_Z, eta_Z);
        
//_________________________________________________________________________________________________________
        
        // PLOT con ROOT particella X
        
         histo_eta_X->Fill(eta_X);
         histo_theta_X->Fill(thetalab_X);
         histo_pt_X->Fill(pt_X);
         histo_pboos_X->Fill(pboos_X);
         histo_ptot_X->Fill(ptot_X);
        
        // PLOT con ROOT particella Z
        
         histo_eta_Z->Fill(eta_Z);
         histo_theta_Z->Fill(thetalab_Z);
         histo_pt_Z->Fill(pt_Z);
         histo_pboos_Z->Fill(pboos_Z);
         histo_ptot_Z->Fill(ptot_Z);
        
//____________________________________________________________________________________________________________
      
  //      if (Xi2>=0.12) i--;      // Taglio sulla Xi2
    }
    
    fclose(fd);
    //__________________________________________________________________________________________________________
    
    
    // Particella X
    
    canv1-> cd(1);
    histo_eta_X->Draw();
    
    canv2-> cd(1);
    histo_theta_X->Draw();
    
    canv3-> cd(1);
    histo_pt_X->Draw();
    
    canv4-> cd(1);
    histo_pboos_X->Draw();
    
    canv5-> cd(1);
    histo_ptot_X->Draw();
    
    // Particella Z
    
    canv6-> cd(1);
    histo_eta_Z->Draw();
    
    canv7-> cd(1);
    histo_theta_Z->Draw();
    
    canv8-> cd(1);
    histo_pt_Z->Draw();
    
    canv9-> cd(1);
    histo_pboos_Z->Draw();
    
    canv10-> cd(1);
    histo_ptot_Z->Draw();
    
    // Salvataggio plots in pdf
 
    canv1->SaveAs("Eta_X.pdf");
    canv2->SaveAs("Theta_X.pdf");
    canv3->SaveAs("pt_X.pdf");
    canv4->SaveAs("pboos_X.pdf");
    canv5->SaveAs("ptot_X.pdf");
    canv6->SaveAs("Eta_Z.pdf");
    canv7->SaveAs("Theta_Z.pdf");
    canv8->SaveAs("pt_Z.pdf");
    canv9->SaveAs("pboos_Z.pdf");
    canv10->SaveAs("ptot_Z.pdf");
  
    
}


