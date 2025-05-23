#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string.h>
#include <numeric>
#include "ConfigFile.tpp"
#include <algorithm>

using namespace std;


double energy(const std::vector<double>& fnow, double dx) {
  double ener = 0.0; // TODO: compute quantity E 
  for (int i = 1; i < fnow.size(); ++i) {
    ener += 0.5 * dx * (fnow[i]*fnow[i] + fnow[i-1]*fnow[i-1]); // regle du trapeze ici 
  }
  return ener;
}


void boundary_condition(vector<double> &fnext, vector<double> &fnow, double const& A, double om, \
		double const& t,double const& dt, \
		vector<double> &beta2, string &bc_l, string &bc_r, int &N)
{
      if (bc_l == "fixe"){
        fnext[0] = fnow[0]; 
        //fnext[0]  = 0.0;
	// NB: on peut aussi utiliser la condition "excitation" et poser A=0
      }else if(bc_l == "libre"){
        fnext[0] = fnext[1]; // TODO : Modifier pour imposer la condition au bord gauche libre
      }else if (bc_l =="sortie"){
        fnext[0] = fnow[0] + (sqrt(beta2[0]))*(fnow[1]-fnow[0]);
        //fnext[0] = fnow[1]- (1.0 /sqrt(beta2[1]))*(fnow[1]-fnow[0]); // TODO : Modifier pour imposer la condition au bord gauche "sortie de l'onde"
      }else if (bc_l == "excitation"){
        fnext[0] = A * sin(om*t); // TODO : Modifier pour imposer la condition au bord gauche sinusoidale
      }else{
        cerr << "Merci de choisir une condition aux bord gauche valide" << endl;
      }
	      
      if (bc_r == "fixe"){
        fnext[N-1] = fnow[N-1]; 
        //fnext[N-1] = 0.0;
	// NB: on peut aussi utiliser la condition "excitation" et poser A=0	
      }else if(bc_r == "libre"){
        fnext[N-1] = fnext[N-2]; 
      }else if (bc_r =="sortie"){
        fnext[N-1] = fnow[N-1] - sqrt(beta2[N-1])*(fnow[N-1]- fnow[N-2]); 
      }else if (bc_r == "excitation"){ 
        fnext[N-1] = A * sin(om*t); 
      }else{
        cerr << "Merci de choisir une condition aux bord droit valide" << endl;
      }
}



double finit(double x, double n_init, double L, double f_hat, double x1, double x2, string initialization)
{
  const double PI = 3.14159265358979323846;
  double finit_ = 0.0;

  // Empêche les valeurs hors du domaine [0, L]
  if (x < 0.0 || x > L) {
      return 0.0;
  }

  if (initialization == "mode") {
      // Mode propre : cosinus avec mode n_init
      finit_ = cos((n_init + 0.5) * PI * x / L);
  } else {
      // Forme donnée par l'équation (3)
      if (x <= x1) {
          finit_ = 0.0;
      } else if (x > x1 && x < x2) {
          finit_ = 0.5 * f_hat * (1.0 - cos(2.0 * PI * (x - x1) / (x2 - x1)));
      } else { // x >= x2
          finit_ = 0.0;
      }
  }

  return finit_;
}

//
// Surcharge de l'operateur pour ecrire les elements d'un tableau
//
template <class T> ostream& operator<< (ostream& o, vector<T> const& v)
{
  unsigned int len(v.size());
  for(unsigned int i(0); i < (len - 1); ++i)
    o << v[i] << " ";
  if(len > 0)
    o << v[len-1];
  return o;
}

//
// Main
//
int main(int argc, char* argv[])
{
  const double PI = 3.1415926535897932384626433832795028841971e0;
  const double g  = 9.81;
  double dx;
  double dt;
  double t;
  double Nsteps;
  int stride(0);
  

  string inputPath("configuration.in"); // Fichier d'input par defaut
  if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice7 config_perso.in")
    inputPath = argv[1];

  ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

  for(int i(2); i<argc; ++i) // Input complementaires ("./Exercice7 config_perso.in input_scan=[valeur]")
    configFile.process(argv[i]);

  // Parametres de simulation :
  double tfin    = configFile.get<double>("tfin");
  int nx         = configFile.get<int>("nx"); // nb intervalles
  double CFL     = configFile.get<double>("CFL");
  double nsteps  = configFile.get<double>("nsteps");
  double A       = configFile.get<double>("A"); // amplitude de l'excitation
  double f_hat   = configFile.get<double>("f_hat");
  double n_init  = configFile.get<double>("n_init");
  double hL      = configFile.get<double>("hL");
  double hR      = configFile.get<double>("hR");
  double h00     = configFile.get<double>("h00"); // profondeur, cas uniforme
  double x1      = configFile.get<double>("x1"); //x_L
  double x2      = configFile.get<double>("x2"); //x_R
  double xa      = configFile.get<double>("xa");
  double xb      = configFile.get<double>("xb");
  double L       = configFile.get<double>("L");
  double om      = configFile.get<double>("om"); // pulsation de l'excitation
  int n_stride(configFile.get<int>("n_stride"));

  int N = nx+1;                                // nb pts de maillage

// Conditions aux bords:
  string bc_l           = configFile.get<string>("cb_gauche");
  string bc_r           = configFile.get<string>("cb_droite");

// Type de forme initiale de la vague: selon donnée Eq.(4) ou mode propre
// ('mode' pour mode propre, autrement Eq.(4))
  string initialization = configFile.get<string>("initialization"); 

// Onde partant vers la gauche ou vers la droite ou statique
// (par exemple 'left', 'right', 'static')
  string initial_state = configFile.get<string>("initial_state");

// Selecteur pour le cas h0 uniforme:
  bool v_uniform        = configFile.get<bool>("v_uniform");

// Selecteur pour choisir le pas de temps:
// true --> dt=tfin/nsteps; t final est exactement tfin
// false --> dt tel que beta_CFL=1; attention, t final n'est pas exactement tfin
  bool impose_nsteps    = configFile.get<bool>("impose_nsteps");

  vector<double> h0(N) ;
  vector<double> vel2(N) ;
  vector<double> x(N) ;
  vector<double> fpast(N), fnow(N), fnext(N), beta2(N);

  dx = L / (N-1);
  bool ecrire_f = configFile.get<bool>("ecrire_f"); // Exporter f(x,t) ou non
 // Eq.(1) ou Eq.(2) [ou Eq.(6) (faculattif)]: Eq1, Eq2 ou Eq6
  string equation_type = configFile.get<string>("equation_type");
  

  for (int i(0); i < N; ++i) { 
    x[i] = i * dx;
    h0[i] = 0.0;

    //if (x[i] < x1 || x[i] > x2) {
    //    throw std::out_of_range("Erreur : x[" + std::to_string(i) + "] = " + std::to_string(x[i]) + " est hors des bornes [" + std::to_string(x1) + ", " + std::to_string(x2) + "].");
    //}
    printf("x1 :%f, x2 : %f, xa : %f, xb : %f\n", x1, x2, xa, xb);
    if (v_uniform) {
        h0[i] = h00;
    } else {
        if (0 <= x[i] && x[i] <= xa) {
            h0[i] = hL;
            printf("je suis dans la premiere condition, h0[i] = %f\n", h0[i]);
        } else if (xa < x[i] && x[i] < xb) {
            h0[i] = 0.5 * (hL + hR) + 0.5 * (hL - hR) * cos(PI * (x[i] - xa) / (xb - xa));
        } else if (xb <= x[i] && x[i] <= L) {
            h0[i] = hR;
        }
    }

    vel2[i] = g * h0[i];
}
printf("h0[0] is %f\n", h0[0]);

  // maiximal value of u^2 (to be used to set dt)
  auto max_vel2 = std::max_element(vel2.begin(), vel2.end());
  double max_vel2_it = *max_vel2;
  // TODO: set dt for given CFL
  dt = CFL * dx / sqrt(max_vel2_it); // MODIFY

  // TODO: define dt and CFL with given nsteps
  if(impose_nsteps){
    dt  = tfin /nsteps; // MODIFY
    CFL = sqrt(max_vel2_it) * dt/dx; // MODIFY
  }

  // Fichiers de sortie :
  string output = configFile.get<string>("output");

  ofstream fichier_x((output + "_x").c_str());
  fichier_x.precision(15);

  ofstream fichier_v((output + "_v").c_str());
  fichier_v.precision(15);

  ofstream fichier_f((output + "_f").c_str());
  fichier_f.precision(15);

  ofstream fichier_en((output + "_en").c_str());
  fichier_en.precision(15);

  // Initialisation des tableaux du schema numerique :

  //TODO initialize f and beta
  for(int i(0); i<N; ++i)
  {
    fpast[i] = 0.;
    fnow[i]  = 0.;
    beta2[i] = vel2[i]*dt*dt/(dx*dx); // TODO: Modifier pour calculer beta^2 aux points de maillage

    fnow[i]  = finit(x[i], n_init,  L, f_hat, x1, x2, initialization);

    if(initial_state =="static"){
      fpast[i] = fnow[i]; // TODO: system is at rest for t<=0, 
    }
    else if(initial_state =="right"){ 
      fpast[i] = finit(x[i]-sqrt(abs(g*h0[i]))*dt, n_init,  L, f_hat, x1, x2, initialization); // TODO: propagation to the right
    }
    else if(initial_state =="left"){
      fpast[i] = finit(x[i]+sqrt(abs(g*h0[i]))*dt, n_init,  L, f_hat, x1, x2, initialization); // TODO: propagation to the left
    }
  }


  cout<<"beta2[0] is "<<beta2[0]<<endl;
  cout<<"dt is "<< dt <<endl;


  // Boucle temporelle :
  for(t=0.; t<tfin-.5*dt; t+=dt)
  {
    // Ecriture :
    if(stride%n_stride == 0)
    {
      if(ecrire_f) fichier_f << t << " " << fnow << endl;
      fichier_en << t << " " << energy(fnow,dx) << endl;
     }
    ++stride;

   
    /*
    // Evolution :
    for(int i(1); i<N-1; ++i)
    { 
        if(equation_type=="A"){
            fnext[i]= (dt*dt/(dx*dx))*(g*h0[i])*(fnow[i+1]-2*fnow[i]+fnow[i-1])+2*fnow[i]-fpast[i] ;
        }
        else if(equation_type=="B"){
           fnext[i]= (dt*dt/(4*dx*dx))*(g*h0[i+1] -g*h0[i-1])*(fnow[i+1]-fnow[i-1]) + beta2[i]*(fnow[i+1]-2*fnow[i]+fnow[i-1]) +2*fnow[i]-fpast[i]; 
           //fnext[i] = 2 * fnow[i] - fpast[i] + (beta2[i+1] * (fnow[i+1] - fnow[i])  -  beta2[i]   * (fnow[i]   - fnow[i-1]));
          }
        else if(equation_type=="C"){
            fnext[i]= (dt*dt/(dx*dx))*(g*h0[i+1]*fnow[i+1]-2*g*h0[i]*fnow[i] + g*h0[i-1]*fnow[i-1]) + 2*fnow[i] - fpast[i];
        }
        else{
            cerr << "Merci de choisir une equation valide" << endl;
        }
    
    }
*/
// Evolution PAS SUREE:
for(int i(1); i<N-1; ++i)
{
  fnext[i] = 0.0; // TODO : Schémas pour les 3 cas, Equation A ou B ou C
  if(equation_type == "A") {
    fnext[i] = 2*(1-beta2[i])*fnow[i] - fpast[i] + beta2[i]*(fnow[i+1]+ fnow[i-1]);
}
else if(equation_type == "B") {
  fnext[i] = 2*(1-beta2[i])*fnow[i] - fpast[i] + beta2[i]*(fnow[i+1]+ fnow[i-1]) //mm que A
  + 0.25*(beta2[i]/vel2[i])*(fnow[i+1]-fnow[i-1])*(vel2[i+1]-vel2[i-1]);
    }
else if(equation_type == "C") {
  //fnext[i] = 2*(1-beta2[i])*fnow[i] - fpast[i] + beta2[i]*(fnow[i+1]+ fnow[i-1]) //mm que A //mmque B
  //+ 0.5*(beta2[i]/vel2[i])*(fnow[i+1]-fnow[i-1])*(vel2[i+1]-vel2[i-1]) // mm que B
  //+ (beta2[i]/vel2[i])*fnow[i]*(vel2[i+1]-2*vel2[i]-vel2[i-1]);
  fnext[i] = fnow[i]*(2.0-4.0*beta2[i]+beta2[i+1]+beta2[i-1]) - fpast[i] + 0.5*(beta2[i+1]-beta2[i-1])*(fnow[i+1]-fnow[i-1])+beta2[i]*(fnow[i+1]+fnow[i-1]);}
  else {
    cerr << "Merci de choisir une equation valide" << endl;
    exit(1);
  }
}


    // Impose boundary conditions
    boundary_condition(fnext, fnow, A, om, t, dt, beta2, bc_l, bc_r, N);

    // Mise a jour et préparer le pas suivant:
    fpast = fnow;
    fnow  = fnext;
  }

  if(ecrire_f) fichier_f << t << " " << fnow << endl; //valeur de f a chaque instant : onde 
  fichier_x << x << endl; // coordonnees du maillage
  fichier_v << vel2 << endl; //vitesse de l'onde 
  fichier_en << t << " " << energy(fnow,dx) << endl; //energie totale a chaque pas de temps 

  fichier_f.close();
  fichier_x.close();
  fichier_v.close();
  fichier_en.close();

  return 0;
}
