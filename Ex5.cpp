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

double energy(const vector<double>& fnow, double dx) {
    double ener = 0.0;
    for (size_t i = 0; i < fnow.size(); ++i) {
        ener += fnow[i] * fnow[i];
    }
    ener *= dx;
    return ener;
}

void boundary_condition(vector<double> &fnext, vector<double> &fnow, double const& A, double om,
                        double const& t, double const& dt,
                        vector<double> &vel2, double dx,
                        string &bc_l, string &bc_r, int &N) {
    if (bc_l == "fixe") {
        fnext[0] = 0.0;
    } else if (bc_l == "libre") {
        fnext[0] = fnext[1];
    } else if (bc_l == "sortie") {
        double c = sqrt(vel2[0]);
        fnext[0] = fnow[1] + (c*dt - dx)/(c*dt + dx) * (fnext[1] - fnow[0]);
    } else if (bc_l == "excitation") {
        fnext[0] = A * sin(om * t);
    } else {
        cerr << "Merci de choisir une condition au bord gauche valide" << endl;
    }

    if (bc_r == "fixe") {
        fnext[N-1] = 0.0;
    } else if (bc_r == "libre") {
        fnext[N-1] = fnext[N-2];
    } else if (bc_r == "sortie") {
        double c = sqrt(vel2[N-1]);
        fnext[N-1] = fnow[N-2] + (c*dt - dx)/(c*dt + dx) * (fnext[N-2] - fnow[N-1]);
    } else if (bc_r == "excitation") {
        fnext[N-1] = A * sin(om * t);
    } else {
        cerr << "Merci de choisir une condition au bord droit valide" << endl;
    }
}

double finit(double x, double n_init, double L, double f_hat, double x1, double x2, string initialization) {
    const double PI = 3.1415926535897932384626433832795028841971e0;
    if (initialization == "mode") {
        return 0.0;
    } else {
        if (x <= x1 || x >= x2) {
            return 0.0;
        } else {
            return 0.5 * f_hat * (1 - cos(2 * PI * (x - x1) / (x2 - x1)));
        }
    }
}

template <class T> ostream& operator<< (ostream& o, vector<T> const& v) {
    unsigned int len(v.size());
    for (unsigned int i = 0; i < (len - 1); ++i)
        o << v[i] << " ";
    if (len > 0)
        o << v[len-1];
    return o;
}

int main(int argc, char* argv[]) {
    const double PI = 3.1415926535897932384626433832795028841971e0;
    const double g  = 9.81;
    double dx;
    double dt;
    double t;
    double Nsteps;
    int stride(0);

    string inputPath("configuration.in");
    if(argc > 1)
        inputPath = argv[1];

    ConfigFile configFile(inputPath);
    for(int i = 2; i < argc; ++i)
        configFile.process(argv[i]);

    double tfin = configFile.get<double>("tfin");
    int nx = configFile.get<int>("nx");
    double CFL = configFile.get<double>("CFL");
    double nsteps = configFile.get<double>("nsteps");
    double A = configFile.get<double>("A");
    double f_hat = configFile.get<double>("f_hat");
    double n_init = configFile.get<double>("n_init");
    double hL = configFile.get<double>("hL");
    double hR = configFile.get<double>("hR");
    double h00 = configFile.get<double>("h00");
    double x1 = configFile.get<double>("x1");
    double x2 = configFile.get<double>("x2");
    double xa = configFile.get<double>("xa");
    double xb = configFile.get<double>("xb");
    double L = configFile.get<double>("L");
    double om = configFile.get<double>("om");
    int n_stride(configFile.get<int>("n_stride"));

    int N = nx + 1;
    string bc_l = configFile.get<string>("cb_gauche");
    string bc_r = configFile.get<string>("cb_droite");
    string initialization = configFile.get<string>("initialization");
    string initial_state = configFile.get<string>("initial_state");
    bool v_uniform = configFile.get<bool>("v_uniform");
    bool impose_nsteps = configFile.get<bool>("impose_nsteps");

    vector<double> h0(N);
    vector<double> vel2(N);
    vector<double> x(N);
    vector<double> fpast(N), fnow(N), fnext(N), beta2(N);

    dx = L / (N - 1);
    bool ecrire_f = configFile.get<bool>("ecrire_f");
    string equation_type = configFile.get<string>("equation_type");

    for(int i = 0; i < N; ++i) {
        x[i] = i * dx;
        if(v_uniform) {
            h0[i] = h00;
        } else {
            if (x[i] <= xa) {
                h0[i] = hL;
            } else if (x[i] >= xb) {
                h0[i] = hR;
            } else {
                double cosine = cos(PI * (x[i] - xa) / (xb - xa));
                h0[i] = 0.5 * (hL + hR) + 0.5 * (hL - hR) * cosine;
            }
        }
        vel2[i] = g * h0[i];
    }

    auto max_vel2 = std::max_element(vel2.begin(), vel2.end());
    if (!impose_nsteps) {
        dt = CFL * dx / sqrt(*max_vel2);
    } else {
        dt = tfin / nsteps;
        CFL = dt * sqrt(*max_vel2) / dx;
    }

    string output = configFile.get<string>("output");
    ofstream fichier_x((output + "_x").c_str()); fichier_x.precision(15);
    ofstream fichier_v((output + "_v").c_str()); fichier_v.precision(15);
    ofstream fichier_f((output + "_f").c_str()); fichier_f.precision(15);
    ofstream fichier_en((output + "_en").c_str()); fichier_en.precision(15);

    for(int i = 0; i < N; ++i) {
        fnow[i] = finit(x[i], n_init, L, f_hat, x1, x2, initialization);
        beta2[i] = dt * dt * vel2[i] / (dx * dx);

        if (initial_state == "static") {
            fpast[i] = fnow[i];
        } else if (initial_state == "right") {
            fpast[i] = fnow[i] - dt * sqrt(vel2[i]) * (fnow[min(N-1, i+1)] - fnow[max(0, i-1)]) / (2 * dx);
        } else if (initial_state == "left") {
            fpast[i] = fnow[i] + dt * sqrt(vel2[i]) * (fnow[min(N-1, i+1)] - fnow[max(0, i-1)]) / (2 * dx);
        }
    }

    for(t = 0.; t < tfin - 0.5 * dt; t += dt) {
        if(stride % n_stride == 0) {
            if(ecrire_f) fichier_f << t << " " << fnow << endl;
            fichier_en << t << " " << energy(fnow, dx) << endl;
        }
        ++stride;

        for(int i = 1; i < N-1; ++i) {
            if (equation_type == "Eq1") {
                fnext[i] = 2*fnow[i] - fpast[i] + beta2[i] * (fnow[i+1] - 2*fnow[i] + fnow[i-1]);
            } else if (equation_type == "Eq2") {
                double dudx_right = (vel2[i+1] + vel2[i]) * (fnow[i+1] - fnow[i]);
                double dudx_left  = (vel2[i]   + vel2[i-1]) * (fnow[i] - fnow[i-1]);
                fnext[i] = 2*fnow[i] - fpast[i] + dt*dt / (4*dx*dx) * (dudx_right - dudx_left);
            } else if (equation_type == "Eq6") {
                double u2f_right = vel2[i+1]*fnow[i+1];
                double u2f_center = vel2[i]*fnow[i];
                double u2f_left = vel2[i-1]*fnow[i-1];
                fnext[i] = 2*fnow[i] - fpast[i] + dt*dt/(dx*dx) * (u2f_right - 2*u2f_center + u2f_left);
            }
        }

        boundary_condition(fnext, fnow, A, om, t, dt, vel2, dx, bc_l, bc_r, N);
        fpast = fnow;
        fnow = fnext;
    }

    if(ecrire_f) fichier_f << t << " " << fnow << endl;
    fichier_x << x << endl;
    fichier_v << vel2 << endl;
    fichier_en << t << " " << energy(fnow, dx) << endl;

    fichier_f.close(); fichier_x.close(); fichier_v.close(); fichier_en.close();
    return 0;
}