// copy from by Szymon Winczewski tersoff potential
// partially based on LAMMPS pair_tersoff library

#ifndef pot_en_atom_h
#define pot_en_atom_h

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>

#include "vec3d.h"

//#define FORCE_FIELD_DEBUG


const double HALF_OF_PI = 2.0 * atan(1.0);
typedef int64_t bigint;


class ForceFieldAtomistic
{
private:
    int number_of_atoms;
    int* n_num;
    int** n_list;
    vec3d** n_bonds;
    double** n_fs;

// parametry (komentarz rzedstawia nazewnictwo w LAMMPS-ie)
    std::string potential_file_name_;
    double param_m, param_n; // powerm, powern
    double param_beta; // beta
    double param_gamma; // gamma
    double param_lambda1; // lam1
    double param_lambda2; // lam2
    double param_lambda3; // lam3
    double param_c, param_d;
    double param_cos_theta_0;
    double param_A, param_B;
    double param_R, param_D;

// parametry pomocnicze
    double cutoff_low;
    double cutoff;
    double cutoff_squared;

    double param_c_squared;
    double param_d_squared;
    double one_plus_param_c_squared_over_param_d_squared;

    double pi_over_4_D;
    double beta_2S_D_0_over_S;

    double coeff_c1;
    double coeff_c2;
    double coeff_c3;
    double coeff_c4;


    double cutoff_func(double r_ij);
    double repulsion_func(double r_ij);
    double angular_func(double cos_theta_jik);
    double attraction_func(double r_ij);

    double cutoff_func_prime(double r_ij);
    double repulsion_func_prime(double r_ij);
    double attraction_func_prime(double r_ij);
    double func_chi_ij(double r_ij, double r_ik, double *r_ij_vec, double *r_ik_vec);
    double angular_func_prime(double cos_theta_jik);


public:
    vec3d** bd_imp(int n_nums, std::string bx, std::string by, std::string bz,
               double Lxl, double Lxh, double Lyl, double Lyh, double Lzl, double Lzh, vec3d **dist_mat);
    int** NN1_implementation(int n_atoms, vec3d **atoms_distance_matrix, double cut_off);
    void read_parameters(std::string file_name);
    void ini_system(std::string atom_file);
    int ini_system_bin(std::string file_name_, int i);
    double get_cutoff();
    void define_system(int number_of_atoms_, int **n_list_, vec3d **n_bonds_);
    void compute_energy();
    void compute_force();
    //constructor and destructor
    ForceFieldAtomistic();
    ~ForceFieldAtomistic();
};

#endif
