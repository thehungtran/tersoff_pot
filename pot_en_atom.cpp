#include "pot_en_atom.h"

using namespace std;

ForceFieldAtomistic::ForceFieldAtomistic()
{
    cout << "object is being created" << endl;
}
ForceFieldAtomistic::~ForceFieldAtomistic()
{
    cout << "object is being deleted" << endl;
}

void ForceFieldAtomistic::read_parameters(string file_name)
{
    ifstream input_file;
    bool end_keyword_found;
    string keywords;

    input_file.open(file_name.c_str());
    if (input_file.is_open() == 0)
    {
        cout << "error: could not locate potential file" << endl;
        exit(0);
    }

    bool param_initialized[14];
    for (int param = 0; param < 14; param++)
        param_initialized[param] = 0;

    end_keyword_found = 0;
    input_file >> keywords;
    do
    {
        if ( keywords == "m" )
        {
            input_file >> param_m;
            param_initialized[0] = 1;
        }
        else if ( keywords == "gamma" )
        {
            input_file >> param_gamma;
            param_initialized[1] = 1;
        }
        else if ( keywords == "lambda3" )
        {
            input_file >> param_lambda3;
            param_initialized[2] = 1;
        }
        else if ( keywords == "c" )
        {
            input_file >> param_c;
            param_initialized[3] = 1;
        }
        else if ( keywords == "d" )
        {
            input_file >> param_d;
            param_initialized[4] = 1;
        }
        else if ( keywords == "cos_theta_0" )
        {
            input_file >> param_cos_theta_0;
            param_initialized[5] = 1;
        }
        else if ( keywords == "n" )
        {
            input_file >> param_n;
            param_initialized[6] = 1;
        }
        else if ( keywords == "beta" )
        {
            input_file >> param_beta;
            param_initialized[7] = 1;
        }
        else if ( keywords == "lambda2" )
        {
            input_file >> param_lambda2;
            param_initialized[8] = 1;
        }
        else if ( keywords == "B" )
        {
            input_file >> param_B;
            param_initialized[9] = 1;
        }
        else if ( keywords == "R" )
        {
            input_file >> param_R;
            param_initialized[10] = 1;
        }
        else if ( keywords == "D" )
        {
            input_file >> param_D;
            param_initialized[11] = 1;
        }
        else if ( keywords == "lambda1" )
        {
            input_file >> param_lambda1;
            param_initialized[12] = 1;
        }
        else if ( keywords == "A" )
        {
            input_file >> param_A;
            param_initialized[13] = 1;
        }
        else if ( keywords == "end" )
            end_keyword_found = 1;
        else
            cout << "unrecognized keyword \"" << keywords << "\" during parsing the config file!" << endl;
        input_file >> keywords;
    }
    while ( ( end_keyword_found == 0 ) && ( input_file.eof() == 0 ) );

    input_file.close();

    if ( ( ( param_m != 1.0 ) && ( param_m != 3.0 ) ) ||
            ( param_n < 0.0 ) ||
            ( param_beta < 0.0 ) ||
            ( param_gamma < 0.0 ) ||
            ( param_lambda1 < 0.0 ) ||
            ( param_lambda2 < 0.0 ) ||
            ( param_c < 0.0 ) ||
            ( param_d < 0.0 ) ||
            ( param_A < 0.0 ) ||
            ( param_B < 0.0 ) ||
            ( param_R < 0.0 ) ||
            ( param_D < 0.0 ) ||
            ( param_D > param_R ) )
    {
        cout << "error: incorrect potential parameters!" << endl;
        exit(0);
    }

    for (int param = 0; param < 14; param++)
        if ( param_initialized[param] == 0 )
        {
            cout << "error: some Tersoff parameters were not initialized" << endl;
            exit(0);
        }

    #ifdef FORCE_FIELD_DEBUG
    cout << "parameters: " << endl;
    cout << "  m           " << param_m << endl;
    cout << "  gamma       " << param_gamma << endl;
    cout << "  lambda3     " << param_lambda3 << endl;
    cout << "  c           " << param_c << endl;
    cout << "  d           " << param_d << endl;
    cout << "  cos_theta_0 " << param_cos_theta_0 << endl;
    cout << "  n           " << param_n << endl;
    cout << "  beta        " << param_beta << endl;
    cout << "  lambda2     " << param_lambda2 << endl;
    cout << "  B           " << param_B << endl;
    cout << "  R           " << param_R << endl;
    cout << "  D           " << param_D << endl;
    cout << "  lambda1     " << param_lambda1 << endl;
    cout << "  A           " << param_A << endl;
    #endif

    cutoff = param_R + param_D;
    if ( cutoff <= 0.0 )
    {
        cout << "error: incorrect cutoff radius" << endl;
        exit(0);
    }

    cutoff_squared = cutoff * cutoff;
    cutoff_low = param_R - param_D;

    param_c_squared = param_c * param_c;
    param_d_squared = param_d * param_d;
    one_plus_param_c_squared_over_param_d_squared = 1.0 + param_c_squared / param_d_squared;

    coeff_c1 = pow( 2.0 * param_n * 1.0e-16, -1.0 / param_n );
    coeff_c2 = pow( 2.0 * param_n * 1.0e-8,  -1.0 / param_n );
    coeff_c3 = 1.0 / coeff_c2;
    coeff_c4 = 1.0 / coeff_c1;

    // constants term for central force
    pi_over_4_D = HALF_OF_PI * 0.5 / param_D;
    //beta_2S_D_0_over_S = param_beta * param_D_0 * sqrt(2 * param_S) /(param_S - 1);

}

double ForceFieldAtomistic::get_cutoff()
{
    return cutoff;
}


void ForceFieldAtomistic::define_system(int number_of_atoms_, int **n_list_, vec3d **n_bonds_)
{
    // n_list is the list of nearest neighbor of each atom
    // with the first column is the number of atoms nearest to atom i

    number_of_atoms = number_of_atoms_;
    n_list = n_list_;
    n_bonds = n_bonds_;
}


void ForceFieldAtomistic::compute_energy()
{
    #ifdef FORCE_FIELD_DEBUG
    cout << endl;
    cout << "*** ForceFieldAtomistic::compute_energy() started! ***" << endl;
    cout << endl;
    #endif

    double total_energy;

    int atom_i;
    int j, k;
    int n_num_i;
//    int *n_list_i;
    vec3d *n_bonds_i;
    double r_ij;
    double r_ij_vec[3];
    double r_ik;
    double r_ik_vec[3];

    double attractive_term;
    double b_ij;
    //double bij[number_of_atoms][number_of_atoms];


    total_energy = 0.0;

    for (atom_i = 0; atom_i < number_of_atoms; atom_i++)
    {
        n_num_i = n_list[atom_i][0];
        n_bonds_i = n_bonds[atom_i];

        for (j = 0; j < n_num_i; j++)
        {
            attractive_term = 0;
            r_ij = n_bonds_i[j].r;
            if ( r_ij > cutoff )
                continue;

            r_ij_vec[0] = n_bonds_i[j].x;
            r_ij_vec[1] = n_bonds_i[j].y;
            r_ij_vec[2] = n_bonds_i[j].z;

// repulsive term
            //total_energy += repulsion_func(r_ij);
            attractive_term = repulsion_func(r_ij);

// attractive term

            b_ij = 0;
            for (k = 0; k < n_num_i; k++)
            {
                if ( k == j )
                    continue;

                r_ik = n_bonds_i[k].r;
                if ( r_ik > cutoff )
                    continue;

                r_ik_vec[0] = n_bonds_i[k].x;
                r_ik_vec[1] = n_bonds_i[k].y;
                r_ik_vec[2] = n_bonds_i[k].z;

                b_ij += func_chi_ij(r_ij, r_ik, r_ij_vec, r_ik_vec);
                //cout << b_ij << endl;
            }

            //bij[atom_i][j] = sqrt((1 + b_ij));
            //attractive_term += 1 / sqrt((1 + b_ij)) * attraction_func(r_ij);
            attractive_term -= attraction_func(r_ij) / sqrt((1 + b_ij));
            attractive_term = attractive_term * 0.5 * cutoff_func(r_ij);
            //total_energy -= 1 / sqrt((1 + b_ij)) * attraction_func(r_ij);
            //total_energy = 0.5 * cutoff_func(r_ij) * total_energy;
            total_energy += attractive_term;
        }
    }
    //cout << bij[100][0] << " " << bij[100][1] << " " << bij[100][2] << " " << bij[100][3] << endl;


    cout << "total_energy = " << setprecision(20) << total_energy << endl;

    #ifdef FORCE_FIELD_DEBUG
    cout << endl;
    cout << "*** ForceFieldAtomistic::compute_energy() done! ***" << endl;
    cout << endl;
    #endif
}

double ForceFieldAtomistic::cutoff_func(double r_ij)
{
    if (r_ij < cutoff_low)
        return 1.0;
    else if (r_ij > cutoff)
        return 0;
    else
        return 0.5 * (1 - sin (HALF_OF_PI * (r_ij - param_R) / param_D));
}

double ForceFieldAtomistic::repulsion_func(double r_ij)
{
    return param_A * exp(- param_lambda1 * r_ij);
}

double ForceFieldAtomistic::attraction_func(double r_ij)
{
    return param_B * exp(- param_lambda2 * r_ij);
}

double ForceFieldAtomistic::angular_func(double cos_theta_jik)
{
    return param_gamma * (one_plus_param_c_squared_over_param_d_squared - param_c_squared /
                       (param_d_squared + (- param_cos_theta_0 + cos_theta_jik) * (- param_cos_theta_0 + cos_theta_jik)));
}

double ForceFieldAtomistic::func_chi_ij (double r_ij, double r_ik, double *r_ij_vec, double *r_ik_vec)
{
    double Chi_ij, cos_theta_jik;
    cos_theta_jik = ( r_ij_vec[0] * r_ik_vec[0] + r_ij_vec[1] * r_ik_vec[1] + r_ij_vec[2] * r_ik_vec[2] ) / ( r_ij * r_ik );

    Chi_ij = cutoff_func(r_ik) * exp(param_lambda3 * (r_ij - r_ik)) * angular_func(cos_theta_jik);
    return Chi_ij;
}

int **NN_implementation(int n_atoms, double **atoms_distance_matrix, double cut_off)
{
    // n_atoms : number of atoms
    // cutoff radius

    int** nn_atoms = new int*[n_atoms];
    int tmp_count;

    for (int i = 0; i < n_atoms; i++)
    {
        tmp_count = 0;
        for (int j = 0; j < n_atoms; j++)
        {
            if ((atoms_distance_matrix[i][j] <= cut_off) && (i != j))
                tmp_count++;
        }

        nn_atoms[i] = new int[tmp_count + 1];
        nn_atoms[i][0] = tmp_count;

        tmp_count = 1;
        for (int j = 0; j < n_atoms; j++)
        {
            if ((atoms_distance_matrix[i][j] <= cut_off) && (i != j))
            {
                nn_atoms[i][tmp_count] = j;
                tmp_count++;
            }
        }
    }
    return nn_atoms;
}


double ForceFieldAtomistic::cutoff_func_prime(double r_ij)
{
    if ((r_ij <= cutoff) && (r_ij >= cutoff_low))
        return - pi_over_4_D * cos (HALF_OF_PI * (r_ij - param_R) / param_D);
    else
        return 0;
}

double ForceFieldAtomistic::repulsion_func_prime(double r_ij)
{
    return  - param_lambda1 * param_A * exp(- param_lambda1 * r_ij);
}

double ForceFieldAtomistic::attraction_func_prime(double r_ij)
{
    return - param_lambda2 * param_B * exp(- param_lambda2 * r_ij);
}

double ForceFieldAtomistic::angular_func_prime(double cos_theta_jik)
{
    return - param_gamma * param_c_squared * 2 * (- param_cos_theta_0 + cos_theta_jik)
                    / ((param_d_squared +(- param_cos_theta_0 + cos_theta_jik) * (- param_cos_theta_0 + cos_theta_jik)) *
                       (param_d_squared +(- param_cos_theta_0 + cos_theta_jik) * (- param_cos_theta_0 + cos_theta_jik)));
}

void ForceFieldAtomistic::compute_force()
{
     double scalar_term;

    int atom_i;
    int k;
    int n_num_i, n_num_j, n_num_k;
    int *n_list_i, *n_list_j, *n_list_k;
    int count_ting;
//    int *n_list_i;

    //vec3d Force_[number_of_atoms][number_of_atoms];
    vec3d *n_bonds_i, *n_bonds_j, *n_bonds_k;
    double r_ij;
    double r_ij_vec[3];
    double r_ji_vec[3];
    double r_ik;
    double r_ik_vec[3];
    double r_jk;
    double r_jk_vec[3];
    double r_il;
    double r_il_vec[3];
    double r_jl;
    double r_jl_vec[3];
    double r_kk1;
    double r_kk1_vec[3];

    double C_3, C_4, C_5, C_6, C_7_8, C_9, C_10;
    double cos_theta_ijk, cos_theta_jik, cos_theta_ikj;

    //double attractive_term;
    double b_ij, b_ji, b_ik, b_jk, b_ki, b_kj;

//    Force_ = new vec3d [number_of_atoms];
count_ting = 0;
ofstream filetype;
filetype.open("result.txt");

    for (atom_i = 0; atom_i < number_of_atoms; atom_i++)
    {
        n_num_i = n_list[atom_i][0]; // number of nn of atoms i
        n_bonds_i = n_bonds[atom_i]; // distance vector of atom i to all other atoms
        n_list_i = n_list[atom_i];
        double resultant_force_x = 0.00; double resultant_force_y = 0.00; double resultant_force_z = 0.00;

        for (int j = 0; j < n_num_i; j++)
        {
            n_num_j = n_list[n_list_i[j+1]][0];
            n_bonds_j = n_bonds[n_list_i[j+1]];
            n_list_j = n_list[n_list_i[j+1]];

            scalar_term = 0.0;
            C_3 = 0; C_4 = 0; C_5 = 0; C_6 = 0; C_7_8 = 0; C_9 = 0; C_10 = 0;

            r_ij = n_bonds_i[j].r ;

            r_ij_vec[0] = n_bonds_i[j].x;
            r_ij_vec[1] = n_bonds_i[j].y;
            r_ij_vec[2] = n_bonds_i[j].z;

            r_ji_vec[0] = - r_ij_vec[0];
            r_ji_vec[1] = - r_ij_vec[1];
            r_ji_vec[2] = - r_ij_vec[2];


            b_ij = 0;
            for (k = 0; k < n_num_i; k++)
            {
                if ( k == j )
                    continue;

                n_bonds_k = n_bonds[k];
                n_num_k   = n_list[k][0];

                r_ik = n_bonds_i[k].r;
                if ( r_ik > cutoff )
                    continue;

                r_ik_vec[0] = n_bonds_i[k].x;
                r_ik_vec[1] = n_bonds_i[k].y;
                r_ik_vec[2] = n_bonds_i[k].z;

                b_ij += func_chi_ij(r_ij, r_ik, r_ij_vec, r_ik_vec);
            }
            b_ij = 1 / sqrt(1 + b_ij);


            b_ji = 0;
            for (k = 0; k < n_num_j; k++)
            {
                if ( n_list_j[k+1] == atom_i )
                    continue;

                r_jk = n_bonds_j[k].r;
                if ( r_jk > cutoff )
                    continue;

                r_jk_vec[0] = n_bonds_j[k].x;
                r_jk_vec[1] = n_bonds_j[k].y;
                r_jk_vec[2] = n_bonds_j[k].z;

                b_ji += func_chi_ij(r_ij, r_jk, r_ji_vec, r_jk_vec);
            }
            b_ji = 1 / sqrt(1 + b_ji);

            // the first two terms of the summation
            scalar_term += cutoff_func_prime(r_ij) *
            (repulsion_func(r_ij) - 0.5 * (b_ij + b_ji) * attraction_func(r_ij)) +
            cutoff_func(r_ij) * (repulsion_func_prime(r_ij) - 0.5 * (b_ij + b_ji) * attraction_func_prime(r_ij));
            // the rest of the summation will be computed here

            // C_3, C_5, C_7
            for (k = 0; k < n_num_i; k++)
            {
                n_num_k = n_list[n_list_i[k+1]][0];
                n_bonds_k = n_bonds[n_list_i[k+1]];
                 if ( k == j )
                    continue;

                r_ik = n_bonds_i[k].r;
                if ( r_ik > cutoff )
                    continue;

                r_ik_vec[0] = n_bonds_i[k].x;
                r_ik_vec[1] = n_bonds_i[k].y;
                r_ik_vec[2] = n_bonds_i[k].z;

                cos_theta_jik = (r_ij_vec[0] * r_ik_vec[0] + r_ij_vec[1] * r_ik_vec[1] + r_ij_vec[2] * r_ik_vec[2])/
                (r_ij * r_ik);

                b_ik = 0;
                for (int l = 0; l < n_num_i; l++)
                {
                    if  (l == k)
                        continue;

                    r_il = n_bonds_i[l].r;
                    if ( r_il > cutoff )
                    continue;

                    r_il_vec[0] = n_bonds_i[l].x;
                    r_il_vec[1] = n_bonds_i[l].y;
                    r_il_vec[2] = n_bonds_i[l].z;


                    b_ik += func_chi_ij(r_ik, r_il, r_ik_vec, r_il_vec);
                }
                b_ik = 1 / (sqrt(1 + b_ik) * (1 + b_ik));


                C_3 += 0.25 * cutoff_func(r_ij) * cutoff_func(r_ik) * attraction_func(r_ij) *
                exp(param_lambda3 * (r_ij - r_ik)) * b_ij * b_ij * b_ij *
                (angular_func_prime(cos_theta_jik) * cos_theta_jik / r_ij -
                 angular_func_prime(cos_theta_jik) / r_ik -
                 param_lambda3 * angular_func((cos_theta_jik)));

                C_5 += 0.25 * cutoff_func(r_ij) * cutoff_func(r_ik) * attraction_func(r_ik) *
                exp(param_lambda3 * (r_ik - r_ij)) * b_ik *
                (param_lambda3 * angular_func(cos_theta_jik) +
                 angular_func_prime(cos_theta_jik) * cos_theta_jik / r_ij -
                 angular_func_prime(cos_theta_jik) / r_ik);

                C_9 += - 0.25 * cutoff_func_prime(r_ij) * cutoff_func(r_ik) * attraction_func(r_ik) *
                exp(param_lambda3 * (r_ik - r_ij)) * angular_func(cos_theta_jik) * b_ik;

            }

            // C_4, C_6, C_8
            for (k = 0; k < n_num_j; k++)
            {
                n_num_k = n_list[n_list_j[k+1]][0];
                n_bonds_k = n_bonds[n_list_j[k+1]];
                n_list_k = n_list[n_list_j[k+1]];
                if ( n_list_j[k+1] == atom_i )
                    continue;

                r_jk = n_bonds_j[k].r;
                if ( r_jk > cutoff )
                    continue;

                r_jk_vec[0] = n_bonds_j[k].x;
                r_jk_vec[1] = n_bonds_j[k].y;
                r_jk_vec[2] = n_bonds_j[k].z;

                cos_theta_ijk = (r_ji_vec[0] * r_jk_vec[0] + r_ji_vec[1] * r_jk_vec[1] + r_ji_vec[2] * r_jk_vec[2]) /
                (r_jk * r_ij);

                b_jk = 0;
                for (int l = 0; l < n_num_j; l++)
                {
                    if (l == k)
                        continue;

                    r_jl = n_bonds_j[l].r;
                    if ( r_jl > cutoff )
                        continue;

                    r_jl_vec[0] = n_bonds_j[l].x;
                    r_jl_vec[1] = n_bonds_j[l].y;
                    r_jl_vec[2] = n_bonds_j[l].z;


                    b_jk += func_chi_ij(r_jk, r_jl, r_jk_vec, r_jl_vec);
                }
                b_jk = 1 / (sqrt(1 + b_jk) * (1 + b_jk));


                C_4 += 0.25 * cutoff_func(r_ij) * cutoff_func(r_jk) * attraction_func(r_ij) *
                exp(param_lambda3 * (r_ij - r_jk)) * b_ji * b_ji * b_ji *
                (angular_func_prime(cos_theta_ijk) * cos_theta_ijk / r_ij -
                 angular_func_prime(cos_theta_ijk) / r_jk -
                 param_lambda3 * angular_func(cos_theta_ijk));

                C_6 += 0.25 * cutoff_func(r_ij) * cutoff_func(r_jk) * attraction_func(r_jk) *
                exp(param_lambda3 * (r_jk - r_ij)) * b_jk *
                (param_lambda3 * angular_func(cos_theta_ijk) +
                 angular_func_prime(cos_theta_ijk) * cos_theta_ijk / r_ij -
                 angular_func_prime(cos_theta_ijk) / r_jk);

                C_10 += - 0.25 * cutoff_func_prime(r_ij) * cutoff_func(r_jk) * attraction_func(r_jk) *
                exp(param_lambda3 * (r_jk - r_ij)) * angular_func(cos_theta_ijk) * b_jk;
            }

            for (int k = 0; k < n_num_i; k++)
            {
                if (k == j)
                continue;
                for (int l = 0; l < n_num_j; l++)
                {
                    if (n_list_i[k+1] != n_list_j[l+1])
                        continue;
                    count_ting++;
                    r_ik = n_bonds_i[k].r;
                    r_ik_vec[0] = - n_bonds_i[k].x;
                    r_ik_vec[1] = - n_bonds_i[k].y;
                    r_ik_vec[2] = - n_bonds_i[k].z;

                    r_jk = n_bonds_j[l].r;
                    r_jk_vec[0] = - n_bonds_j[l].x;
                    r_jk_vec[1] = - n_bonds_j[l].y;
                    r_jk_vec[2] = - n_bonds_j[l].z;

                    cos_theta_ikj = (r_ik_vec[0] * r_jk_vec[0] + r_ik_vec[1] * r_jk_vec[1] + r_ik_vec[2] * r_jk_vec[2]) /
                    (r_ik * r_jk);

                    b_ki = 0;
                    b_kj = 0;
                    n_list_k = n_list[n_list_i[k+1]];
                    n_bonds_k = n_bonds[n_list_i[k+1]];
                    for (int k1 = 0; k1 < n_list_k[0]; k1++)
                    {
                        if (n_list_k[k1+1] == atom_i)
                            continue;

                        r_kk1 = n_bonds_k[k1].r;
                        if (r_kk1 > cutoff)
                            continue;
                        r_kk1_vec[0] = n_bonds_k[k1].x;
                        r_kk1_vec[1] = n_bonds_k[k1].y;
                        r_kk1_vec[2] = n_bonds_k[k1].z;

                        b_ki += func_chi_ij(r_ik, r_kk1, r_ik_vec, r_kk1_vec);
                    }
                    for (int k1 = 0; k1 < n_list_k[0]; k1++)
                    {
                        if (n_list_k[k1+1] == n_list_i[j+1])
                            continue;

                        r_kk1 = n_bonds_k[k1].r;
                        if (r_kk1 > cutoff)
                            continue;
                        r_kk1_vec[0] = n_bonds_k[k1].x;
                        r_kk1_vec[1] = n_bonds_k[k1].y;
                        r_kk1_vec[2] = n_bonds_k[k1].z;

                        b_kj += func_chi_ij(r_jk, r_kk1, r_jk_vec, r_kk1_vec);
                    }
                    b_ki = 1 / (sqrt(1 + b_ki) * (1 + b_ki));
                    b_kj = 1 / (sqrt(1 + b_kj) * (1 + b_kj));

                    C_7_8 += 0.25 * cutoff_func(r_ik) * cutoff_func(r_jk) * angular_func_prime(cos_theta_ikj) *
                    r_ij / (r_ik * r_jk) *
                    (attraction_func(r_ik) * exp(param_lambda3 * (r_ik - r_jk)) * b_ki +
                     attraction_func(r_jk) * exp(param_lambda3 * (r_jk - r_ik)) * b_kj);
                }
            }

            scalar_term += C_3 + C_4 + C_5 + C_6 + C_7_8 + C_9 + C_10;
            // finalize the central force
            //if (atom_i < n_list_i[j+1])
            resultant_force_x += scalar_term * r_ij_vec[0] / r_ij;
            resultant_force_y += scalar_term * r_ij_vec[1] / r_ij;
            resultant_force_z += scalar_term * r_ij_vec[2] / r_ij;
        }
        filetype << atom_i + 1 <<" " << resultant_force_x - n_fs[atom_i][0]<<" " << resultant_force_y - n_fs[atom_i][1]<< " "
        << resultant_force_z - n_fs[atom_i][2] <<";"<< endl;
       // filetype << atom_i + 1 <<" " << resultant_force_x <<" " << resultant_force_y << " "
       // << resultant_force_z <<";"<< endl;

    }

    filetype.close();
}
void ForceFieldAtomistic::ini_system(string atom_file)
{
  ifstream input_sample;
  input_sample.open(atom_file.c_str());

  int n_nums;
  int Ct, idx;
  double** n_atoms;
  double Lxl, Lxh, Lyl, Lyh, Lzl, Lzh; // box bounds
  string bx, by, bz, x1, x2, y1, y2, z1, z2;


  // read from sample (works)
  string keys, tmps;

  Ct = 0;
  do
  {
      input_sample >> keys;
      cout << keys << endl;

      if (keys == "NUMBER")
      {
          input_sample >> tmps;
          input_sample >> tmps;
          input_sample >> n_nums;
          n_atoms = new double*[n_nums];
          for (int i = 0; i < n_nums; i++)
          {
              n_atoms[i] = new double [3];
          }
      }
      else if (keys =="z")
      {
          for (int i = 0; i < n_nums; i++)
          {
            input_sample >> idx;
            input_sample >> tmps;
            input_sample >> n_atoms[idx - 1][0];
            input_sample >> n_atoms[idx - 1][1];
            input_sample >> n_atoms[idx - 1][2];
            Ct++;
          }
      }
      else if (keys == "BOX")
      {
              input_sample >> tmps;
              input_sample >> bx;
              input_sample >> by;
              input_sample >> bz;

             /* input_sample >> Lxl;
              input_sample >> Lxh;
              input_sample >> Lyl;
              input_sample >> Lyh;
              input_sample >> Lzl;
              input_sample >> Lzh;
              input_sample >> x1;
              input_sample >> x2;
              input_sample >> y1;
              input_sample >> y2;
              input_sample >> z1;
              input_sample >> z2;*/
      }
  }
  while (!input_sample.eof() && Ct < n_nums);

  Lxl = - 0.29210298;
  Lxh = 17.33210298;
  Lyl = -0.32336347;
  Lyh = 17.52236347;
  Lzl = -50.00;
  Lzh = 50.00;
  /*Lxl = 0.00;
  Lxh = 17.0399999999;
  Lyl = 0.00;
  Lyh = 17.199000002;
  Lzl = -50.00;
  Lzh = 50.00;*/

  cout << Lxl << " " << n_nums << endl;

  // check if any particle is outside the box
  int abn = 0;
  for (int i = 0; i < n_nums; i++)
  {
      if (n_atoms[i][0] > Lxh)
      {
          n_atoms[i][0] = n_atoms[i][0] - (Lxh - Lxl);
          abn++;
      }
      else if (n_atoms[i][0] < Lxl)
      {
          n_atoms[i][0] = n_atoms[i][0] + (Lxh - Lxl);
          abn++;
      }

      if (n_atoms[i][1] > Lyh)
      {
          n_atoms[i][1] = n_atoms[i][1] - (Lyh - Lyl);
          abn++;
      }
      else if (n_atoms[i][1] < Lyl)
      {
          n_atoms[i][1] = n_atoms[i][1] + (Lyh - Lyl);
          abn++;
      }

      if (n_atoms[i][2] > Lzh)
      {
          n_atoms[i][2] = n_atoms[i][2] - (Lzh - Lzl);
          abn++;
      }
      else if (n_atoms[i][2] < Lzl)
      {
          n_atoms[i][2] = n_atoms[i][2] + (Lzh - Lzl);
          abn++;
      }
  }

  int** near_neigh; // n_list

  vec3d** dist_mat = new vec3d*[n_nums];
  vec3d** nn_dist = new vec3d*[n_nums];

  for (int i = 0; i < n_nums; i++)
  {
      dist_mat[i] = new vec3d[n_nums];
      for (int j = 0; j < n_nums; j++)
      {
          dist_mat[i][j].x = n_atoms[j][0] - n_atoms[i][0];
          dist_mat[i][j].y = n_atoms[j][1] - n_atoms[i][1];
          dist_mat[i][j].z = n_atoms[j][2] - n_atoms[i][2];
          dist_mat[i][j].r = sqrt((n_atoms[j][0] - n_atoms[i][0]) * (n_atoms[j][0] - n_atoms[i][0]) +
                                  (n_atoms[j][1] - n_atoms[i][1]) * (n_atoms[j][1] - n_atoms[i][1]) +
                                  (n_atoms[j][2] - n_atoms[i][2]) * (n_atoms[j][2] - n_atoms[i][2]));
      }
  }
  dist_mat = bd_imp(n_nums, bx, by, bz, Lxl, Lxh, Lyl, Lyh, Lzl, Lzh, dist_mat);

  near_neigh = NN1_implementation(n_nums, dist_mat, cutoff);

    for (int i = 0; i < n_nums; i++)
    {
        nn_dist[i] = new vec3d [near_neigh[i][0]];
        for (int j = 0; j < near_neigh[i][0]; j++)
        {
            nn_dist[i][j].x = dist_mat[i][near_neigh[i][j + 1]].x;
            nn_dist[i][j].y = dist_mat[i][near_neigh[i][j + 1]].y;
            nn_dist[i][j].z = dist_mat[i][near_neigh[i][j + 1]].z;
            nn_dist[i][j].r = dist_mat[i][near_neigh[i][j + 1]].r;
        }
    }
    //define_system(n_nums, near_neigh , nn_dist);
    number_of_atoms = n_nums;
    n_list = near_neigh;
    n_bonds = nn_dist;
    compute_energy();
//    compute_force(filetype);

    for (int i = 0; i < n_nums; i++)
    {
        delete [] dist_mat[i];
        delete [] n_atoms[i];
        delete [] nn_dist[i];
        delete [] near_neigh[i];
    }
    delete [] n_atoms;
    delete [] dist_mat;
    delete [] nn_dist;
    delete [] near_neigh;
}
vec3d** ForceFieldAtomistic::bd_imp(int n_nums, string bx, string by, string bz,
               double Lxl, double Lxh, double Lyl, double Lyh, double Lzl, double Lzh, vec3d **dist_mat)
{
     // check boundary conditions

  if ((bx == "pp") && (by == "pp") && (bz == "ff"))
  {
        double Lx, Ly;
        Lx = Lxh - Lxl;
        Ly = Lyh - Lyl;

        for (int i = 0; i < n_nums; i++)
        {
            for (int j = 0; j < n_nums; j++)
            {
                if (dist_mat[i][j].x < - Lx * 0.5)
                    dist_mat[i][j].x = dist_mat[i][j].x + Lx;
                else if (dist_mat[i][j].x >= Lx * 0.5)
                    dist_mat[i][j].x = dist_mat[i][j].x - Lx;

                if (dist_mat[i][j].y < - Ly * 0.5)
                    dist_mat[i][j].y = dist_mat[i][j].y + Ly;
                else if (dist_mat[i][j].y >= Ly * 0.5)
                    dist_mat[i][j].y = dist_mat[i][j].y - Ly;
            }
        }
  }
  else if ((bx == "pp") && (by == "pp") && (bz == "pp"))
  {
        double Lx, Ly, Lz;
        Lx = Lxh - Lxl;
        Ly = Lyh - Lyl;
        Lz = Lzh - Lzl;

        for (int i = 0; i < n_nums; i++)
        {
            for (int j = 0; j < n_nums; j++)
            {
                if (dist_mat[i][j].x < - Lx * 0.5)
                    dist_mat[i][j].x = dist_mat[i][j].x + Lx;
                else if (dist_mat[i][j].x >= Lx * 0.5)
                    dist_mat[i][j].x = dist_mat[i][j].x - Lx;

                if (dist_mat[i][j].y < - Ly * 0.5)
                    dist_mat[i][j].y = dist_mat[i][j].y + Ly;
                else if (dist_mat[i][j].y >= Ly * 0.5)
                    dist_mat[i][j].y = dist_mat[i][j].y - Ly;

                if (dist_mat[i][j].z < - Lz * 0.5)
                    dist_mat[i][j].z = dist_mat[i][j].z + Lz;
                else if (dist_mat[i][j].z >= Lz * 0.5)
                    dist_mat[i][j].z = dist_mat[i][j].z - Lz;
            }
        }
  }

  for (int i = 0; i < n_nums; i++)
  {
      for (int j = 0; j < n_nums; j++)
      {
          dist_mat[i][j].r = sqrt(dist_mat[i][j].x * dist_mat[i][j].x +
                                  dist_mat[i][j].y * dist_mat[i][j].y +
                                  dist_mat[i][j].z * dist_mat[i][j].z);
      }
  }
  return dist_mat;
}

int** ForceFieldAtomistic::NN1_implementation(int n_atoms, vec3d **atoms_distance_matrix, double cut_off)
{
    // n_atoms : number of atoms
    // cutoff radius

    int tmp_count;
    int** nn_atoms = new int*[n_atoms];
    double t1, t2;
    for (int i = 0; i < n_atoms; i++)
    {
        tmp_count = 0;
        for (int j = 0; j < n_atoms; j++)
        {
            if (i == j)
                continue;
            else if (atoms_distance_matrix[i][j].r <= cut_off)
                tmp_count++;
            for (int k = 0; k < n_atoms; k++)
            {
                if ((k == j) || (k ==i))
                    continue;
                t1 = atoms_distance_matrix[i][k].r; t2 = atoms_distance_matrix[j][k].r;
                if ((t1 <= cut_off) && (t2 <= cut_off))
                    tmp_count++;
            }
        }

        nn_atoms[i] = new int [tmp_count + 1];
        nn_atoms[i][0] = tmp_count;
        tmp_count = 1;
        for (int j = 0; j < n_atoms; j++)
        {
            if (i == j)
                continue;
            else if (atoms_distance_matrix[i][j].r <= cut_off)
            {
                nn_atoms[i][tmp_count] = j;
                tmp_count++;
            }
            for (int k = 0; k < n_atoms; k++)
            {
                if ((k == j) || (k ==i))
                    continue;
                t1 = atoms_distance_matrix[i][k].r; t2 = atoms_distance_matrix[j][k].r;
                if ((t1 <= cut_off) && (t2 <= cut_off))
                {
                    nn_atoms[i][tmp_count] = j;
                    tmp_count++;
                }
            }
        }
    }
    return nn_atoms;
}
int ForceFieldAtomistic::ini_system_bin(string file_name_, int i)
{
    ifstream file_;

    file_.open(file_name_.c_str(), ios::binary | ios::in);

    bigint timestep_bigint;
    int timestep;
    bigint number_of_atoms_bigint;
    int number_of_atoms1;
    int tmp_bounding_box_orthogonal_triclinic;
    int tmp_bounding_box_boundary_type[3][2];
    double tmp_bounding_box_xlo, tmp_bounding_box_xhi;
    double tmp_bounding_box_ylo, tmp_bounding_box_yhi;
    double tmp_bounding_box_zlo, tmp_bounding_box_zhi;
    double tmp_bounding_box_xy, tmp_bounding_box_xz, tmp_bounding_box_yz;
    int size_one; // liczba wartosci w linii
    int nchunk; // liczba segmentow poszczegolnych procesorow
    int n;
    int maxbuf = 0;
    double *buf = NULL;
    int atom_number;
    int atom_type;
    stringstream atom_name_ss;
    string atom_name;

//   int countertmp = 2;
//   for(int tstps = 0; tstps <= countertmp; tstps++)
//   {
    file_.ignore(7272*i);
// pierwszy rekord - TIMESTEP
    file_.read((char*)&timestep_bigint, sizeof(bigint));
    if ( file_.eof() == 1 )
        return 1;
    timestep = int(timestep_bigint);

// drugi rekord - NUMBER OF ATOMS

    file_.read((char*)&number_of_atoms_bigint, sizeof(bigint));
    if ( file_.eof() == 1 )
        return 2;
    number_of_atoms1 = int(number_of_atoms_bigint);

    if ( number_of_atoms1 < 1 )
    {
        #ifdef FILE_READER_EXIT_ON_ERRORS
        raiseError(ERR_FILE_READER, 13, "readLammpstrjbinFrame", "incorrect line", line_number_);
        #endif
        return 2;
    }

// trzeci rekord - TRICLINIC
    file_.read((char*)&tmp_bounding_box_orthogonal_triclinic, sizeof(int));
    if ( file_.eof() == 1 )
        return 2;

// czwarty rekord - BOUNDARY
    file_.read((char*)&tmp_bounding_box_boundary_type[0][0], sizeof(int));
    file_.read((char*)&tmp_bounding_box_boundary_type[0][1], sizeof(int));
    file_.read((char*)&tmp_bounding_box_boundary_type[1][0], sizeof(int));
    file_.read((char*)&tmp_bounding_box_boundary_type[1][1], sizeof(int));
    file_.read((char*)&tmp_bounding_box_boundary_type[2][0], sizeof(int));
    file_.read((char*)&tmp_bounding_box_boundary_type[2][1], sizeof(int));
    if ( file_.eof() == 1 )
        return 2;

// piaty rekord - BOX BOUNDS
    file_.read((char*)&tmp_bounding_box_xlo, sizeof(double));
    file_.read((char*)&tmp_bounding_box_xhi, sizeof(double));
    file_.read((char*)&tmp_bounding_box_ylo, sizeof(double));
    file_.read((char*)&tmp_bounding_box_yhi, sizeof(double));
    file_.read((char*)&tmp_bounding_box_zlo, sizeof(double));
    file_.read((char*)&tmp_bounding_box_zhi, sizeof(double));
    if ( file_.eof() == 1 )
        return 2;

    if ( tmp_bounding_box_orthogonal_triclinic == 1 )
    {
// rekord opcjonalny - nachylenie pudla symulacyjnego
        file_.read((char*)&tmp_bounding_box_xy, sizeof(double));
        file_.read((char*)&tmp_bounding_box_xz, sizeof(double));
        file_.read((char*)&tmp_bounding_box_yz, sizeof(double));
        if ( file_.eof() == 1 )
            return 2;
    }
    else
    {
        tmp_bounding_box_xy = 0.0;
        tmp_bounding_box_xz = 0.0;
        tmp_bounding_box_yz = 0.0;
    }

// szosty rekord - size_one
    file_.read((char*)&size_one, sizeof(int));
    if ( file_.eof() == 1 )
        return 2;

// wymagamy aby dostepne byly: id type x y z
    if ( size_one < 8 )
        return 2;

// siodmy rekord - nchunk
    file_.read((char*)&nchunk, sizeof(int));
    if ( file_.eof() == 1 )
        return 2;

//testing adasdsdasd
    double** n_atoms = new double*[number_of_atoms1];
    double** n_forces = new double*[number_of_atoms1];

// glowny blok, sklada sie z nchunk blokow zapisanych przez poszczegolne procesory
    for (int chunk = 0; chunk < nchunk; chunk++)
    {
// rozmiar bloku
        file_.read((char*)&n, sizeof(int));
        if ( file_.eof() == 1 )
            return 2;

        if ( n > maxbuf )
        {
            if ( buf != NULL )
                delete [] buf;
            maxbuf = n;
            buf = new double [maxbuf];
        }

        file_.read((char*)buf, n * sizeof(double));
        if ( file_.eof() == 1 )
            return 2;

            n /= size_one;
            int i = 0;
            for (int j = 0; j < n; j++)
            {
                atom_number = int(buf[i++]);
                atom_type = int(buf[i++]);
                atom_name_ss.clear();
                atom_name_ss.str("");
                atom_name_ss << atom_type;
                atom_name = atom_name_ss.str();

                n_atoms[atom_number-1] = new double[3];
                n_atoms[atom_number-1][0] = buf[i++];
                n_atoms[atom_number-1][1] = buf[i++];
                n_atoms[atom_number-1][2] = buf[i++];
                n_forces[atom_number-1] = new double[3];
                n_forces[atom_number-1][0] = buf[i++];
                n_forces[atom_number-1][1] = buf[i++];
                n_forces[atom_number-1][2] = buf[i++];
// w przypadku gdy potrzebne inne dane to nalezy parsowac kolejne elementy
                //for (int k = 8; k < size_one; k++)
                //    i++;
            }
    }

    if ( buf != NULL )
        delete [] buf;
//}

  int abn = 0;
  for (int i = 0; i < number_of_atoms1; i++)
  {
      if (n_atoms[i][0] > tmp_bounding_box_xhi)
      {
          n_atoms[i][0] = n_atoms[i][0] - (tmp_bounding_box_xhi - tmp_bounding_box_xlo);
          abn++;
      }
      else if (n_atoms[i][0] < tmp_bounding_box_xlo)
      {
          n_atoms[i][0] = n_atoms[i][0] + (tmp_bounding_box_xhi - tmp_bounding_box_xlo);
          abn++;
      }

      if (n_atoms[i][1] > tmp_bounding_box_yhi)
      {
          n_atoms[i][1] = n_atoms[i][1] - (tmp_bounding_box_yhi - tmp_bounding_box_ylo);
          abn++;
      }
      else if (n_atoms[i][1] < tmp_bounding_box_ylo)
      {
          n_atoms[i][1] = n_atoms[i][1] + (tmp_bounding_box_yhi - tmp_bounding_box_ylo);
          abn++;
      }

      if (n_atoms[i][2] > tmp_bounding_box_zhi)
      {
          n_atoms[i][2] = n_atoms[i][2] - (tmp_bounding_box_zhi - tmp_bounding_box_zlo);
          abn++;
      }
      else if (n_atoms[i][2] < tmp_bounding_box_zlo)
      {
          n_atoms[i][2] = n_atoms[i][2] + (tmp_bounding_box_zhi - tmp_bounding_box_zlo);
          abn++;
      }
  }

  int** near_neigh; // n_list

  vec3d** dist_mat = new vec3d*[number_of_atoms1];
  vec3d** nn_dist = new vec3d*[number_of_atoms1];

  for (int i1 = 0; i1 < number_of_atoms1; i1++)
  {
      dist_mat[i1] = new vec3d [number_of_atoms1];
      for (int j1 = 0; j1 < number_of_atoms1; j1++)
      {
          dist_mat[i1][j1].x = n_atoms[j1][0] - n_atoms[i1][0];
          dist_mat[i1][j1].y = n_atoms[j1][1] - n_atoms[i1][1];
          dist_mat[i1][j1].z = n_atoms[j1][2] - n_atoms[i1][2];
          dist_mat[i1][j1].r = sqrt((n_atoms[j1][0] - n_atoms[i1][0]) * (n_atoms[j1][0] - n_atoms[i1][0]) +
                                  (n_atoms[j1][1] - n_atoms[i1][1]) * (n_atoms[j1][1] - n_atoms[i1][1]) +
                                  (n_atoms[j1][2] - n_atoms[i1][2]) * (n_atoms[j1][2] - n_atoms[i1][2]));
      }
  }
  dist_mat = bd_imp(number_of_atoms1, "pp", "pp", "ff", tmp_bounding_box_xlo, tmp_bounding_box_xhi
                    , tmp_bounding_box_ylo, tmp_bounding_box_yhi, tmp_bounding_box_zlo
                    , tmp_bounding_box_zhi, dist_mat);


  near_neigh = NN1_implementation(number_of_atoms1, dist_mat, cutoff);

    for (int i1 = 0; i1 < number_of_atoms1; i1++)
    {
        nn_dist[i1] = new vec3d [near_neigh[i1][0]];
        for (int j1 = 0; j1 < near_neigh[i1][0]; j1++)
        {
            nn_dist[i1][j1].x = dist_mat[i1][near_neigh[i1][j1 + 1]].x;
            nn_dist[i1][j1].y = dist_mat[i1][near_neigh[i1][j1 + 1]].y;
            nn_dist[i1][j1].z = dist_mat[i1][near_neigh[i1][j1 + 1]].z;
            nn_dist[i1][j1].r = dist_mat[i1][near_neigh[i1][j1 + 1]].r;
        }
    }

    //define_system(n_nums, near_neigh , nn_dist);

        n_fs = n_forces;
        number_of_atoms = number_of_atoms1;
        n_list = near_neigh;
        n_bonds = nn_dist;
        compute_force();
        compute_energy();

    for (int i1 = 0; i1 < n; i1++)
    {
        delete [] dist_mat[i1];
        delete [] nn_dist[i1];
        delete [] near_neigh[i1];
        delete [] n_atoms[i1];
        delete [] n_forces[i1];
    }

    delete [] dist_mat;
    delete [] nn_dist;
    delete [] near_neigh;
    delete [] n_atoms;
    delete [] n_forces;

    return 0;
}
