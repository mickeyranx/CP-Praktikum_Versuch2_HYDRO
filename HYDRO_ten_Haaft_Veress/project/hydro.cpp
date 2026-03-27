

#include <iostream>
#include <vector>
#include <tuple>
#include <iomanip>
#include <fstream>
#include <cstdio>
#include <ctime>
#include <cmath>

using namespace std;

//b-type
static vector<double> psi_prev;
static vector<double> psi_next;
static vector<double> rho_prev;
static vector<double> rho_next;
static vector<double> eps_prev;
static vector<double> eps_next;

//a-type
static vector<double> vel_prev;
static vector<double> vel_next;


//determines the next densities for each cell with upwind
vector<double> next_density(vector<double> &rho_prev, vector<double> &rho_next, vector<double> &vel_prev,int N,double step) {
    
    
    //calculate the first del_rho_j_b using the ghost cells used in the extenden upwind 
    double check_value = (rho_prev[2] - rho_prev[1]) * (rho_prev[1] - rho_prev[0]);
    double del_rho_j_b = check_value > 0 ? 2 * (check_value/(rho_prev[2] - rho_prev[0])) : 0;


    vector<double> F_m(N+4,0); //initialize mass fluxes
    
    //caclulation of the mass fluxes
    for (int j = 2; j < N+3; j++) //loop throug non ghost cells only
    {
        double rho_j = rho_prev[j];
        
        double rho_j_f = rho_prev[j+1]; //f is forward
     
        double rho_j_b = rho_prev[j-1]; //b is backward

        double u_j = vel_prev[j];

        
        bool velocity_forward = vel_prev[j] > 0;

        //calculation of rho_adv_j with extended upwind 
        double check_value = (rho_j_f - rho_j) * (rho_j - rho_j_b); //geometric mean, van Leer
        double del_rho_j = check_value > 0 ? 2 * (check_value / (rho_j_f - rho_j_b)) : 0;
        double rho_adv_j = 0;
        if(velocity_forward){
            rho_adv_j = rho_j_b + 0.5 * (1 - u_j * step) * del_rho_j_b;
        }
        else {
            rho_adv_j = rho_j - 0.5 * (1 + u_j * step) * del_rho_j;
        }

        F_m[j] = rho_adv_j * u_j; //calculate mass flux for cell j
        del_rho_j_b = del_rho_j; //pass on the current delta del_rho_j for the next iteration
        
    }

    for (int j = 2; j < N+2; j++)//calcualte the denities with the previously calculated mass fluxes
    {
        rho_next[j] = rho_prev[j] - step * (F_m[j + 1] - F_m[j]);
    }
    
    //dont need to return the densities because they are changed within this method by reference
    return F_m;

}


//calculates the next velocity profile with advection
vector<double> next_velocity_adv(int N,double step ,vector<double> F_m){
    vector<double> u_1(N + 4, 0); //note that u[0] is not used at all and only exists to prevent confusion
   
    double del_u_j = 0;
    vector<double> F_I(N + 4, 0);

    for (int j = 2; j < N + 2; j++) //loop throug all non ghost cells 
    {

        double u_j_b = vel_prev[j - 1];
        double u_j = vel_prev[j];
        double u_j_f = vel_prev[j + 1];
        double u_j_f_f = vel_prev[j + 2]; //_f_f = 2 times forward
        
        if (j == 2) { //for the first iteration we calculate del_u_j, for the next iterations it is passed in by del_u_j_f from the previous iteration
            double check_value = (u_j_f - u_j) * (u_j - u_j_b);
            del_u_j = check_value > 0 ? 2 * (check_value / (u_j_f - u_j_b)) : 0;
        }

       
        double av_u_j = 0.5 * (u_j + u_j_f); //av = averaged value

        //calculate u_adv_j with extension of upwind
        double check_value = (u_j_f_f - u_j_f) * (u_j_f - u_j);
        double del_u_j_f = check_value > 0 ? 2 * (check_value / (u_j_f_f - u_j)) : 0;
        
        double u_adv_j = 0;
        if (av_u_j > 0) {
            u_adv_j = u_j + 0.5 * (1 - av_u_j * step) * del_u_j;
        }
        else {
            u_adv_j = u_j_f - 0.5 * (1 + av_u_j * step) * del_u_j_f;
        }

        //calculate F_I profile (which is technically not needed, because we will only need the previous and the current F_I)
        double F_I_j = 0.5 * (F_m[j] + F_m[j + 1]) * u_adv_j;
        F_I[j] = F_I_j;

        double rho_av_prev = 0.5 * (rho_prev[j - 1] + rho_prev[j]);
        double rho_av_next = 0.5 * (rho_next[j - 1] + rho_next[j]);

        //calculate u_1
        if (j > 2) { //not even sure if this is needed, dont want to touch it thoug because the code works. u_1[2] is set to 0 due to the boundary conditions anyway
            u_1[j] = (u_j * rho_av_prev - step * (F_I_j - F_I[j - 1])) / rho_av_next;
        }
        del_u_j = del_u_j_f; //pass on the current delta del_u_j_f for the next iteration
    }

    return u_1;

}

//calculates the next energy profile with advection
vector<double> next_energy_adv(int N ,double step, vector<double> F_m) {
    vector<double> e_1(N + 4, 0);
    vector<double> F_e(N+4,0);
    
    double check_value = (eps_prev[2] - eps_prev[1]) * (eps_prev[1] - eps_prev[0]);
    double del_e_j_b = check_value > 0 ? 2 * check_value/(eps_prev[2] - eps_prev[0]) : 0;
    for (int j = 2; j < N + 3; j++)
    {
        double e_j_b = eps_prev[j - 1];
        double e_j = eps_prev[j];
        double e_j_f = eps_prev[j + 1];

        double u_j = vel_prev[j];


        bool velocity_forward = u_j > 0;
        double check_value = (e_j_f - e_j) * (e_j - e_j_b);
        double del_e_j = check_value > 0 ? 2 * (check_value / (e_j_f - e_j_b)) : 0;
  
        double e_adv_j = 0;
        if (velocity_forward) {
            e_adv_j = e_j_b + 0.5 * (1 - u_j * step) * del_e_j_b;
        }
        else {
            e_adv_j = e_j - 0.5 * (1 + u_j * step) * del_e_j;
        }

        F_e[j] = F_m[j] * e_adv_j; //calcualte the energy flux


        del_e_j_b = del_e_j;
    }
    for (int j = 2; j < N+2; j++)
    {

        e_1[j] = (eps_prev[j] * rho_prev[j] - step * (F_e[j + 1] - F_e[j])) / rho_next[j];
    }
    return e_1;
}

void simulate_next_euler_step(int N, double step) {

    //1. advection step
    //boundary conditions for _prev vectors are  already applied through first time or end of this method

    //start with density which also return F_m profile for the calculation of the next step
    //the new density profile is saved in rho_next
    vector<double> F_m = next_density(rho_prev, rho_next, vel_prev,N, step);

    //apply boundary conditions for rho_next
    rho_next[1] = rho_next[2];
    rho_next[0] = rho_next[3];
    rho_next[N + 2] = rho_next[N + 1];
    rho_next[N + 3] = rho_next[N];

    
    vector<double> u_1 = next_velocity_adv(N, step, F_m); //calcualte u_1 profile with advection

    //apply boundary conditions for u_1
    u_1[2] = 0;
    u_1[1] = -u_1[3];
    u_1[N + 2] = 0;
    u_1[N + 3] = -u_1[N + 1];

    vector<double> e_1 = next_energy_adv(N, step , F_m);  //calculate e_1 profile with advection

    //apply boundary conditions fpr e_1
    e_1[1] = e_1[2];
    e_1[0] = e_1[3];
    e_1[N + 2] = e_1[N + 1];
    e_1[N + 3] = e_1[N];

    

    //2. apply work done by pressure
    
    //calculate pressure
    vector<double> pressure(N+4, 0);
    for (int j = 2; j < N+2; j++)
    {
        pressure[j] = 0.4 * rho_next[j] * e_1[j];
        double p_j = pressure[j];
        if (j > 2) {
            vel_next[j] = u_1[j] - 2 * step * (p_j - pressure[j - 1]) / (rho_next[j - 1] + rho_next[j]);
        }
        
        eps_next[j] = e_1[j] - step * p_j / rho_next[j] * (u_1[j + 1] - u_1[j]);
    }

    
    //apply boundary conditions after forces
    vel_next[2] = 0;
    vel_next[1] = -vel_next[3];
    vel_next[N + 2] = 0;
    vel_next[N + 3] = -vel_next[N + 1];

    eps_next[1] = eps_next[2];
    eps_next[0] = eps_next[3];
    eps_next[N + 2] = eps_next[N + 1];
    eps_next[N + 3] = eps_next[N];



    //prepare for the next step
    rho_prev = rho_next;
    vel_prev = vel_next;
    eps_prev = eps_next;

}


//performs exercise 1 input params are: number of lattice points N, simulation time t, velocity u, courant number sigma, calculation domain [xmin, xmax] and the name of the output file
void linearAdvection(int N, double time, double u ,double sigma,int xmin, int xmax, string filename) {
    
    // 1. initialize vectors
    vector<double> psi_prev(N + 4, 0); //N+4 because of 2 ghost cells on each side
    vector<double> psi_next(N + 4, 0);
    vector<double> lattice(N + 4, 0);
    vector<double> vel_prev(N + 4, u);
    vector<double> vel_next(N + 4, u);

    //calculate timestep and step size throug courant number and Delta_x
    double del_x = (xmax - xmin) / (double) N;
    double time_step = sigma * del_x / u;
    double step = sigma / u;

    int k = time/time_step; //number of iterations for the respective time 
    
    for (int i = 0; i < lattice.size(); i++)//create lattice x axis profile
    {
        lattice[i] = -1 - 2 * del_x + i * del_x;
    }
    
    for (int i = 0; i < lattice.size(); i++)//invoke psi_0 (starting profile at t=0)
    {
        if (abs(lattice[i]) <= 1.0 / 3.0) psi_prev[i] = 1.0;
         
    }

    
    for (int i = 0; i < k; i++) //performing iterations
    {
        
        //the densities and velocities are passed in by reference, so the method "next_density" already changes these values
        next_density(psi_prev, psi_next, vel_prev, N, step); //calculate next density profile, the return value ist not important 
        psi_prev = psi_next;
        //boundary conditions for the next iteration
        psi_prev[0] = psi_next[N];
        psi_prev[1] = psi_next[N + 1];
        psi_prev[N + 2] = psi_next[2];
        psi_prev[N + 3] = psi_next[3];
    }
    
    
    //write data to file
    ofstream data_file;
    data_file.open(filename);
    data_file << "x" << "\t" << "psi" << "\n";
    for (int i = 2; i < N + 3; i++)
    {
        data_file << fixed << setprecision(3) << lattice[i] << "\t";
        data_file << fixed << setprecision(6) << psi_prev[i] << "\n";
    }

    data_file.close();
}


//contains exercise 2: 1D-shock-tube
void shockTube(string output_filename) {
    //no. iterations
    int k = 228;

    //starting values
    double rho_l = 1.0;
    double rho_r = 0.125;
    double p_l = 1.0;
    double p_r = 0.1;
    double u_l = 0.0;
    double u_r = 0.0;
    double eps_l = 2.5;
    double eps_r = 2.0;
    double x_0 = 0.5;
    //computing area
    double x_min = 0.0;
    double x_max = 1.0;
    //no. lattice points
    double N = 100;
    double del_x = 0.01;
    double time_step = 0.001;
    double step = time_step / del_x;

    //initialize vectors
    rho_prev = vector<double>(N+4,0);
    rho_next = vector<double>(N + 4, 0);
    vel_prev = vector<double>(N + 4, 0);
    vel_next = vector<double>(N + 4, 0);
    eps_prev = vector<double>(N + 4, 0);
    eps_next = vector<double>(N + 4, 0);
    
    vector<double> lattice = vector<double>(N + 4, 0);

   
    for (int i = 0; i < lattice.size(); i++) //invoke calculating area
    {
        lattice[i] = - 2 * del_x + i * del_x;
    }

    
    for (int i = 0; i < N+4; i++) //apply starting condition (starting profile)
    {
        if (lattice[i] <= x_0) {
            rho_prev[i] = rho_l;
            
            eps_prev[i] = eps_l;
        }
        else {
            rho_prev[i] = rho_r;
            
            eps_prev[i] = eps_r;
        }

    }
    //apply boundary conditions for starting conditions
    eps_prev[1] = eps_prev[2];
    eps_prev[0] = eps_prev[3];
    eps_prev[N + 2] = eps_prev[N + 1];
    eps_prev[N + 3] = eps_prev[N];

    rho_prev[1] = rho_prev[2];
    rho_prev[0] = rho_prev[3];
    rho_prev[N + 2] = rho_prev[N + 1];
    rho_prev[N + 3] = rho_prev[N];
    
    //simulation k steps
    for (int i = 0; i < k; i++)
    {
        simulate_next_euler_step(N, step);
    }


    //write data to file
    ofstream data_file;
    data_file.open(output_filename);
    data_file << "x_B" << "\t" << "rho" << "\t" << "eps" << "\t" << "p" << "\t" << "T" << "\t" << "C" << "\t";
    data_file << "x_A" << "\t" << "u" << "\n";

    for (int i = 2; i < N + 2; i++) //write all data of the cells and other computable quantieties to a file
    {
        double eps = eps_prev[i];
        double rho = rho_prev[i];
        double vel = vel_prev[i];
        //a-type
        
        data_file << fixed << setprecision(5);
        data_file << lattice[i] + 0.5 * del_x << "\t";
        data_file << rho << "\t";
        data_file << eps << "\t";
        //pressure
        double p_j = 0.4 * eps * rho; //the 0.4 comes from 1 - gamma, with gamma being the adiabat index of 1.4 in this experiment
        data_file << p_j << "\t";
        //temperature
        data_file << 0.4 * eps << "\t";
        //courant-number
        data_file << (sqrt(1.4 * p_j / rho) + abs(vel)) * step << "\t";
        data_file << lattice[i] << "\t" << vel << "\n";
    }

    data_file.close();

}

int main()
{
    clock_t start;
    double runtime;
    start = clock();
    //--------------------------------------------
    //                 exercise 1
    //--------------------------------------------
    cout << "doing exercise 1 : linear advection" << endl;
    //filename of generated data file
    string filename_1 = "results/exercise_1_N=40_t=4_large_courant.txt";
    string filename_2 = "results/exercise_1_N=40_t=4.txt";
    string filename_3 = "results/exercise_1_N=400_t=400.txt";
    int xmin = -1;
    int xmax = 1;
    double velocity = 1.0;
    linearAdvection(40, 4.0, velocity, 1.2, xmin, xmax, filename_1);
    linearAdvection(40, 4.0, velocity, 0.8, xmin, xmax, filename_2);
    linearAdvection(400, 400.0, velocity, 0.8, xmin, xmax, filename_3);
   
    cout << "exercise 1 done " << endl;
    //--------------------------------------------
    //                 exercise 2
    //--------------------------------------------
    
    cout << "doing exercise 2 : 1D-Shocktube " << endl;
    
    //filename of generated data file
    string filename = "results/exercise_2_t=.228";
    //inputparams can be changed inside the method
    shockTube(filename);
    
    
    cout << "exercise 2 done " << endl;
    
    runtime = (clock() - start) / (double) CLOCKS_PER_SEC;
   
    cout << "runtime: " << runtime << "s" << "\n";
    cout << "code ran succesfully";
    
}

