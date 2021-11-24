#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <time.h>
#include <random>
#include <string>
#include "functions.hpp"
#include <iomanip>
using namespace std;

int main(int argc, char **argv)
{
    clock_t StartTime, EndTime;
    StartTime = clock();
    // File input part
    string str = argv[2];
    int NumOfThreads;
    istringstream iss(str);
    iss >> NumOfThreads;
    omp_set_num_threads(NumOfThreads);
    fstream cin(argv[3], ios::in);
    string NameOfInputFile = argv[3];
    string NameOfOutputFile;
    string dot = ".";
    for (int i = 0; i < NameOfInputFile.length(); i++)
    {
        string temp(1, NameOfInputFile[i]);
        if (temp.compare(dot) == 0)
            break;
        else
            NameOfOutputFile += temp;
    }
    string TrajectoryFile = NameOfOutputFile + "_Trj.xyz";
    NameOfOutputFile += ".out";
    fstream fileout(NameOfOutputFile, ios::out);
    fstream fout(TrajectoryFile, ios::out);

    // Parameters and statements
    fileout << "This program is based on rigid water molecule model and SPC force field." << endl;
    fileout << "Steepest descents method is used to optimize the initial configuration" << endl;
    fileout << "The bond length is set as 0.1nm, and the bond angle is 109.47 degree." << endl;
    fileout << "The temperature is set as 277K." << endl;
    fileout << "The integral step-size is set as 1 femtosecond." << endl;
    fileout << "Note: The centre of mass is set as the oxygen atom." << endl;
    for (int i = 0; i < 5; i++)
    {
        fileout << endl;
    }
    int N, N_steps;
    double temperature = 277, Pressure = 1000000;
    double Rcut = 1;
    cin >> N;
    cin >> N_steps;

    // System Initialization
    // This part will initialize parameters other than velocity
    // fileout << "Now the program is initializing the system" << endl;
    for (int i = 0; i < 5; i++)
    {
        fileout << endl;
    }
    double Volume, L;                // L is the length of the edge of box, and I dun want to explain Volume
    Volume = cal_size(N);            // Unit: nm^3
    L = pow(Volume, 1.0 / 3);        // Unir: nm
    int xy = floor(pow(N, 1.0 / 3)); // Variable xy represents the number of cells on x and y directions
    if (pow(xy, 3) < N)
    {
        xy++;
    }
    double cell_length = L / xy; // Unit: nm
    int count_x = 1, count_y = 1, count_z = 1;
    srand((unsigned int)time(0));
    matrix ROH1, ROH2; // The vector from O to H in principal axis
    ROH1.row = 3, ROH1.column = 1, ROH2.row = 3, ROH2.column = 1;
    vector<double> Rrow;
    Rrow.push_back(0.0);
    ROH1.mat.push_back(Rrow), ROH2.mat.push_back(Rrow);
    Rrow[0] = -0.081649;
    ROH1.mat.push_back(Rrow);
    Rrow[0] *= -1;
    ROH2.mat.push_back(Rrow);
    Rrow[0] = 0.057736;
    ROH1.mat.push_back(Rrow), ROH2.mat.push_back(Rrow);
    Rrow.clear();
    for (int i = 1; i <= N; i++) // This loop will initialize coordinates, acceleration
    {
        water_molecule temp;
        temp.Q.row = 4, temp.Q.column = 1;
        temp.Angular_momentum_space.row = 3, temp.Angular_momentum_space.column = 1;
        temp.rOH1.column = 1, temp.rOH1.row = 3, temp.rOH2.column = 1, temp.rOH2.row = 3;
        vector<double> O_xyz, H1_xyz, H2_xyz; // Unit: nm
        O_xyz.push_back(16);
        H1_xyz.push_back(1);
        H2_xyz.push_back(1);
        O_xyz.push_back(count_x * cell_length - cell_length / 2);
        O_xyz.push_back(count_y * cell_length - cell_length / 2);
        O_xyz.push_back(count_z * cell_length - cell_length / 2);
        temp.coordinates.push_back(O_xyz);
        int rz1, rx, rz2;
        rz1 = rand() % 180;
        rx = rand() % 180;
        rz2 = rand() % 180;
        vector<double> q0, q1, q2, q3;
        q0.push_back(cos(Deg_to_rad(rx / 2)) * cos(0.5 * Deg_to_rad(rz1 + rz2)));
        q1.push_back(sin(Deg_to_rad(rx / 2)) * cos(0.5 * Deg_to_rad(rz1 - rz2)));
        q2.push_back(sin(Deg_to_rad(rx / 2)) * sin(0.5 * Deg_to_rad(rz1 - rz2)));
        q3.push_back(cos(Deg_to_rad(rx / 2)) * sin(0.5 * Deg_to_rad(rz1 + rz2)));
        temp.Q.mat.push_back(q0);
        temp.Q.mat.push_back(q1);
        temp.Q.mat.push_back(q2);
        temp.Q.mat.push_back(q3);
        temp.Rotation_Matrix = ProduceRotationMatrix(temp.Q);
        temp.rOH1 = Matrix_multiplication(Matrix_Transformation(temp.Rotation_Matrix), ROH1);
        temp.rOH2 = Matrix_multiplication(Matrix_Transformation(temp.Rotation_Matrix), ROH2);
        for (int j = 1; j <= 3; j++)
        {
            H1_xyz.push_back(temp.coordinates[0][j] + temp.rOH1.mat[j - 1][0]);
            H2_xyz.push_back(temp.coordinates[0][j] + temp.rOH2.mat[j - 1][0]);
        }
        temp.coordinates.push_back(H1_xyz);
        temp.coordinates.push_back(H2_xyz);
        vector<double> temp_accel;
        for (int j = 1; j <= 3; j++)
        {
            temp_accel.push_back(0.0);
        }
        vector<double> temp_AMs;
        temp_AMs.push_back(0.0);
        for (int j = 1; j <= 3; j++)
        {
            temp.Acceleration.push_back(temp_accel);
            temp.Angular_momentum_space.mat.push_back(temp_AMs);
            temp.Acceleration_total.push_back(0.0);
        }
        water_system.push_back(temp);
        count_x++;
        if (count_x > xy)
        {
            count_y++;
            count_x = 1;
        }
        if (count_y > xy)
        {
            count_z++;
            count_y = 1;
        }
    }
    // Optimization with steepest descents
    int OptSwitch = 0;
    double LimitForce = 55;
    int steps = 0;
    double A_System = 0;
    double FormerMaxForce = 0;
    double OptCoefficient = 1;
    for (int j = 0; j < N; j++) // Calculate acceleration and potential energy
    {
        for (int k = j + 1; k < N; k++)
        {
            for (int m = 0; m < 3; m++)
            {
                for (int n = 0; n < 3; n++)
                {
                    double distance_x = distance_calculator(water_system[j].coordinates[m][1], water_system[k].coordinates[n][1], Rcut, L);
                    double distance_y = distance_calculator(water_system[j].coordinates[m][2], water_system[k].coordinates[n][2], Rcut, L);
                    double distance_z = distance_calculator(water_system[j].coordinates[m][3], water_system[k].coordinates[n][3], Rcut, L);
                    double distance_j_k = sqrt(distance_x * distance_x + distance_y * distance_y + distance_z * distance_z);
                    double total_mass = water_system[j].coordinates[m][0] + water_system[k].coordinates[n][0];
                    double ajk, ajk_x, ajk_y, ajk_z;
                    if (total_mass == 32.0)
                    {
                        ajk = O_O_acceleration(distance_j_k, Rcut);
                    }
                    else if (total_mass == 17.0)
                    {
                        ajk = O_H_acceleration(distance_j_k, Rcut);
                    }
                    else if (total_mass == 2.0)
                    {
                        ajk = H_H_acceleration(distance_j_k, Rcut);
                    }
                    ajk_x = ajk * distance_x / distance_j_k;
                    ajk_y = ajk * distance_y / distance_j_k;
                    ajk_z = ajk * distance_z / distance_j_k;
                    water_system[j].Acceleration[m][0] += ajk_x;
                    water_system[j].Acceleration[m][1] += ajk_y;
                    water_system[j].Acceleration[m][2] += ajk_z;
                    water_system[k].Acceleration[n][0] -= ajk_x;
                    water_system[k].Acceleration[n][1] -= ajk_y;
                    water_system[k].Acceleration[n][2] -= ajk_z;
                }
            }
        }
    }
    for (int i = 0; i < N; i++)
    {
        double a_temp = 0;
        for (int j = 0; j < 3; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                water_system[i].Acceleration_total[j] = water_system[i].Acceleration[k][j];
            }
            a_temp += pow(water_system[i].Acceleration_total[j], 2);
            if (sqrt(a_temp) > FormerMaxForce)
            {
                FormerMaxForce = sqrt(a_temp);
            }
        }
    }
    while (OptSwitch == 0)
    {
        steps++;
        for (int i = 0; i < N; i++)
        {
            double a_scalar = sqrt(pow(water_system[i].Acceleration_total[0], 2.0) + pow(water_system[i].Acceleration_total[1], 2.0) + pow(water_system[i].Acceleration_total[2], 2.0));
            for (int j = 0; j < 3; j++)
            {
                for (int k = 1; k <= 3; k++)
                {
                    water_system[i].FormerCoor[j][k - 1] = water_system[i].coordinates[j][k];
                    water_system[i].coordinates[j][k] += 0.01 * OptCoefficient * water_system[i].Acceleration_total[k - 1] / a_scalar;
                }
            }
        }
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    water_system[i].Acceleration[j][k] = 0;
                }
                water_system[i].Acceleration_total[j] = 0;
            }
        }
        double CurrentMaxForce = 0;
        for (int j = 0; j < N; j++) // Calculate acceleration and potential energy
        {
            for (int k = j + 1; k < N; k++)
            {
                for (int m = 0; m < 3; m++)
                {
                    for (int n = 0; n < 3; n++)
                    {
                        double distance_x = distance_calculator(water_system[j].coordinates[m][1], water_system[k].coordinates[n][1], Rcut, L);
                        double distance_y = distance_calculator(water_system[j].coordinates[m][2], water_system[k].coordinates[n][2], Rcut, L);
                        double distance_z = distance_calculator(water_system[j].coordinates[m][3], water_system[k].coordinates[n][3], Rcut, L);
                        double distance_j_k = sqrt(distance_x * distance_x + distance_y * distance_y + distance_z * distance_z);
                        double total_mass = water_system[j].coordinates[m][0] + water_system[k].coordinates[n][0];
                        double ajk, ajk_x, ajk_y, ajk_z;
                        if (total_mass == 32.0)
                        {
                            ajk = O_O_acceleration(distance_j_k, Rcut);
                        }
                        else if (total_mass == 17.0)
                        {
                            ajk = O_H_acceleration(distance_j_k, Rcut);
                        }
                        else if (total_mass == 2.0)
                        {
                            ajk = H_H_acceleration(distance_j_k, Rcut);
                        }
                        ajk_x = ajk * distance_x / distance_j_k;
                        ajk_y = ajk * distance_y / distance_j_k;
                        ajk_z = ajk * distance_z / distance_j_k;
                        water_system[j].Acceleration[m][0] += ajk_x;
                        water_system[j].Acceleration[m][1] += ajk_y;
                        water_system[j].Acceleration[m][2] += ajk_z;
                        water_system[k].Acceleration[n][0] -= ajk_x;
                        water_system[k].Acceleration[n][1] -= ajk_y;
                        water_system[k].Acceleration[n][2] -= ajk_z;
                    }
                }
            }
        }
        for (int i = 0; i < N; i++)
        {
            double a_temp = 0;
            for (int j = 0; j < 3; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    water_system[i].Acceleration_total[j] = water_system[i].Acceleration[k][j];
                }
                a_temp += pow(water_system[i].Acceleration_total[j], 2);
                if (sqrt(a_temp) > CurrentMaxForce)
                {
                    CurrentMaxForce = sqrt(a_temp);
                }
            }
        }
        if (CurrentMaxForce < FormerMaxForce)
        {
            FormerMaxForce = CurrentMaxForce;
            if (steps == 10000 || CurrentMaxForce < LimitForce)
            {
                OptSwitch = 1;
                break;
            }
        }
        else if (CurrentMaxForce >= FormerMaxForce)
        {
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    for (int k = 0; k < 3; k++)
                    {
                        water_system[i].coordinates[j][k + 1] = water_system[i].FormerCoor[j][k];
                    }
                }
            }
            if (OptCoefficient == 1)
            {
                OptCoefficient = 0.1;
            }
            else if (OptCoefficient == 0.1)
            {
                OptCoefficient = 0.01;
            }
            else if (OptCoefficient == 0.01)
            {
                OptSwitch = 1;
                break;
            }
        }
    }

    matrix moi_tensor_space_inversed; // moi stands for moment of inertia
    moi_tensor_space_inversed.column = 3, moi_tensor_space_inversed.row = 3;
    for (int i = 1; i <= 3; i++)
    {
        vector<double> moi_row;
        for (int j = 1; j <= 3; j++)
        {
            if (i == j)
            {
                if (i == 1)
                    moi_row.push_back(50);
                else if (i == 2)
                    moi_row.push_back(149.9948);
                else if (i == 3)
                    moi_row.push_back(75.0019);
            }
            else
            {
                moi_row.push_back(0.0);
            }
        }
        moi_tensor_space_inversed.mat.push_back(moi_row);
        moi_row.clear();
    }

    // This part will initialize the velocity part of the system
    double total_vx = 0.0, total_vy = 0.0, total_vz = 0.0;
    default_random_engine generator((unsigned int)(time(0)));
    double V_bound = 357.69; // 357.69 m/s is the root-mean-square speed of each direction under 277K
    normal_distribution<double> MaxwellDistribution(0, 1);
    for (int i = 0; i < N; i++)
    {
        double vx = 0, vy = 0, vz = 0;
        if (i == N - 1)
        {
            vx = 0 - total_vx;
            vy = 0 - total_vy;
            vz = 0 - total_vz;
        }
        else
        {
            double coefficient_x = 2, coefficient_y = 2, coefficient_z = 2;
            while (coefficient_x > 1 || coefficient_x < -1)
            {
                coefficient_x = MaxwellDistribution(generator);
            }
            while (coefficient_y > 1 || coefficient_y < -1)
            {
                coefficient_y = MaxwellDistribution(generator);
            }
            while (coefficient_z > 1 || coefficient_z < -1)
            {
                coefficient_z = MaxwellDistribution(generator);
            }
            vx = coefficient_x * V_bound;
            vy = coefficient_y * V_bound;
            vz = coefficient_z * V_bound;
            if (total_vx > 0)
                vx = -fabs(vx);
            else if (total_vx < 0)
                vx = fabs(vx);
            if (total_vy > 0)
                vy = -fabs(vy);
            else if (total_vy < 0)
                vy = fabs(vy);
            if (total_vz > 0)
                vz = -fabs(vz);
            else if (total_vz < 0)
                vz = fabs(vz);
            total_vx += vx;
            total_vy += vy;
            total_vz += vz;
        }
        water_system[i].velocity_xyz.push_back(vx);
        water_system[i].velocity_xyz.push_back(vy);
        water_system[i].velocity_xyz.push_back(vz);
        double v_scalar = sqrt(vx * vx + vy * vy + vz * vz);
        water_system[i].speed = v_scalar;
    }
    double f = cal_deviation_molecule(water_system, N); // This part will calibrate the deviation of velocity
    if (f > 1.1 || f < 0.9)
    {
        double a = sqrt(1.0 / f);
        for (int i = 0; i < N; i++)
        {
            water_system[i].velocity_xyz[0] *= a;
            water_system[i].velocity_xyz[1] *= a;
            water_system[i].velocity_xyz[2] *= a;
            water_system[i].speed *= a;
        }
    }

    // Start simulations
    // Rcut = 1nm, switch function is not used
    // Leap Frog Method is used to calculate the coordinates, velocity and angular momentum of the particle
    fileout << "**********************************************************" << endl;
    fileout << "Starting simulation" << endl;
    fileout.setf(ios::fixed);

    double potential_energy_sum = 0;
    double sum_times = 0;
    fileout << "ΔH(vap) list" << endl;
    fileout << "Simulation Step                   Enthalpy of vaporization" << endl;
    for (int i = 1; i <= N_steps; i++)
    {
        printf("Step %d\n", i);
#pragma omp parallel
        {
#pragma omp for
            for (int j = 0; j < N; j++)
            {
                for (int p = 0; p < 3; p++)
                {
                    water_system[j].Acceleration_total[p] = 0;
                    for (int q = 0; q < 3; q++)
                    {
                        water_system[j].Acceleration[p][q] = 0;
                    }
                }
            }
        }
        double pe_temp = 0;         // pe stands for potential energy
        for (int j = 0; j < N; j++) // Calculate acceleration and potential energy
        {

            for (int k = j + 1; k < N; k++)
            {
                for (int m = 0; m < 3; m++)
                {
                    for (int n = 0; n < 3; n++)
                    {
                        double distance_x = distance_calculator(water_system[j].coordinates[m][1], water_system[k].coordinates[n][1], Rcut, L);
                        double distance_y = distance_calculator(water_system[j].coordinates[m][2], water_system[k].coordinates[n][2], Rcut, L);
                        double distance_z = distance_calculator(water_system[j].coordinates[m][3], water_system[k].coordinates[n][3], Rcut, L);
                        double distance_j_k = sqrt(distance_x * distance_x + distance_y * distance_y + distance_z * distance_z);
                        double total_mass = water_system[j].coordinates[m][0] + water_system[k].coordinates[n][0];
                        double ajk, ajk_x, ajk_y, ajk_z;
                        if (total_mass == 32.0)
                        {
                            ajk = O_O_acceleration(distance_j_k, Rcut);
                            if (i % 10 == 0) // Print ΔH(vap) every 1000 steps
                            {
                                pe_temp += O_O_U(distance_j_k, Rcut);
                            }
                        }
                        else if (total_mass == 17.0)
                        {
                            ajk = O_H_acceleration(distance_j_k, Rcut);
                            if (i % 10 == 0)
                            {
                                pe_temp += O_H_U(distance_j_k, Rcut);
                            }
                        }
                        else if (total_mass == 2.0)
                        {
                            ajk = H_H_acceleration(distance_j_k, Rcut);
                            if (i % 10 == 0)
                            {
                                pe_temp += H_H_U(distance_j_k, Rcut);
                            }
                        }
                        ajk_x = ajk * distance_x / distance_j_k;
                        ajk_y = ajk * distance_y / distance_j_k;
                        ajk_z = ajk * distance_z / distance_j_k;
                        water_system[j].Acceleration[m][0] += ajk_x;
                        water_system[j].Acceleration[m][1] += ajk_y;
                        water_system[j].Acceleration[m][2] += ajk_z;
                        water_system[k].Acceleration[n][0] -= ajk_x;
                        water_system[k].Acceleration[n][1] -= ajk_y;
                        water_system[k].Acceleration[n][2] -= ajk_z;
                    }
                }
            }
        }
        if (pe_temp != 0)
        {
            pe_temp /= N;
            fileout << i << "           " << -pe_temp << endl;
            if (i >= 10)
            {
                potential_energy_sum += pe_temp;
                sum_times++;
            }
        }
// Modifying the coordinates and all other stuff related to them
#pragma omp parallel
        {
#pragma omp for
            for (int j = 0; j < N; j++)
            {
                matrix aH1, aH2;
                aH1.row = 3, aH2.row = 3, aH1.column = 1, aH2.column = 1;
                for (int k = 1; k <= 3; k++)
                {
                    vector<double> a1, a2;
                    a1.push_back(water_system[j].Acceleration[1][k - 1]);
                    a2.push_back(water_system[j].Acceleration[2][k - 1]);
                    aH1.mat.push_back(a1);
                    aH2.mat.push_back(a2);
                    a1.clear();
                    a2.clear();
                }
                matrix Angular_momentum_principal;
                // Principal angular momentum
                water_system[j].torque =
                    Matrix_addition(matrix_times_number(vector_cross_product(water_system[j].rOH1, aH1), 18), matrix_times_number(vector_cross_product(water_system[j].rOH2, aH2), 18));
                // M=rxF
                water_system[j].Angular_momentum_space = Matrix_addition(water_system[j].Angular_momentum_space, matrix_times_number(water_system[j].torque, 0.5));
                // Ls(t)=Ls(t-δt/2)+M*δt/2
                Angular_momentum_principal = Matrix_multiplication(water_system[j].Rotation_Matrix, water_system[j].Angular_momentum_space);
                // Lp=A*Ls
                water_system[j].Angular_velocity = matrix_times_number(Matrix_multiplication(moi_tensor_space_inversed, Angular_momentum_principal), pow(10.0, -6.0));
                //ω=I^(-1)*Lp
                matrix dQ_t = dQ(water_system[j].Q, water_system[j].Angular_velocity);
                // dQ=dq*ω/2
                matrix Q_t_middle = Matrix_addition(water_system[j].Q, matrix_times_number(dQ_t, 0.5)); // The quarternion at (t+δt/2)
                water_system[j].Angular_momentum_space = Matrix_addition(water_system[j].Angular_momentum_space, matrix_times_number(water_system[j].torque, 0.5));
                // The space angular momentum at (t+δt/2)
                water_system[j].Rotation_Matrix = ProduceRotationMatrix(Q_t_middle); // Rotation matrix at (t+δt/2)
                Angular_momentum_principal = Matrix_multiplication(water_system[j].Rotation_Matrix, water_system[j].Angular_momentum_space);
                // Principal angular momentum at (t+δt/2)
                water_system[j].Angular_velocity = matrix_times_number(Matrix_multiplication(moi_tensor_space_inversed, Angular_momentum_principal), pow(10.0, -6.0));
                // Angular speed at (t+δt/2)
                matrix dQ_tplushalf = dQ(Q_t_middle, water_system[j].Angular_velocity);
                // dQ at (t+δt/2)
                water_system[j].Q = Matrix_addition(water_system[j].Q, matrix_times_number(dQ_tplushalf, 1));
                // Quarternion at (t+δt)
                water_system[j].Rotation_Matrix = ProduceRotationMatrix(water_system[j].Q);
                // Rotation Matrix at (t+δt)
                water_system[j].rOH1 = Matrix_multiplication(Matrix_Transformation(water_system[j].Rotation_Matrix), ROH1);
                water_system[j].rOH2 = Matrix_multiplication(Matrix_Transformation(water_system[j].Rotation_Matrix), ROH2);
            }
        }
#pragma omp parallel
        {
#pragma omp for
            for (int j = 0; j < N; j++)
            {
                for (int p = 0; p < 3; p++)
                {
                    for (int q = 0; q < 3; q++)
                    {
                        water_system[j].Acceleration_total[p] += water_system[j].Acceleration[q][p];
                    }
                }
                double speed_temp = 0;
                for (int k = 0; k < 3; k++) // This loop calculate the velocity of next step
                {
                    water_system[j].velocity_xyz[k] += 1 * water_system[j].Acceleration_total[k];
                    speed_temp += pow(water_system[j].velocity_xyz[k], 2.0);
                }
                water_system[j].speed = sqrt(speed_temp);
            }
        }
        if (i <= 100)
        {
            f = cal_deviation_molecule(water_system, N);
            if (f > 1.1 || f < 0.9)
            {
                double a = sqrt(1.0 / f);
#pragma omp parallel
                {
#pragma omp for
                    for (int j = 0; j < N; j++)
                    {
                        water_system[j].velocity_xyz[0] *= a;
                        water_system[j].velocity_xyz[1] *= a;
                        water_system[j].velocity_xyz[2] *= a;
                        water_system[j].speed *= a;
                    }
                }
            }
        }
        else if (i > 100 && i <= 1000)
        {
            if (i % 100 == 0)
            {
                f = cal_deviation_molecule(water_system, N);
                if (f > 1.1 || f < 0.9)
                {
                    double a = sqrt(1.0 / f);
#pragma omp parallel
                    {
#pragma omp for
                        for (int j = 0; j < N; j++)
                        {
                            water_system[j].velocity_xyz[0] *= a;
                            water_system[j].velocity_xyz[1] *= a;
                            water_system[j].velocity_xyz[2] *= a;
                            water_system[j].speed *= a;
                        }
                    }
                }
            }
        }
        else if (i > 1000)
        {
            if (i % 1000 == 0)
            {
                f = cal_deviation_molecule(water_system, N);
                if (f > 1.1 || f < 0.9)
                {
                    double a = sqrt(1.0 / f);
#pragma omp parallel
                    {
#pragma omp for
                        for (int j = 0; j < N; j++)
                        {
                            water_system[j].velocity_xyz[0] *= a;
                            water_system[j].velocity_xyz[1] *= a;
                            water_system[j].velocity_xyz[2] *= a;
                            water_system[j].speed *= a;
                        }
                    }
                }
            }
        }
#pragma omp parallel
        {
#pragma omp for
            for (int j = 0; j < N; j++)
            {
                for (int k = 1; k <= 3; k++)
                {
                    water_system[j].coordinates[0][k] += 1 * water_system[j].velocity_xyz[k - 1] / 1000000;
                    water_system[j].coordinates[1][k] = water_system[j].coordinates[0][k] + water_system[j].rOH1.mat[k - 1][0];
                    water_system[j].coordinates[2][k] = water_system[j].coordinates[0][k] + water_system[j].rOH2.mat[k - 1][0];
                }
                for (int p = 0; p < 3; p++)
                {
                    for (int q = 1; q <= 3; q++)
                    {
                        if (water_system[j].coordinates[p][q] > L)
                        {
                            water_system[j].coordinates[p][q] -= L;
                        }
                        else if (water_system[j].coordinates[p][q] < 0)
                        {
                            water_system[j].coordinates[p][q] += L;
                        }
                    }
                }
            }
        }
        if (i % 10 == 0)
        {
            fout << "                    " << N * 3 << endl;
            fout << "          Step " << i << endl;
            fout.setf(ios::fixed);
            fout << setprecision(8);
            for (int j = 0; j < N; j++)
            {
                fout << "O             " << water_system[j].coordinates[0][1] * 10 << "        " << water_system[j].coordinates[0][2] * 10 << "        " << water_system[j].coordinates[0][3] * 10
                     << endl;
                fout << "H             " << water_system[j].coordinates[1][1] * 10 << "        " << water_system[j].coordinates[1][2] * 10 << "        " << water_system[j].coordinates[1][3] * 10
                     << endl;
                fout << "H             " << water_system[j].coordinates[2][1] * 10 << "        " << water_system[j].coordinates[2][2] * 10 << "        " << water_system[j].coordinates[2][3] * 10
                     << endl;
            }
        }
    }
    double enthalpy = potential_energy_sum / sum_times;
    fileout << "The average enthalpy of vaporization is: " << -enthalpy << " J/mol" << endl;
    if (enthalpy != 0)
    {
        fileout << "Normal Termination" << endl;
    }
    else
    {
        fileout << "Error Termination" << endl;
    }
    EndTime = clock();
    double PastTime = (double)(EndTime - StartTime) / CLOCKS_PER_SEC;
    PastTime = round(PastTime);
    int hours, minutes, seconds, RemainTime;
    hours = floor(PastTime / 3600);
    RemainTime = PastTime - hours * 3600;
    minutes = RemainTime / 60;
    RemainTime = RemainTime - 60 * minutes;
    fileout << setprecision(0);
    fileout << "Program running duration: " << hours << ":" << minutes << ":" << RemainTime << endl;
    for (int i = 0; i <= 5; i++)
    {
        fileout << endl;
    }
    fileout << "This program is supported by cola and green tea";
    return 0;
}