#include <iostream>
#include <vector>
#include "water.hpp"
#include <cmath>
using namespace std;

double Deg_to_rad(double deg)
{
    double rad = deg * 3.1415926535897932 / 180;
    return rad;
}

double cal_size(int Num) // Calculate the size of the box
{
    double V;
    double N = Num;
    V = 18 * N / 602;
    return V; // The unit of V is (nm^3)
}

double cal_deviation_molecule(vector<water_molecule> sample, int Num) // This function calculates the deviation of temperature
{
    double v_square_sum = 0;
    for (int serial = 0; serial < Num; serial++)
    {
        double a = pow(sample[serial].speed, 2.0);
        v_square_sum += a;
    }
    double v_square_avg = v_square_sum / (Num * 1.0);
    double T_cal = 0.018 * v_square_avg / (3 * 8.314);
    double factor = T_cal / 277.0;
    return factor;
}

double distance_calculator(double a, double b, double c, double d) // This function calculates the distances between two atoms
{
    double result = b - a;
    if (b - a > c)
        result = b - a - d;
    else if (b - a < -c)
        result = b - a + d;
    return result;
}

// The unit of acceleration is m/(s*fs), and the unit of double "d" is nanometer
// These functions calculates the accelerations between atoms
double O_O_acceleration(double d, double R_cut)
{
    double a;
    if (d <= R_cut && d >= -R_cut)
        a = 5.1812 / (d * d) + 1.7564 * pow(10.0, -6) / pow(d, 13) - 8.7278 * pow(10.0, -4) / pow(d, 7);
    else
        a = 0;
    a *= -1;
    return a;
}

double O_H_acceleration(double d, double R_cut)
{
    double a;
    if (d <= R_cut && d >= -R_cut)
        a = 2.5906 / (d * d);
    else
        a = 0;
    return a;
}

double H_H_acceleration(double d, double R_cut)
{
    double a;
    if (d <= R_cut && d >= -R_cut)
        a = -1.296 / (d * d);
    else
        a = 0;
    return a;
}

// This functions calculate the potential energy between atoms
// PS: the unit of energy is J/mol
double O_O_U(double d, double R_cut)
{
    double u;
    if (d <= R_cut && d >= -R_cut)
        u = 93262 / d + 0.00263467 / pow(d, 12) - 2.618343 / pow(d, 6);
    else
        u = 0;
    return u;
}

double O_H_U(double d, double R_cut)
{
    double u;
    if (d <= R_cut && d >= -R_cut)
        u = -46631 / d;
    else
        u = 0;
    return u;
}

double H_H_U(double d, double R_cut)
{
    double u;
    if (d <= R_cut && d >= -R_cut)
        u = 23316 / d;
    else
        u = 0;
    return u;
}

void UpdateAngularMomentum(int a, matrix b, matrix c, matrix d)
{
    matrix aH1, aH2;
    aH1.row = 3, aH2.row = 3, aH1.column = 1, aH2.column = 1;
    for (int k = 1; k <= 3; k++)
    {
        vector<double> a1, a2;
        a1.push_back(water_system[a].Acceleration[1][k - 1]);
        a2.push_back(water_system[a].Acceleration[2][k - 1]);
        aH1.mat.push_back(a1);
        aH2.mat.push_back(a2);
        a1.clear();
        a2.clear();
    }
    matrix Angular_momentum_principal;
    // Principal angular momentum
    water_system[a].torque = Matrix_addition(matrix_times_number(vector_cross_product(water_system[a].rOH1, aH1), 18), matrix_times_number(vector_cross_product(water_system[a].rOH2, aH2), 18));
    // M=rxF
    water_system[a].Angular_momentum_space = Matrix_addition(water_system[a].Angular_momentum_space, matrix_times_number(water_system[a].torque, 0.5));
    // Ls(t)=Ls(t-δt/2)+M*δt/2
    Angular_momentum_principal = Matrix_multiplication(water_system[a].Rotation_Matrix, water_system[a].Angular_momentum_space);
    // Lp=A*Ls
    water_system[a].Angular_velocity = matrix_times_number(Matrix_multiplication(d, Angular_momentum_principal), pow(10.0, -6.0));
    //ω=I^(-1)*Lp
    matrix dQ_t = dQ(water_system[a].Q, water_system[a].Angular_velocity);
    // dQ=dq*ω/2
    matrix Q_t_middle = Matrix_addition(water_system[a].Q, matrix_times_number(dQ_t, 0.5)); // The quarternion at (t+δt/2)
    water_system[a].Angular_momentum_space = Matrix_addition(water_system[a].Angular_momentum_space, matrix_times_number(water_system[a].torque, 0.5));
    // The space angular momentum at (t+δt/2)
    water_system[a].Rotation_Matrix = ProduceRotationMatrix(Q_t_middle); // Rotation matrix at (t+δt/2)
    Angular_momentum_principal = Matrix_multiplication(water_system[a].Rotation_Matrix, water_system[a].Angular_momentum_space);
    // Principal angular momentum at (t+δt/2)
    water_system[a].Angular_velocity = matrix_times_number(Matrix_multiplication(d, Angular_momentum_principal), pow(10.0, -6.0));
    // Angular speed at (t+δt/2)
    matrix dQ_tplushalf = dQ(Q_t_middle, water_system[a].Angular_velocity);
    // dQ at (t+δt/2)
    water_system[a].Q = Matrix_addition(water_system[a].Q, matrix_times_number(dQ_tplushalf, 1));
    // Quarternion at (t+δt)
    water_system[a].Rotation_Matrix = ProduceRotationMatrix(water_system[a].Q);
    // Rotation Matrix at (t+δt)
    water_system[a].rOH1 = Matrix_multiplication(Matrix_Transformation(water_system[a].Rotation_Matrix), b);
    water_system[a].rOH2 = Matrix_multiplication(Matrix_Transformation(water_system[a].Rotation_Matrix), c);
}