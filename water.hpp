#include <vector>
#include "matrix.hpp"
using namespace std;

struct water_molecule
{
    vector<vector<double>> coordinates;
    double FormerCoor[3][3];
    // The 1st line is O's coordinates, the 2nd line H1 and the 3rd line H2, the first element each line is the relative atomic mass
    // PS: The length unit is nanometer
    vector<double> velocity_xyz; // The speed unit is meter per second, and this is set as v(t-δt/2)
    double speed;
    matrix Angular_momentum_space;       // This is set as L(t-δt/2)
    vector<vector<double>> Acceleration; // Unit: m*(s*fs)^-1
    vector<double> Acceleration_total;
    matrix Q;
    matrix Rotation_Matrix;
    matrix torque; // It's actually torque/m
    matrix Angular_velocity;
    matrix rOH1, rOH2;
};

vector<water_molecule> water_system;