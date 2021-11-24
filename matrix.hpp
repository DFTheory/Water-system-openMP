#include <iostream>
#include <vector>
#include <omp.h>
using namespace std;

struct matrix
{
    vector<vector<double>> mat;
    int row;
    int column;
};

matrix Matrix_multiplication(matrix a, matrix b)
{
    if (a.column != b.row)
    {
        cout << "Wrong Input for multiplication!";
        abort();
    }
    else
    {
        matrix temp;
        temp.row = a.row;
        temp.column = b.column;
        vector<double> row_temp;
        for (int i = 0; i < a.row; i++)
        {
            for (int j = 0; j < b.column; j++)
            {
                double matrix_element = 0;
                for (int p = 0; p < a.column; p++)
                {
                    matrix_element += a.mat[i][p] * b.mat[p][j];
                }
                row_temp.push_back(matrix_element);
            }
            temp.mat.push_back(row_temp);
            row_temp.clear();
        }
        return temp;
    }
}

matrix Matrix_Transformation(matrix a)
{
    matrix transformed;
    transformed.row = a.column;
    transformed.column = a.row;
    vector<double> temp;
    for (int i = 0; i < a.column; i++)
    {
        for (int j = 0; j < a.row; j++)
        {
            temp.push_back(a.mat[j][i]);
        }
        transformed.mat.push_back(temp);
        temp.clear();
    }
    return transformed;
}

void print_matrix(matrix a)
{
    for (int i = 0; i < a.row; i++)
    {
        for (int j = 0; j < a.column; j++)
        {
            cout << a.mat[i][j] << " ";
        }
        cout << endl;
    }
}

matrix vector_cross_product(matrix a, matrix b) // a and b are column vectors
{
    matrix cross_product;
    cross_product.column = 1, cross_product.row = 3;
    vector<double> temp;
    temp.push_back(a.mat[1][0] * b.mat[2][0] - a.mat[2][0] * b.mat[1][0]);
    cross_product.mat.push_back(temp);
    temp.clear();
    temp.push_back(a.mat[2][0] * b.mat[0][0] - a.mat[0][0] * b.mat[2][0]);
    cross_product.mat.push_back(temp);
    temp.clear();
    temp.push_back(a.mat[0][0] * b.mat[1][0] - a.mat[1][0] * b.mat[0][0]);
    cross_product.mat.push_back(temp);
    temp.clear();
    return cross_product;
}

matrix Matrix_addition(matrix a, matrix b)
{
    matrix sum_of_matrix;
    sum_of_matrix.row = a.row;
    sum_of_matrix.column = a.column;
    if (a.row != b.row || a.column != b.column)
    {
        cout << "Wrong Input for addition!";
        abort();
    }
    else
    {
        vector<double> row_temp;
        for (int i = 0; i < a.row; i++)
        {
            for (int j = 0; j < a.column; j++)
            {
                row_temp.push_back(a.mat[i][j] + b.mat[i][j]);
            }
            sum_of_matrix.mat.push_back(row_temp);
            row_temp.clear();
        }
        return sum_of_matrix;
    }
}

matrix matrix_times_number(matrix a, double b)
{
    matrix c;
    c.row = a.row, c.column = a.column;
    for (int i = 0; i < a.row; i++)
    {
        vector<double> row_temp;
        for (int j = 0; j < a.column; j++)
        {
            row_temp.push_back(b * a.mat[i][j]);
        }
        c.mat.push_back(row_temp);
        row_temp.clear();
    }
    return c;
}

matrix ProduceRotationMatrix(matrix q)
{
    matrix rotation;
    rotation.row = 3, rotation.column = 3;
    vector<double> rotation_row;
    rotation_row.push_back(q.mat[0][0] * q.mat[0][0] + q.mat[1][0] * q.mat[1][0] - q.mat[2][0] * q.mat[2][0] - q.mat[3][0] * q.mat[3][0]);
    rotation_row.push_back(2 * (q.mat[1][0] * q.mat[2][0] + q.mat[0][0] * q.mat[3][0]));
    rotation_row.push_back(2 * (q.mat[1][0] * q.mat[3][0] - q.mat[0][0] * q.mat[2][0]));
    rotation.mat.push_back(rotation_row);
    rotation_row.clear();
    rotation_row.push_back(2 * (q.mat[1][0] * q.mat[2][0] - q.mat[0][0] * q.mat[3][0]));
    rotation_row.push_back(q.mat[0][0] * q.mat[0][0] - q.mat[1][0] * q.mat[1][0] + q.mat[2][0] * q.mat[2][0] - q.mat[3][0] * q.mat[3][0]);
    rotation_row.push_back(2 * (q.mat[2][0] * q.mat[3][0] + q.mat[0][0] * q.mat[1][0]));
    rotation.mat.push_back(rotation_row);
    rotation_row.clear();
    rotation_row.push_back(2 * (q.mat[1][0] * q.mat[3][0] + q.mat[0][0] * q.mat[2][0]));
    rotation_row.push_back(2 * (q.mat[2][0] * q.mat[3][0] - q.mat[0][0] * q.mat[1][0]));
    rotation_row.push_back(q.mat[0][0] * q.mat[0][0] - q.mat[1][0] * q.mat[1][0] - q.mat[2][0] * q.mat[2][0] + q.mat[3][0] * q.mat[3][0]);
    rotation.mat.push_back(rotation_row);
    rotation_row.clear();
    return rotation;
}

matrix dQ(matrix q, matrix w)
{
    matrix dQuarter;
    dQuarter.column = 1, dQuarter.row = 4;
    vector<double> derivative;
    derivative.push_back(0.5 * (-q.mat[1][0] * w.mat[0][0] - q.mat[2][0] * w.mat[1][0]) - q.mat[3][0] * w.mat[2][0]);
    dQuarter.mat.push_back(derivative);
    derivative.clear();
    derivative.push_back(0.5 * (q.mat[0][0] * w.mat[0][0] - q.mat[3][0] * w.mat[1][0]) + q.mat[2][0] * w.mat[2][0]);
    dQuarter.mat.push_back(derivative);
    derivative.clear();
    derivative.push_back(0.5 * (q.mat[3][0] * w.mat[0][0] + q.mat[0][0] * w.mat[1][0]) - q.mat[1][0] * w.mat[2][0]);
    dQuarter.mat.push_back(derivative);
    derivative.clear();
    derivative.push_back(0.5 * (-q.mat[2][0] * w.mat[0][0] + q.mat[1][0] * w.mat[1][0]) + q.mat[0][0] * w.mat[2][0]);
    dQuarter.mat.push_back(derivative);
    derivative.clear();
    return dQuarter;
}