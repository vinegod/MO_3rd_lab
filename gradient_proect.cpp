#include <iostream>
#include <vector>
#include <cmath>
#include <tuple>

#include "matrix.h"


double count_alpha(const auto& func, const Matrix& parameters, const Matrix& d_parameters) {
    double alpha = 0.5;
    int i = 0;
    for (;func(parameters + (-1)*alpha*d_parameters) > func(parameters);i++)
        alpha /= 2;
    std::cout << "Iterations for alpha: " << i << '\n';
    return alpha;
}

bool epsilon(Matrix lhs, const Matrix& rhs) {
    for (auto i = 0; i < lhs.GetNumColumns(); i++)
        lhs[0][i] = lhs[0][i] - rhs[0][i];
    lhs = lhs*Transp(lhs);
    return lhs[0][0] > 0.0000001;
}


class Proection {
private:
    double _x1, _y1, _z1;
public:
    Proection() : _x1(1), _y1(1), _z1(1) {}
    Proection(double x1, double y1, double z1) : _x1(x1), _y1(y1), _z1(z1) {}
    auto ProectMatrix(Matrix m) {
        m[0][0] /= _x1;
        m[0][1] /= _y1;
        m[0][2] /= _z1;
        return m;
    }

    auto AntiProect(Matrix m) {
        m[0][0] *= _x1;
        m[0][1] *= _y1;
        m[0][2] *= _z1;
        return m;
    }
    auto ProectFunc(const auto& func) {
        return func(_x1, _y1, _z1);
    }
    auto ProectBall (Matrix&& centre, double radius) {
        return [centre, radius](Matrix vector) {
            return centre + radius*(vector + (-1)*centre)/(Norm(vector + (-1)*centre));
            };
    };
};


void GradientMethod( Matrix&& parameters, const auto& function, const auto& dfunc, auto&& proection)
{
    auto func = proection.ProectFunc(function);
    auto ProectBall = proection.ProectBall(Matrix({{0, 0, 0}}), 1);
    int i = 0;
    Matrix previous_parameters;
    do
    {
        previous_parameters = parameters;
        auto d_parameters = dfunc(parameters, func);
        auto alpha = count_alpha(func, parameters, d_parameters);
        parameters = parameters + (-1)*alpha * d_parameters;
        parameters = ProectBall(parameters);
        std::cout << "Parameters on iteration " << ++i << ": " << proection.AntiProect(parameters) << '\n';
    } while( epsilon(parameters, previous_parameters) /*&& (func(parameters) - func(previous_parameters)) > 0.000001*/);

}


int main() {

    auto Func = [](double x1=1, double y1=1, double z1=1) {
        return [x1, y1, z1](const Matrix& m) -> double {
            auto x = m[0][0]/x1;
            auto y = m[0][1]/y1;
            auto z = m[0][2]/z1;
            return x + 4*y + z;
        };
    };

    
    auto dFunc = [](Matrix m, const auto& Function) -> Matrix {
        double h = 0.0000001;
        Matrix result(m.GetNumRows(), m.GetNumColumns());
        for (auto i = 0; i < m.GetNumColumns(); i++) {
            m[0][i] += h;
            result[0][i] = Function(m);
            m[0][i] -= 2*h;
            result[0][i] -= Function(m);
            m[0][i] += h;
        }
        return result / (2 * h);
    };

    GradientMethod(Matrix({{-110, 1000, 100}}), Func, dFunc, Proection(1, sqrt(3), sqrt(2)));
    return 0;
}