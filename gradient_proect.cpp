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
    std::vector<double> _parameters;
    Matrix _centre;
    double _radius;
public:
    Proection() = default;
    Proection(std::vector<double>&& parameters,
        Matrix&& centre, double r) :
        _parameters(parameters), _centre(centre), _radius(r)
        {}
    auto ProectMatrix(Matrix m) {
        for (int i = 0; i < m.GetNumColumns(); i++)
        m[0][i] /= _parameters[i];
        return m;
    }

    auto AntiProect(Matrix m) {
        for (int i =0; i < m.GetNumColumns(); i++)
        m[0][i] *= _parameters[i];
        return m;
    }
    auto ProectFunc(const auto& func) {
        return func(_parameters);
    }
    auto ProectBall () {
        Matrix& centre = _centre;
        double& radius = _radius;
        return [&centre, &radius](Matrix vector) {
            return centre + radius*(vector + (-1)*centre)/(Norm(vector + (-1)*centre));
            };
    };
};


void GradientMethod( Matrix&& parameters, const auto& function, const auto& dfunc, auto&& proection)
{
    auto func = proection.ProectFunc(function);
    auto ProectBall = proection.ProectBall();
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

    auto Func = [](const std::vector<double>& params) {
        return [params](const Matrix& m) -> double {
            auto x = m[0][0]/params[0];
            auto y = m[0][1]/params[1];
            auto z = m[0][2]/params[2];
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

    GradientMethod(Matrix({{-10000, 10000, 10000}}), Func, dFunc, Proection({1, sqrt(3), sqrt(2)}, Matrix({{0, 0, 0}}), 1));
    return 0;
}