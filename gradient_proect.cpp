#include <iostream>
#include <vector>
#include <cmath>
#include <tuple>

#include "matrix.h"

#define h 0.00000001

double count_alpha(const auto& func, const Matrix& parameters, const Matrix& d_parameters, int& alpha_iterations) {
    double alpha = 0.5;
    for (;func(parameters + (-1)*alpha*d_parameters) > func(parameters);alpha_iterations++)
        alpha /= 2;
    return alpha;
}

bool epsilon(Matrix lhs, const Matrix& rhs) {
    for (auto i = 0; i < lhs.GetNumColumns(); i++)
        lhs[0][i] = lhs[0][i] - rhs[0][i];
    lhs = lhs*Transp(lhs);
    return lhs[0][0] > h;
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
        for(int i = 0; i < m.GetNumColumns(); i++)
            m[0][i] = m[0][i] / _parameters[i] + _centre[0][i];
        return m;
    }

    /* dont needed
    auto AntiProect(Matrix m) {
        for (int i =0; i < m.GetNumColumns(); i++)
        m[0][i] *= _parameters[i];
        return m;
    }*/

    double getRadius() const { return _radius; } 
    
    bool CheckBall(const Matrix& vector) {
        double sum = 0;
        for (int i = 0; i < vector.GetNumColumns(); i++) {
                sum += vector[0][i] * vector[0][i];
            }
        std::cout << sum << '\n';
        return fabs(sum - _radius) < h;
    }

    auto ProectFunc(const auto& func) {
        return func(_parameters, _centre);
    }
    auto ProectBall () {
        double& radius = _radius;
        auto& parameters = _parameters;
        return [&radius, &parameters](Matrix vector) {
            double sum = 0;
            for (int i = 0; i < vector.GetNumColumns(); i++) {
                sum += vector[0][i] * vector[0][i] / parameters[i];
            }
            if (sum < radius)
                return vector;
            return radius* vector /Norm(vector);
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
        int alpha_iterations = 0;
        auto alpha = count_alpha(func, parameters, d_parameters, alpha_iterations);
        parameters = parameters + (-1)*alpha * d_parameters;
        parameters = ProectBall(parameters);
        std::cout << "Parameters on iteration " << ++i << ": " << proection.ProectMatrix(parameters);
        std::cout << "Iterations for alpha: " << alpha_iterations << '\n';
        std::cout << "Function: " << func(parameters) << "\n\n";
    } while( epsilon(parameters, previous_parameters) && fabs(func(parameters) - func(previous_parameters)) > h);
    /*if (proection.CheckBall(parameters))
        std::cout << proection.ProectMatrix(parameters) << "is on the balls edge\n";*/
}


int main() {

    auto Func = [](const std::vector<double>& params, const Matrix& centre) {
        return [params, centre](const Matrix& m) -> double {
            auto x = m[0][0]/params[0] + centre[0][0];
            auto y = m[0][1]/params[1] + centre[0][1];
            auto z = m[0][2]/params[2] + centre[0][2];
            return x +4*y + z;
        };
    };

    
    auto dFunc = [](Matrix m, const auto& Function) -> Matrix {
        //double h = 0.0000001;
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
    //x^2 + 3y^2 +2z^2 <= 1
    //Matrix({{begin from?}}), Our Func, derivate of func
    for (int i = 1; i <= 100'000; i=i*10) {
        std::cout << "### For (" << i << ", " << i << ", " << i << ") ###\n"; 
        GradientMethod(Matrix({{static_cast<double>(i), static_cast<double>(i), static_cast<double>(i)}}), Func, dFunc,
                    //Proect({{Elispe -> ball}}, Matrix{{coordinates of centre}}, radius of ball)
                    Proection({sqrt(1), sqrt(3), sqrt(2)}, Matrix({{0, 0, 0}}), 1));
    }
    return 0;
}