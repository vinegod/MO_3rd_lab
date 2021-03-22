#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <cmath>
#include <iomanip>

using std::vector;

class Matrix {
private:
  int num_rows;
  int num_columns;

  vector<vector<double>> elements;

public:
  Matrix();
  Matrix(vector<vector<double>>&& elements_);
  Matrix(int num_rows, int num_columns);
  void Reset(int num_rows_, int num_columns_);
  int GetNumRows() const;
  int GetNumColumns() const;
  const std::vector<double>& operator[](size_t i) const;
  vector<double>& operator[](size_t i);
};

bool operator==(const Matrix& one, const Matrix& two);
Matrix operator+(const Matrix& one, const Matrix& two);
Matrix operator+(const Matrix& one, double two);
Matrix operator*(const Matrix& one, const Matrix& two);
Matrix Transp(const Matrix& matrix);
Matrix MatrixE(const Matrix&);
double Norm(const Matrix& lhs);
Matrix operator*(double one, const Matrix& two);

Matrix row_matrix(const vector<double>& v);
Matrix column_matrix(const vector<double>& v);

Matrix operator/(const Matrix& matrix, double x);

std::istream& operator>>(std::istream& in, Matrix& matrix);
std::ostream& operator<<(std::ostream& out, const Matrix& matrix);


#endif