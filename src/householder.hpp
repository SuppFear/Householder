#include "../Dense-CSR/dense_matrix.h"
#include "../Dense-CSR/vector_operations.h"
#include <cmath>

std::vector<Dense_Matrix> QR (const Dense_Matrix& A){
    std::vector<Dense_Matrix> QR = {A, A};
    size_t n  = A.columns();
    std::vector<double> a = A.get_matrixElemenets();
    size_t r = a.size() / n;
    std::vector<std::vector<double>> vecs = std::vector<std::vector<double>>(n-1);
    for(size_t i = 0; i < n - 1; i++){
        std::vector<double> x = std::vector<double>(n-i);
        for(size_t j = 0; j < r - i; j++){
            x[j] = a[(i+j)*n+i];
        }
        std::vector<double> v = x;
        if (x[0] >= 0){
            v[0] += std::sqrt(x*x);
        }
        else{
            v[0] -= std::sqrt(x*x);
        }
        vecs[i] = v;
        double vx = v * x;
        double vv = v * v;
        for(size_t k = 0; k < r - i; k++){
            a[(k+i)*n+i] = x[k] - (2.0 * v[k] * vx) / vv;
        }
        for(size_t m = i + 1; m < n; m++){
            for(size_t j = 0; j < r - i; j++){
                x[j] = a[(j+i)*n+m];
            }
            vx = v * x;
            for(size_t k = 0; k < r - i; k++){
                a[(k+i)*n+m] = x[k] - (2.0 * v[k] * vx) / vv;
            }
        }
    }

    std::vector<double> q = std::vector<double>(r*n);
    for(size_t i = 0; i < n; i++){
        for(size_t j = 0; j < r; j++){
            if (i == j){
                q[i*n+j] = 1.0;
            }
            else{
                q[i*n+j] = 0.0;
            }
        }
    }
    std::vector<double> y = std::vector<double>(n*r);
    for(size_t i = 0; i < n; i++){
        for (size_t j = 0; j < r; j++){
            y[i*n+j] = vecs[0][i] * vecs[0][j];
        }
    }
    q = q - (2.0 / (vecs[0] * vecs[0])) * y;
    for(size_t i = 1; i < n - 1; i++){
        for(size_t j = 0; j < r; j++){
            std::vector<double> y = std::vector<double>(n-i);
            for (size_t k = 0; k < r - i; k++){
                y[k] = q[j*n+i+k];
            }
            double vv = vecs[i] * vecs[i];
            double vx = vecs[i] * y;
            for (size_t k = 0; k < r - i; k++){
                q[j*n+i+k] = y[k] - (2.0 * vx * vecs[i][k]) / vv;
            }
        }
    }

    QR[0] = Dense_Matrix(n, q);
    QR[1] = Dense_Matrix(n, a);
    return QR;
}