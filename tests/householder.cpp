#include "../src/householder.hpp"

int main(){
    Dense_Matrix A(3, {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0});
    A.print();
    std::vector<Dense_Matrix> qr = QR(A);
    qr[0].print();
    qr[1].print();
    Dense_Matrix a(3, {2.0, -4.0, 3.0, 1.0, -2.0, 4.0, 3.0, -1.0, 5.0});
    a.print();
    std::vector<double> b = {1.0, 3.0, 2.0};
    std::vector<double> x = Solve_with_QR(a, b);
    for (size_t i = 0; i < x.size(); i++){
        std::cout<<x[i]<<" ";
    }
}