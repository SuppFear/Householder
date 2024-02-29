#include "householder.hpp"

int main(){
    Dense_Matrix A(3, {1, 2, 3, 4, 5, 6, 7, 8, 9});
    A.print();
    std::vector<Dense_Matrix> qr = QR(A);
    qr[0].print();
    qr[1].print();
}