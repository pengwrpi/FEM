#include <iostream>
#include <cstdlib>

void allocate_mem(double (&x)[2])
{
    x[0] = 1.0;
    x[1] = 2.0;
}

int main(int argc, char* argv[])
{
    double n[2];
    n[0] = 0.0;
    n[1] = 3.0;
    allocate_mem(n);
    std::cout << n[0] << " " << n[1] << std::endl;
#if 0
    double n[][2];
    n[0][0] = 10.0;
    n[0][1] = 20.0;
    n[1][0] = 30.0;
    n[1][1] = 40.0;
    std::cout << n[0][0] << " " << n[0][1] << std::endl; 
    std::cout << n[1][0] << " " << n[1][1] << std::endl; 
    allocate_mem(&n[0]);
    std::cout << n[0][0] << " " << n[0][1] << std::endl; 
    std::cout << n[1][0] << " " << n[1][1] << std::endl; 
#endif
    return 0;
}
