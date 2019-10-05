
#include "angleSolver.hpp"

int main( int argc, const char** argv )
{      

    angleSolver angle_solver;
    double alpha = angle_solver.get_ALPHA();
    double thet = angle_solver.get_THET();
    cout << "alpha is: "<<alpha*180/M_PI<<endl;
    cout << "thet is: "<<thet*180/M_PI<<endl;
    
}