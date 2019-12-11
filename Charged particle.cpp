#include <iostream>
#include <fstream>

using namespace std;
#include <math.h>

double calcOmx(double x, double y, double z){
    double d = 7.94e24, q = 1.9e-19, m = 9.1e-31;
    double bx = -3*d*x*z/(pow((pow(x,2.0)+pow(y,2.0)+pow(z,2.0)),2.5));
    return bx;
}
double calcOmy(double x, double y, double z){
    double d = 7.94e24, q = 1.9e-19, m = 9.1e-31;
    double by = -3*d*y*z/(pow((pow(x,2.0)+pow(y,2.0)+pow(z,2.0)),2.5));
    return by;
}
double calcOmz(double x, double y, double z){
    double d = 7.94e24, q = 1.9e-19, m = 9.1e-31;
    double bz = d*(pow(x,2.0)+pow(y,2.0)-2*pow(z,2.0))/(pow((pow(x,2.0)+pow(y,2.0)+pow(z,2.0)),2.5));
    return bz;
}
double calcVT(double x, double y, double z, double vx, double vy, double vz){
    double t[3] = {-3*x*z/((pow((pow(x,2.0)+pow(y,2.0)+pow(z,2.0)),0.5))*(pow((pow(x,2.0)+pow(y,2.0)+4*pow(z,2.0)),0.5))),
                    -3*y*z/((pow((pow(x,2.0)+pow(y,2.0)+pow(z,2.0)),0.5))*(pow((pow(x,2.0)+pow(y,2.0)+4*pow(z,2.0)),0.5))),
                    (pow(x,2.0)+pow(y,2.0)-2*pow(z,2.0))/((pow((pow(x,2.0)+pow(y,2.0)+pow(z,2.0)),0.5))*(pow((pow(x,2.0)+pow(y,2.0)+4*pow(z,2.0)),0.5)))};
    double vt = vx*t[0]+vy*t[1]+vz*t[2];
    return vt;
}
double calcVF(double x, double y, double z, double vx, double vy, double vz){
    double f[3] = {-y/(pow((pow(x,2.0)+pow(y,2.0)),0.5)), x/(pow((pow(x,2.0)+pow(y,2.0)),0.5)), 0};
    double vf = vx*f[0]+vy*f[1]+vz*f[2];
    return vf;
}

int main(int argc, char* argv[]){
    double y[6];
    double r = 6378000;
    double dt = 1e-8, om= 6.6e5;
    y[0]= 2*r;
    y[1]= 100;
    y[2]= 100;
    double k1[6],k2[6];

    ofstream myfile;
  myfile.open ("vfff.txt");
 // for (int j = 0; j<6; j++){
 // myfile << y[j]<<"\t";
 // }
 // myfile << "\n";
int c = 0;

    cin >> y[3] >> y[4]>> y[5];
    for (int i=0; i<1e6; i++)
    {
        k1[0]=dt*y[3];
        k1[1]=dt*y[4];
        k1[2]=dt*y[5];
        k1[3]=dt*(calcOmz(y[0], y[1], y[2])*y[4]-calcOmy(y[0], y[1], y[2])*y[5]);
        k1[4]=dt*(calcOmx(y[0], y[1], y[2])*y[5]-calcOmz(y[0], y[1], y[2])*y[3]);
        k1[5]=dt*(calcOmy(y[0], y[1], y[2])*y[3]-calcOmx(y[0], y[1], y[2])*y[4]);

        k2[0]=dt*(y[3]+0.5*k1[3]);
        k2[1]=dt*(y[4]+0.5*k1[4]);
        k2[2]=dt*(y[5]+0.5*k1[5]);
        k2[3]=dt*(calcOmz(y[0]+0.5*k1[0], y[1]+0.5*k1[1], y[2]+0.5*k1[2])*(y[4]+0.5*k1[4])-calcOmy(y[0]+0.5*k1[0], y[1]+0.5*k1[1], y[2]+0.5*k1[2])*(y[5]+0.5*k1[5]));
        k2[4]=dt*(calcOmx(y[0]+0.5*k1[0], y[1]+0.5*k1[1], y[2]+0.5*k1[2])*(y[5]+0.5*k1[5])-calcOmz(y[0]+0.5*k1[0], y[1]+0.5*k1[1], y[2]+0.5*k1[2])*(y[3]+0.5*k1[3]));
        k2[5]=dt*(calcOmy(y[0]+0.5*k1[0], y[1]+0.5*k1[1], y[2]+0.5*k1[2])*(y[3]+0.5*k1[3])-calcOmx(y[0]+0.5*k1[0], y[1]+0.5*k1[1], y[2]+0.5*k1[2])*(y[4]+0.5*k1[4]));

        y[0]+=k2[0];
        y[1]+=k2[1];
        y[2]+=k2[2];
        y[3]+=k2[3];
        y[4]+=k2[4];
        y[5]+=k2[5];

c++;
//  for (int j = 0; j<6; j++){
//  myfile << y[j]<<"\t"<<c;
//  }
//  myfile << "\n";
myfile << c*dt <<"\t"<<calcVF(y[0], y[1], y[2], y[3], y[4], y[5])<<"\n";
}

myfile.close();


return 0;
};
