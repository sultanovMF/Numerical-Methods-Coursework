#include <iostream>
#include <span>
#include <fstream>
#include <string>
#include <iterator>
#include <vector>
#include <functional>
#include <indicators/indicators.hpp>
#include <filesystem>
#include <chrono>
#include <ctime>  
#include <omp.h>

class FDHandler {
public:
    // handler must return the accurance
    virtual double handle(std::vector<double> u, int n) = 0;
};

class FDHandlerW2F  : public FDHandler{
private:
    int Nt;
    int Nx;
    double dx;
    double dt;
    int i = 0;
public:
    FDHandlerW2F(int Nx, double L, int Nt, double T) 
        : Nx(Nx), Nt(Nt), dt(T / (Nt - 1)), dx(L / (Nx - 1)){
    }
    double handle(std::vector<double> u, int n) {
        std::cout << "Layer: " << i + 1 << " / " << Nt << std::endl;
        static std::filesystem::path path{
            u8"C:\\Users\\mur-m\\OneDrive\\Documents\\Courseworks\\Coursework V\\Numerical Methods Coursework\\numerical_result.txt"
        };
        static std::ofstream output_file(path);
        static std::ostream_iterator<std::string> output_iterator(output_file, "\n");
        int j;
        for (j = 0; j < Nx; ++j) {
            output_file << j * dx << "\t" << n * dt << "\t" << u[j] << "\n";
        }

        ++i;
        return 0;
    }
};

class FDWaveSolver {
public:
    FDWaveSolver(double a, std::function<double(double x, double t)> f, FDHandler* handler) : a(a), handler(handler), f(f) {
    }
protected:
    virtual void set_ic() = 0;
    virtual void set_bc() = 0;

    std::vector<double> u;
    std::vector<double> u_prev;
    std::vector<double> u_prev_prev;

    std::function<double(double x, double t)> f; // heterogenity
    double a;  // wave eq coeff

    int Nt; // number of time layers
    int Nx; // number of time layers

    double dt; // step by time
    double dx; // step by x coordinate

    FDHandler* handler;
public:
    void solve() {
        double C = a * dt / dx;

        if (C > 1.) throw std::invalid_argument("Courant number C > 1, so scheme is unstable!");

        double C2 = C * C;
        for (int n = 1; n < Nt - 1; ++n) {
            double result = 0;
//#pragma omp parallel for reduction(+:result) schedule(auto)
            for (int i = 1; i < Nx - 1; ++i) {
                u[i] = 2. * u_prev[i] - u_prev_prev[i]
                    + C2 * (u_prev[i + 1] - 2. * u_prev[i] + u_prev[i - 1])
                    + dt * dt * f(i * dx, n * dt);
            }
            set_bc();

            // обработка слоя пользователем
            this->handler->handle(u, n + 1);

            u_prev_prev = u_prev;
            u_prev = u;
        }

    }
};

class Wave1d1bc : public FDWaveSolver{
public:
    Wave1d1bc(double a, int Nx, double L, int Nt, double T,
        std::function<double(double)> psi_1,
        std::function<double(double)> psi_2,
        std::function<double(double)> phi_0,
        std::function<double(double)> phi_L,
        std::function<double(double x, double t)> f,
        FDHandler* handler = 0) : FDWaveSolver(a, f, handler) {
        this->Nt = Nt;
        this->Nx = Nx;

        this->dx = L / (Nx - 1);
        this->dt = T / (Nt - 1);

        this->psi_1 = psi_1;
        this->psi_2 = psi_2;
        this->phi_0 = phi_0;
        this->phi_L = phi_L;

        u_prev_prev.resize(Nx);
        u_prev.resize(Nx);
        u.resize(Nx);

        // set initial condition and special firmula for first and second layers
        set_ic();
        // enforce boundary condition for second layer
        set_bc();

        this->handler->handle(u_prev, 0);
        this->handler->handle(u, 1);

        u_prev_prev = u_prev;
        u_prev = u;
    }
private:
    std::function<double(double)> psi_1;
    std::function<double(double)> psi_2;
    std::function<double(double)> phi_0;
    std::function<double(double)> phi_L;

    void set_ic() {
        for (int i = 0; i < Nx; ++i) {
            u_prev[i] = psi_1(dx * i);
        }

        // Special formula for derrivative condition
        for (int i = 1; i < Nx - 1; ++i) {
            u[i] = psi_1(i * dx) + dt * psi_2(i * dx);
        }
    };
    // n is number of time layer
    void set_bc() {
        u[0] = phi_0(0);
        u[Nx - 1] = phi_L(Nt * dt);
    }
};

double zero(double) { return 0; }
void write_vec_to_file(const std::vector<double>& u, std::string spath) {
    // auto start = std::chrono::system_clock::now();

    //static std::filesystem::path path{ 
    //    u8"C:\\Users\\mur-m\\OneDrive\\Documents\\Courseworks\\Coursework V\\Numerical Methods Coursework\\result.txt"
    //};
    static std::filesystem::path path{ spath };
    static std::ofstream output_file(path);
    static std::ostream_iterator<std::string> output_iterator(output_file, "\n");
    //std::copy(u.begin(), u.end(), output_iterator);
    for (const auto& e : u) {
        output_file << e << "\t";
    }
    output_file << "\n";
}

double foo(double x) {
    return 4 - x;
}

double zero2(double, double) {
    return 0;
}
int main() {
    //indicators::show_console_cursor(false);
    double a = 1;
    int Nx = 100;
    int Nt = 1000;
    double L = 4;
    double T = 1;
   
    //FDHandlerW2F handler(Nx, L, Nt, T);
    //Wave1d1bc solver(a, Nx, L, Nt, T, zero, foo, zero, zero, zero2, &handler);
    //solver.solve();
    //Wave1d1bc solver(1000, 4., 500, 2., zero, foo, zero, zero);
   const double pi = 3.14159265358979323846;

    std::filesystem::path path{ u8"C:\\Users\\mur-m\\OneDrive\\Documents\\Courseworks\\Coursework V\\Numerical Methods Coursework\\target_result.txt" };
    std::ofstream output_file(path);
    std::ostream_iterator<std::string> output_iterator(output_file, "\n");

    double dx = L / (Nx - 1);
    double dt = T / (Nt - 1);
    for (int i = 0; i < Nt; ++i) {
        double t = i * dt;
        for (int j = 0; j < Nx; ++j) {
            double x = j * dx;
            double u = 0.;
            for (int k = 1; k < 100000; ++k) {
                u += 32 * std::sin(k * pi * x / 4.) * std::sin(k * pi * t / 4) / (pi * pi * k * k);
            }
            output_file << x << "\t" << t << "\t" << u << "\n";
        }
        std::cout << "Layer " << i + 1 << " / " << Nt << std::endl;
    }
;   return 0;
}