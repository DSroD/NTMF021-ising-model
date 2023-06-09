#ifndef ISING_H
#define ISING_H

#include <bitset>
#include <random>

// Number of steps for thermalization
#define THERMALIZE_SWEEP_COUNT 100
// Number of steps for MC averaging
#define MC_SWEEP_COUNT 3000

/**
 * Result of the MC simulation of one system
 */
typedef struct MCResult {
    double energy;
    double magnetisation;
    double susceptibility;
    double specific_heat;
} MCResult;

template<int N>
class Ising {
using Grid = signed char[N * N];
public:
    explicit Ising(unsigned int seed);
    int compute_energy();
    int compute_magnetization();
    void thermalize(double beta);
    MCResult mc_sweep(double beta);
    void draw(std::basic_ostream<char>& out);

private:
    Grid grid{};
    std::minstd_rand gen;
    std::uniform_real_distribution<double> dist;
    void sweep(const double* exps);
    bool metropolisAccept(signed char energy_delta_index, const double* exps);
    void precompute_exp(double beta, double* exps);
};

#endif //ISING_H
