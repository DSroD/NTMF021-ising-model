#include "ising.h"
#include <cmath>
#include <iostream>

template <int N>
const double norm1 = 1.0 / (MC_SWEEP_COUNT * N * N);

template<int N>
int Ising<N>::compute_energy() {
    int energy = 0;
    for (int x = 0; x < N; x++) {
        for (int y = 0; y < N; y++) {
            signed char s_e = grid.test(x + y * N) - !grid.test(x + y * N);
            signed char de = 0;

            // Right neighbour
            int rx = (x + 1) % N;
            de += grid.test(rx + y * N) - !grid.test(rx + y * N);
            // Left neighbour
            int lx = ((x - 1) % N + N) % N;
            de += grid.test(lx + y * N) - !grid.test(lx + y * N);
            // Top neighbour
            int ty = ((y - 1) % N + N) % N;
            de += grid.test(x + ty * N) - !grid.test(x + ty * N);
            // Bottom neighbour
            int by = (y + 1) % N;
            de += grid.test(x + by * N) - !grid.test(x + by * N);

            energy -= s_e * de;
        }
    }
    return energy;
}

template<int N>
signed char Ising<N>::compute_energy_change_index(int x, int y) {
    signed char de = 0;
    int ty = ((y - 1) % N + N) % N;
    int by = (y + 1) % N;
    int lx = ((x - 1) % N + N) % N;
    int rx = (x + 1) % N;

    signed char sp_c = grid.test(x + y * N) - !grid.test(x + y * N);
    de += grid.test(rx + y * N) - !grid.test(rx + y * N);
    de += grid.test(lx + y * N) - !grid.test(lx + y * N);
    de += grid.test(x + ty * N) - !grid.test(x + ty * N);
    de += grid.test(x + by * N) - !grid.test(x + by * N);

    // We are computing index of the precomputed exponential, therefore this is not
    // multiplied by two
    return sp_c * de;
}

template<int N>
int Ising<N>::compute_magnetization() {
    // Number of up - number of down
    int mg = 2 * grid.count() - (N*N);
    return mg;
}

template<int N>
void Ising<N>::sweep(const double* exps) {
    for (int pos = 0; pos < N * N; pos++) {
        // try change spin for this position
        signed char delta_index = compute_energy_change_index(pos % N, pos / N);
        bool isAccepted = metropolisAccept(delta_index, exps);
        grid.set(pos, !grid.test(pos) * isAccepted + grid.test(pos) * !isAccepted);
    }
}

template<int N>
bool Ising<N>::metropolisAccept(signed char energy_delta_index, const double* exps) {
    return (energy_delta_index <= 0.0) || (dist(gen) < exps[energy_delta_index - 1]);
}

template<int N>
void Ising<N>::thermalize(double beta) {
    double exps[4];
    precompute_exp(beta, exps);
    for (int swp = 0; swp < THERMALIZE_SWEEP_COUNT; swp++) {
        sweep(exps);
    }
}

template<int N>
MCResult Ising<N>::mc_sweep(double beta) {
    double exps[4];
    precompute_exp(beta, exps);
    long energy_acc = 0;
    // This can get huge - normalize on the way!
    double energy_sqr_acc = 0;
    long magnetisation_acc = 0;
    long magnetisation_sqr_acc = 0;

    for (int swp = 0; swp < MC_SWEEP_COUNT; swp++) {
        sweep(exps);
        // Rounding errors here are OK'ish
        int energy = compute_energy() / 4;
        energy_acc += energy;
        energy_sqr_acc += energy * energy * norm1<N>;
        int mag = compute_magnetization();
        magnetisation_acc += mag;
        magnetisation_sqr_acc += (mag * mag);
    }

    return {
        (double) energy_acc * norm1<N>,
        (double) magnetisation_acc * norm1<N>,
        beta * ((double) magnetisation_sqr_acc - (double) magnetisation_acc * (double) magnetisation_acc / MC_SWEEP_COUNT) * norm1<N>,
        beta * beta * (energy_sqr_acc - (double) energy_acc * (double) energy_acc * norm1<N> / MC_SWEEP_COUNT),
    };
}

template<int N>
Ising<N>::Ising(unsigned int seed)  : gen(seed), dist(0.0, 1.0) {
    // Set all spins up for initial configuration
    grid.set();
}

template<int N>
void Ising<N>::draw(std::basic_ostream<char>& out) {
    for (int y = 0; y < N; y++) {
        for (int x = 0; x < N; x++) {
            out << (grid.test(x + y * N) ? "O " : "X ");
        }
        out << std::endl;
    }
}

template<int N>
void Ising<N>::precompute_exp(double beta, double *exps) {
    // This could even be reduced to two possibilities (4 and 8 energy delta)...
    for (int i = 1; i < 5; i++) {
        exps[i-1] = exp(-beta * 2 * i);
    }
}
