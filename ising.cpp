#include <iostream>
#include <cmath>
#include "ising.h"

template <int N>
constexpr double norm1 = 1.0 / (MC_SWEEP_COUNT * N * N);

/**
 * Initialize Ising model system
 * @tparam N Grid size (side)
 * @param seed Random seed for the Metropolis-Hastings algorithm
 */
template<int N>
Ising<N>::Ising(unsigned int seed)  : gen(seed), dist(0.0, 1.0) {
    // Set all spins up for initial configuration
    memset(grid, 1, N * N * sizeof(char));
}

/**
 * Draw current state of the system
 * @tparam N Grid size (side)
 * @param out Output stream
 */
template<int N>
void Ising<N>::draw(std::basic_ostream<char>& out) {
    for (int y = 0; y < N; y++) {
        for (int x = 0; x < N; x++) {
            out << (grid[x + N * y] == -1 ? "O " : "X ");
        }
        out << std::endl;
    }
}

/**
 * Precompute exponentials for the Metropolis-Hastings algorithm
 * @tparam N Grid size (side)
 * @param beta Inverse temperature
 * @param exps Resulting array of length 4
 */
template<int N>
void Ising<N>::precompute_exp(double beta, double *exps) {
    // This could even be reduced to two possibilities (4 and 8 energy delta)...
    for (int i = 1; i < 5; i++) {
        exps[i-1] = exp(-beta * 2 * i);
    }
}

/**
 * Computes internal energy of the spin system
 * @tparam N Grid size (side)
 * @return Energy of the spin system
 */
template<int N>
int Ising<N>::compute_energy() {
    int energy = 0;
    for (int y = 0; y < N; y++) {
        for (int x = 0; x < N; x++) {
            // Neighbour indices
            int ty = ((y - 1) % N + N) % N;
            int by = (y + 1) % N;
            int lx = ((x - 1) % N + N) % N;
            int rx = (x + 1) % N;
            energy -= grid[y * N + x] * (grid[y * N + lx] + grid[y * N + rx] + grid[ty * N + x] + grid[by * N + x]);
        }
    }
    return energy;
}
/**
 * Compute magnetization of the grid
 * @tparam N Grid size (side)
 * @return Magnetization of the spin grid
 */
template<int N>
int Ising<N>::compute_magnetization() {
    // Number of up - number of down
    int mg = 0;
    for (int i = 0; i < N*N; i++) {
        mg += grid[i];
    }
    return mg;
}

/**
 * Determines if Metropolis-Hastings accepts the spin change
 * @tparam N Grid size (side)
 * @param energy_delta_index Index of the exponential factor (energy delta / 2) or negative number
 * @param exps Precomputed exponential factors
 * @return If spin should be changed
 */
template<int N>
bool Ising<N>::metropolisAccept(signed char energy_delta_index, const double* exps) {
    return (energy_delta_index <= 0.0) || (dist(gen) < exps[energy_delta_index - 1]);
}
/**
 * Perform one sweep through the grid.
 * This method tries to change every spin in the grid and accepts the change
 * based on Metropolis-Hastings algorithm
 * @tparam N Grid size (side)
 * @param exps Precomputed exponential factors for the Metropolis-Hastings
 */
template<int N>
void Ising<N>::sweep(const double* exps) {
    for (int y = 0; y < N; y++) {
        for (int x = 0; x < N; x++) {
            // Neighbouring indices
            int ty = ((y - 1) % N + N) % N;
            int by = (y + 1) % N;
            int lx = ((x - 1) % N + N) % N;
            int rx = (x + 1) % N;

            // Index of the exponential factor (energy delta / 2), can be negative
            signed char change_idx = grid[y * N + x] * (
                    grid[y * N + lx] + grid[y * N + rx] + grid[ty * N + x] + grid[by * N + x]);

            bool isAccepted = metropolisAccept(change_idx, exps);

            grid[y * N + x] = - isAccepted * grid[y * N + x] + !isAccepted * grid[y * N + x];
        }
    }
}

/**
 * Thermalize system by doing THERMALIZE_SWEEP_COUNT sweeps
 * @tparam N Grid size (side)
 * @param beta Reverse temperature of the system
 */
template<int N>
void Ising<N>::thermalize(double beta) {
    double exps[4];
    precompute_exp(beta, exps); // Precompute exponential factors for this temperature
    for (int swp = 0; swp < THERMALIZE_SWEEP_COUNT; swp++) {
        sweep(exps);
    }
}

/**
 * Perform MC_SWEEP_COUNT sweeps of the grid and calculate MCResult values
 * @tparam N Grid size (side)
 * @param beta Reverse temperature of the system
 * @return Result containing energy, magnetisation, susceptibility and specific heat
 */
template<int N>
MCResult Ising<N>::mc_sweep(double beta) {
    double exps[4];
    // Precompute exponential factors for this temperature
    precompute_exp(beta, exps);
    long energy_acc = 0;
    long magnetisation_acc = 0;
    // This can get huge - normalize on the way!
    double energy_sqr_acc = 0;
    double magnetisation_sqr_acc = 0;

    for (int swp = 0; swp < MC_SWEEP_COUNT; swp++) {
        sweep(exps);
        // Rounding errors here are OK'ish
        long energy = compute_energy() / 4;
        energy_acc += energy;
        energy_sqr_acc += energy * energy * norm1<N>;
        long mag = compute_magnetization();
        magnetisation_acc += mag;
        magnetisation_sqr_acc += (mag * mag) * norm1<N>;
    }

    return {
            (double) energy_acc * norm1<N>,
            (double) magnetisation_acc * norm1<N>,
            beta * ((double) magnetisation_sqr_acc - (double) magnetisation_acc * (double) magnetisation_acc * norm1<N> / MC_SWEEP_COUNT),
            beta * beta * (energy_sqr_acc - (double) energy_acc * (double) energy_acc * norm1<N> / MC_SWEEP_COUNT),
    };
}