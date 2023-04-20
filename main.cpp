#include <string>
#include <iostream>
#include <fstream>
#include <cstring>
#include "ising.h"
#include "ising.cpp"

#define T_STEPS 100
#define T_MIN 1
#define T_MAX 3.0

typedef struct IsingResult {
    MCResult result;
    double temperature;
} IsingResult;

template<int N>
void generate_result(IsingResult* results) {
    #pragma omp parallel for
    for (int t_step = 0; t_step < T_STEPS; t_step++) {
        double t = (T_MAX - T_MIN) / (T_STEPS - 1) * t_step + T_MIN;
        double beta = 1.0 / t;
        auto i = new Ising<N>(time(nullptr));
        i->thermalize(beta);
        MCResult res = i->mc_sweep(beta);
        results[t_step] = {res, t};
        delete i;
    }
}

void stream_results_csv(std::basic_ostream<char>& stream, IsingResult* results) {
    stream << "temperature,energy,magnetisation,susceptibility,specific_heat" << std::endl;
    for (int row = 0; row < T_STEPS; row++) {
        stream << results[row].temperature
                  << ", " << results[row].result.energy
                  << ", " << results[row].result.magnetisation
                  << ", " << results[row].result.susceptibility
                  << ", " << results[row].result.specific_heat
                  << std::endl;
    }
}

void print_results_csv(IsingResult* results) {
    stream_results_csv(std::cout, results);
}

void save_to_file_csv(std::string& file_name, IsingResult* results) {
    std::ofstream out_file(file_name);
    stream_results_csv(out_file, results);
}

typedef void (*ResCaller)(IsingResult*);

typedef struct MenuItem {
    std::string menuText;
    ResCaller caller;
} MenuItem;

#define IMENU_NUM_ITEMS 5
#define GEN_RES_ROW(i) {#i, [](IsingResult* res) {generate_result<i>(res);}}
static const MenuItem i_menu_items[] = {
        GEN_RES_ROW(16),
        GEN_RES_ROW(32),
        GEN_RES_ROW(64),
        GEN_RES_ROW(128),
        GEN_RES_ROW(256),
};

#define PMENU_NUM_ITEMS 2
static const MenuItem p_menu_items[] = {
        {"Print to stdout", [](IsingResult* res) { print_results_csv(res);}},
        {"Save", [](IsingResult* res) {
            std::string filename;
            std::cout << "Please enter output file name: ";
            std::cin >> filename;
            filename = "out/" + filename + ".csv";
            std::cout << std::endl;
            save_to_file_csv(filename, res);
        }}
};

unsigned char show_menu(const MenuItem* menu, const std::string& menu_title, unsigned int menu_length, bool printSel = true) {
    unsigned int selection;
    std::cout << menu_title << std::endl;
    for (unsigned int ln = 0; ln < menu_length; ln++) {
        std::cout << ln << " > " << menu[ln].menuText << std::endl;
    }
    std::cin >> selection;
    std::cout << std::endl;
    if (selection < 0 || selection >= menu_length) {
        std::cout << "Wrong selection!" << std::endl;
        return show_menu(menu, menu_title, menu_length, true);
    }
    return selection;
}

int main(int argc, char* argv[]) {
    bool test_run = false;
    unsigned char gs = 3;
    for (int i = 0; i < argc; i++) {
        if (strcmp(argv[i], "test") == 0) {
            test_run = true;
        }
    }
    IsingResult isingResult[T_STEPS];
    if (!test_run)
        gs = show_menu(i_menu_items, "Select grid size.", IMENU_NUM_ITEMS);
    (i_menu_items[gs].caller)(isingResult);
    if (!test_run) {
        unsigned char ps = show_menu(p_menu_items, "What should be done with results", PMENU_NUM_ITEMS);
        (p_menu_items[ps].caller)(isingResult);
    }
    return 0;
}
