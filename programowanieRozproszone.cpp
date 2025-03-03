#include <iostream>
#include <vector>
#include <complex>
#include <random>
#include <cmath>
#include <windows.h>
#include <thread>
#include "MetodaJacobiego.h"

int main() {
    LARGE_INTEGER cpuFreqHz, t1, t2;
    long double policzonyCzas;

    int liczbaProb = 1;

    int maxThreads = std::thread::hardware_concurrency();
    if (maxThreads == 0) maxThreads = 8;

    QueryPerformanceFrequency(&cpuFreqHz);
    MetodaJacobiego system = MetodaJacobiego();
	system.setRozmiar(100);

    std::cout << " =========== " << system.getRozmiar() << " =========== " << std::endl;
    for (int proba = 1; proba <= liczbaProb; proba++) {
        std::cout << " ========== PROBA  " << proba << " ========== " << std::endl;

        QueryPerformanceCounter(&t1);
        system.generujMacierz();
        QueryPerformanceCounter(&t2);

        policzonyCzas = (static_cast<long double>(t2.QuadPart - t1.QuadPart)) / cpuFreqHz.QuadPart * 1000;
        std::cout << "macierze wygenerowane: " << policzonyCzas << " ms" << std::endl;

        QueryPerformanceCounter(&t1);
        std::vector<std::complex<double>> wynikSekwencyjny = system.oblicz();
        QueryPerformanceCounter(&t2);

        policzonyCzas = (static_cast<long double>(t2.QuadPart - t1.QuadPart)) / cpuFreqHz.QuadPart * 1000;
        std::cout << "Czas wykonania sekwencyjnie: " << policzonyCzas << " ms" << std::endl;
        std::cout << " ... " << std::endl;


        QueryPerformanceCounter(&t1);
        std::vector<std::complex<double>> wynikOpenMP = system.obliczOpenMP();
        QueryPerformanceCounter(&t2);

        policzonyCzas = (static_cast<long double>(t2.QuadPart - t1.QuadPart)) / cpuFreqHz.QuadPart * 1000;
        std::cout << "Czas wykonania OpenMP: " << policzonyCzas << " ms. ";

        if (system.porownajWektory(wynikSekwencyjny, wynikOpenMP)) {
            std::cout << "Wektory sa takie same" << std::endl;
        }
        else {
            std::cout << "Wektory sa rozne" << std::endl;
        }
        std::cout << " ... " << std::endl;

        std::vector<std::complex<double>> wynikThreads;
        for (int i = 6; i <= maxThreads; i += 4) {
            system.setWatki(i);
            wynikThreads.clear();

            QueryPerformanceCounter(&t1);
            wynikThreads = system.obliczWatki();
            QueryPerformanceCounter(&t2);

            policzonyCzas = (static_cast<long double>(t2.QuadPart - t1.QuadPart)) / cpuFreqHz.QuadPart * 1000;
            std::cout << "Czas wykonania " << system.getWatki() << " threads: " << policzonyCzas << " ms. ";

            if (system.porownajWektory(wynikSekwencyjny, wynikThreads)) {
                std::cout << "Wektory sa takie same" << std::endl;
            }
            else {
                std::cout << "Wektory sa rożne" << std::endl;
            }
            std::cout << " ... " << std::endl;
        }

        std::cout << " ========== KONIEC " << proba << " ========== " << std::endl;
    }

    return 0;
}

