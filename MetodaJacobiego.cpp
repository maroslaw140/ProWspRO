#include "MetodaJacobiego.h"
#include <iostream>
#include <random>
#include <cmath>
#include <limits>
#include <thread>
#include <omp.h>

MetodaJacobiego::MetodaJacobiego() {
    //this->liczbaWatkow = std::thread::hardware_concurrency();
    this->minLimit = std::numeric_limits<double>::min();
    this->maxLimit = std::numeric_limits<double>::max();
}

int MetodaJacobiego::getWatki() {
	return this->liczbaWatkow;
}

void MetodaJacobiego::setWatki(int liczbaWatki) {
	this->liczbaWatkow = liczbaWatki;
}

int MetodaJacobiego::getRozmiar() {
	return this->rozmiar;
}

void MetodaJacobiego::setRozmiar(int rozmiar) {
	this->rozmiar = rozmiar;
}

void MetodaJacobiego::generujMacierz() {
    try {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(this->minLimit, this->maxLimit);

        this->macierzA.clear();
        this->macierzA.resize(this->rozmiar, std::vector<std::complex<double>>(this->rozmiar));

        this->wektorB.clear();
        this->wektorB.resize(this->rozmiar);


        #pragma omp parallel for
        for (int i = 0; i < this->rozmiar; ++i) {
            this->wektorB[i] = std::complex<double>(dis(gen), dis(gen));
        }

        #pragma omp parallel for
        for (int i = 0; i < this->rozmiar; ++i) {
            for (int j = 0; j < this->rozmiar; ++j) {
                this->macierzA[i][j] = std::complex<double>(dis(gen), dis(gen));
            }
        }

        // Zapewnienie, ¿e macierz A jest dobrze uwarunkowana
        #pragma omp parallel for
        for (int i = 0; i < this->rozmiar; ++i) {
            this->macierzA[i][i] = std::complex<double>(dis(gen) + 1.0, dis(gen) + 1.0);
        }
    }
    catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        return;
    }
}


std::vector<std::complex<double>> MetodaJacobiego::oblicz() {
    std::vector<std::complex<double>> X(this->rozmiar, { 0.0, 0.0 });
    std::vector<std::complex<double>> noweX(this->rozmiar, { 0.0, 0.0 });

    while (true) {
        for (int i = 0; i < this->rozmiar; ++i) {
            std::complex<double> sumaa = { 0.0, 0.0 };
            for (int j = 0; j < this->rozmiar; ++j) {
                if (i != j) {
                    sumaa += this->macierzA[i][j] * X[j];
                }
            }
            noweX[i] = (this->wektorB[i] - sumaa) / this->macierzA[i][i];
        }

        double maxError = 0.0;
        for (int i = 0; i < this->rozmiar; ++i) {
            maxError = std::max(maxError, std::abs(noweX[i] - X[i]));
        }

        if (maxError < this->tolerancja) {
            return noweX;
        }
        X = noweX;
    }
}


std::vector<std::complex<double>> MetodaJacobiego::obliczOpenMP() {
    std::vector<std::complex<double>> X(this->rozmiar, { 0.0, 0.0 });
    std::vector<std::complex<double>> noweX(this->rozmiar, { 0.0, 0.0 });

    while (true) {
        #pragma omp parallel for
        for (int i = 0; i < this->rozmiar; ++i) {
            std::complex<double> suma = { 0.0, 0.0 };
            for (int j = 0; j < this->rozmiar; ++j) {
                if (i != j) {
                    suma += macierzA[i][j] * X[j];
                }
            }
            noweX[i] = (wektorB[i] - suma) / macierzA[i][i];
        }

        double maxError = 0.0;
        for (int i = 0; i < this->rozmiar; ++i) {
            maxError = std::max(maxError, std::abs(noweX[i] - X[i]));
        }

        if (maxError < this->tolerancja) {
            return noweX;
        }
        X = noweX;
    }
}




void MetodaJacobiego::policzRownolegle(
    int start, int koniec, std::vector<std::complex<double>>& x, 
    std::vector<std::complex<double>>& noweX
) {
    for (int i = start; i < koniec; ++i) {
        std::complex<double> suma = { 0.0, 0.0 };
        for (int j = 0; j < this->rozmiar; ++j) {
            if (i != j) {
                suma += this->macierzA[i][j] * x[j];
            }
        }
        noweX[i] = (this->wektorB[i] - suma) / this->macierzA[i][i];
    }
}

std::vector<std::complex<double>> MetodaJacobiego::obliczWatki() {
    std::vector<std::complex<double>> x(this->rozmiar, { 0.0, 0.0 });
    std::vector<std::complex<double>> noweX(this->rozmiar, { 0.0, 0.0 });

    std::vector<std::thread> threads;

    while (true) {
        threads.clear();

        for (int i = 0; i < this->liczbaWatkow; ++i) {
            int start = i * (this->rozmiar / this->liczbaWatkow);
            int koniec = (i == this->liczbaWatkow - 1) ? this->rozmiar : (i + 1) * (this->rozmiar / this->liczbaWatkow);
            threads.push_back(std::thread(&MetodaJacobiego::policzRownolegle, this, start, koniec, std::ref(x), std::ref(noweX)));
        }

        for (auto& t : threads) {
            t.join();
        }

        double maxError = 0.0;
        for (int i = 0; i < this->rozmiar; ++i) {
            maxError = std::max(maxError, std::abs(noweX[i] - x[i]));
        }

        if (maxError < this->tolerancja) {
            return noweX;
        }

        x = noweX;
    }
}

bool MetodaJacobiego::porownajWektory(const std::vector<std::complex<double>>& v1, const std::vector<std::complex<double>>& v2) {
    double epsilon = this->tolerancja;

    if (v1.size() != v2.size()) {
        return false;
    }

    for (size_t i = 0; i < v1.size(); ++i) {
        if (std::abs(v1[i] - v2[i]) > epsilon) {
            return false;
        }
    }

    return true;
}


