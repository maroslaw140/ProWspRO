#ifndef METODAJACOBIEGO_H
#define METODAJACOBIEGO_H

#include <vector>
#include <complex>

class MetodaJacobiego {
private:
    int rozmiar = 40000;
    double minLimit;
    double maxLimit;

    double tolerancja = 1e-15;
    int liczbaWatkow = 4;

    std::vector<std::vector<std::complex<double>>> macierzA;
    std::vector<std::complex<double>> wektorB;

public:
    MetodaJacobiego();
    void generujMacierz();
    std::vector<std::complex<double>> oblicz();
    std::vector<std::complex<double>> obliczOpenMP();

	void policzRownolegle(int start, int end, std::vector<std::complex<double>>& x, std::vector<std::complex<double>>& x_new);
	std::vector<std::complex<double>> obliczWatki();

    bool porownajWektory(const std::vector<std::complex<double>>& v1, const std::vector<std::complex<double>>& v2);

    int getWatki();
	void setWatki(int liczbaWatki);

	int getRozmiar();
	void setRozmiar(int rozmiar);
};

#endif
