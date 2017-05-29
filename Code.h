#pragma once
#ifndef _C2_creodfactore
#define _C2_creodfactore
#include <vector>
#include <cmath>
#include <iostream>
#include <functional>

#define _MYDEBUG_ false

namespace C2 {
	typedef std::pair <double, double> Size2D;
	typedef std::pair <int, int> Size2I;

	class Newdon_it
	{
	public:
		Newdon_it(std::vector<Size2D>&inputpoint);
		~Newdon_it(void);
		double res(const double & xt, const int & _order = 1);
		std::vector<double> GetAns(void)const;
	private:
		class _Data;
		_Data *_data;
	};
	class LongBerg {
	public:
		LongBerg();
		~LongBerg();
		double Integral(std::function<double(const double &)> && infunction, const double && a, const double &&b, const double && error = 1e-6);
	private:
		class _Data;
		_Data *_data;
	};
	//  启动为全主元高斯消去
	//  关闭为列主元高斯消去
#define GaissEliminationType false
	class GaissElimination {
	public:
		GaissElimination(void);
		~GaissElimination(void);
		std::vector<double> Calculation(const int &size, double **matrix)const;
	private:
		class _Data;
		_Data *_data;
	};
	class SOR {
	public:

	private:

	};
}
#endif // !_C2_creodfactore
