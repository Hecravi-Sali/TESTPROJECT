#include "Code.h"
namespace C2 {
	class Newdon_it::_Data {
	public:
		std::vector<Size2D> _array;
		std::vector<double> _ans;
		void calcaute()
		{
			auto pointsize = _array.size();
			_ans.clear();
			for (auto i = 0; i < pointsize; i++)
				_ans.push_back(_array[i].second);
			for (auto i = 0; i < pointsize; i++)
				for (int j = pointsize - 1; j > i; j--)
					_ans[j] = (_ans[j] - _ans[j - 1]) / (_array[j].first - _array[j - 1 - i].first);
#if _MYDEBUG_
			for each(auto i in _ans)
				std::cout << i << "  ";
			std::cout << std::endl;
#endif
		}
	};
	Newdon_it::Newdon_it(std::vector<Size2D>&inputpoint) {
		_data = new _Data;
		_data->_array = inputpoint;
		_data->calcaute();
	}
	Newdon_it::~Newdon_it(void) {
		delete _data;
	}
	double Newdon_it::res(const double & xt, const int & _order)
	{
		double Res = _data->_ans[0];
		double iteration_part = 1;
		for (auto i = 1; i<_order; i++) {
			iteration_part *= (xt - _data->_array[i - 1].first);
			Res += _data->_ans[i] * iteration_part;
		}
		return Res;
	}
	std::vector<double> Newdon_it::GetAns(void) const {
		return _data->_ans;
	}

	class LongBerg::_Data {
	public:
		double coutleanswer(const double & error) {
			auto coutlenext = [=](const double & Tn, const double & T2n, const int & m) {
				return T2n + 1 / (pow(4, m) - 1) * (T2n - Tn);
			};
			auto coutlet2n = [=](const double & Tn, const int64_t & n) {
				auto truen = pow(2, n - 1);
				double step = _size / truen;
				double interval = step * 0.5;
				double xx = 0.;
				double point = _start + interval;
				for (auto i = 0; i < truen; i++, point += step)
					xx += _function(point);
				return Tn * 0.5 + interval*xx;
			};
			double up = 0., pas = 0.;
			while (true) {
				for (auto i = 0; i < _Tn.size(); i++) {
					pas = _Tn[i];
					if (i == 0)
						_Tn[i] = coutlet2n(_Tn[0], _Tn.size());
					else
						_Tn[i] = coutlenext(up, _Tn[i - 1], i);
					up = pas;
				}
#if _MYDEBUG_
				for each(auto i in _Tn)
					std::cout << i << "  ";
				std::cout << std::endl;
#endif
				pas = coutlenext(up, _Tn[_Tn.size() - 1], _Tn.size());
				if (abs(up - _Tn[_Tn.size() - 1]) < error) {
					return pas;
				}
				else
					_Tn.push_back(pas);
			}
		}
		float _start;
		float _size;
		std::function<double(double)>_function;
		std::vector<double>_Tn;
	};
	LongBerg::LongBerg(void) {
		_data = new _Data;
	}
	LongBerg::~LongBerg(void) {
		delete _data;
	}
	double LongBerg::Integral(std::function<double(const double &)> && infunction, const double && a, const double &&b, const double && error) {
		_data->_start = a;
		_data->_size = b - a;
		_data->_function = infunction;
		_data->_Tn.push_back(_data->_size*(infunction(a) + infunction(b))*0.5);
		return _data->coutleanswer(error);
	}

	class GaissElimination::_Data {
	public:
#if _MYDEBUG_
		void show(double **c, const int &size) {
			for (auto i = 0; i < size; i++) {
				for (auto k = 0; k < size + 1; k++)
					std::cout << c[i][k] << "      ";
				std::cout << std::endl;
			}
			std::cout << std::endl;
		}
#endif
		inline void huanhang(const Size2I &begin, const int &move, const int &size, double **matrix) {
			double pas = 0.;
			for (auto i = begin.second; i < size + 1; i++) {
				pas = matrix[begin.first][i];
				matrix[begin.first][i] = matrix[move][i];
				matrix[move][i] = pas;
			}
		}
		inline void huanlie(const Size2I &begin, const int &move, const int &size, double **matrix) {
			double pas = 0.;
			for (auto i = 0; i < size; i++) {
				pas = matrix[i][begin.second];
				matrix[i][begin.second] = matrix[i][move];
				matrix[i][move] = pas;
			}
		}
		inline Size2I max(const Size2I &begin, const int &size, double **matrix) {
			double max = 0.;
			Size2I point;
			for (auto i = begin.first; i < size; i++)
				for (auto k = begin.second; k < size; k++)
					if (fabs(matrix[i][k]) > max) {
						max = fabs(matrix[i][k]);
						point.first = i;
						point.second = k;
					}
			return point;
		}
		inline void minusline(const int &begin, const int &size, double **matrix) {
			for (auto i = begin + 1; i < size; i++) {
				double cc = static_cast<double>(matrix[i][begin]) / static_cast<double>(matrix[begin][begin]);
				cc = -cc;
				matrix[i][begin] -= matrix[i][begin];
				for (auto k = begin + 1; k < size + 1; k++)
					matrix[i][k] += matrix[begin][k] * cc;
			}
		}

		std::vector<double> sca(double **matrix_x, const int &size) {
			std::vector<double> ans;
#if GaissEliminationType
			int *oo = new int[size];
#endif
			for (auto i = 0; i < size; i++) {
#if GaissEliminationType
				oo[i] = i;
#endif
				ans.push_back(0.);
			}
			for (auto i = 0; i < size - 1; i++) {
				Size2I &nowpoint = Size2I(i, i);
				Size2I &maxpoint = max(nowpoint, size, matrix_x);
				if (maxpoint.first != nowpoint.first) {
					huanhang(nowpoint, maxpoint.first, size, matrix_x);
#if _MYDEBUG_
					std::cout << std::endl;
					show(matrix_x, size);
					std::cout << std::endl;
#endif
				}
#if GaissEliminationType
				if (maxpoint.second != nowpoint.second) {
					huanlie(nowpoint, maxpoint.second, size, matrix_x);
					int pas = oo[nowpoint.second];
					oo[nowpoint.second] = oo[maxpoint.second];
					oo[maxpoint.second] = pas;
#if _MYDEBUG_
					std::cout << "----jioahuanlie----" << std::endl;
					show(matrix_x, size);
					std::cout << std::endl;
#endif
				}
#endif
				if (nowpoint.first == nowpoint.second == 0.)
					return ans;
				minusline(i, size, matrix_x);
#if _MYDEBUG_
				std::cout << std::endl;
				for (auto i = 0; i < size; i++)
#if GaissEliminationType

					std::cout << oo[i] << " ";
#else
					std::cout << i << " ";
#endif
				std::cout << std::endl;
				show(matrix_x, size);
				system("pause");
#endif
			}
#if GaissEliminationType
			for (auto i = size - 1; i >= 0; i--) {
				double sum = 0.;
				for (auto k = i + 1; k < size; k++)
					sum += matrix_x[i][k] * ans[oo[k]];
				ans[oo[i]] = (matrix_x[i][size] - sum) / matrix_x[i][i];
			}
#else
			for (auto i = size - 1; i >= 0; i--) {
				double sum = 0.;
				for (auto k = i + 1; k < size; k++)
					sum += matrix_x[i][k] * ans[k];
				ans[i] = (matrix_x[i][size] - sum) / matrix_x[i][i];
			}
#endif
			return ans;
		}
	};
	GaissElimination::GaissElimination(void) {
		_data = new _Data;
	}
	GaissElimination::~GaissElimination(void) {
		delete _data;
	}
	std::vector<double> GaissElimination::Calculation(const int &size, double **matrix)const
	{
		double **mat = nullptr;
		mat = new double*[size];
		for (auto i = 0; i < size; i++) {
			mat[i] = new double[size + 1];
			for (auto k = 0; k < size + 1; k++)
				mat[i][k] = matrix[i][k];
		}
		std::vector<double>& ans = _data->sca(mat, size);
#if _MYDEBUG_
		for each(auto i in ans)
			std::cout << i << std::endl;
		for (auto i = 0; i < size; i++) {
			double sum = 0.;
			for (auto k = 0; k < size; k++)
				sum += matrix[i][k] * ans[k];
			std::cout << sum << std::endl;
		}
#endif
		for (auto i = 0; i < size; i++)
			delete mat[i];
		delete mat;
		return ans;
	}
}
