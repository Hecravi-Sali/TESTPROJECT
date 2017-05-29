#include "Code.h"
using namespace std;
using namespace C2;
void NewdonTest(void) {
	int num_of_inter_point; 
	cout << "输入插值点个数: ";
	cin >> num_of_inter_point;
	std::vector<Size2D> data;
	for (int i = 0; i<num_of_inter_point; i++){
		double a1, a2;
		cin >> a1 >> a2;
		data.push_back(Size2D(a1, a2));
	}
	Newdon_it Ni(data);
	double xt;
	int Order;
	while (true) {
		cin >> xt;
		cin >> Order;
		printf("%.21f", Ni.res(xt, Order));
	}
}
void LoneTest(void) {
	LongBerg a;
	auto f = [=](const double & x)->double {
		if (x == 0)
			return 1;
		return sin(x) / x; };
	cout.precision(16);
	cout << a.Integral(f, 0, 3, 1e-6) << endl;
 }
void GaissEliminationTest(void) {
	int num = 0;
	cin >> num;
	double **c = nullptr;
	c = new double*[num];
	for (auto i = 0; i < num; i++) {
		c[i] = new double[num + 1];
		for (auto k = 0; k < num + 1; k++)
			cin >> c[i][k];
	}
	GaissElimination count;
	vector<double>& ans = count.Calculation(num, c);
	for each(auto i in ans)
		cout << i << "  ";
	cout << endl;
	/*
	3
	-0.002 2 2 0.4
	1 0.78125 0 1.3816
	3.996 5.5625 4 7.4178

	4
	1	   0.333    1.5   -0.333  3
	-2.01  1.45     0.5   2.95    5.4
	4.32   -1.95    0     2.08    0.13
	5.11   -4       3.33  -1.11   3.77
	*//*
	3
	1 1 1  3
	0 4 -1 5
	2 -2 1 2

	3
	1 2 3 14
	2 5 2 18
	3 1 5 20

	*/
}
int main(int argc, char *argv[]) {
#if false
	NewdonTest();
#endif
#if false
	LoneTest();
#endif
#if true
	GaissEliminationTest();
#endif
	system("pause");
	return 0;
}

