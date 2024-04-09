#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

double f(vector<double> x) {
    return 2 * x[0] * x[0] - 4 * x[0] * x[1] + 3 * x[1] * x[1] + x[0] - x[1];
}

vector<double> gradient(vector<double> x) {
    return {4 * x[0] - 4 * x[1] + 1, -4 * x[0] + 6 * x[1] - 1};
}

vector<vector<double>> hessian() {
    return {{4, -4}, {-4, 6}};
}

vector<double> next_s(vector<double> grad1, vector<double> grad2, vector<double> s) {
    return {-grad2[0] + (pow(grad2[0]*grad2[0] + grad2[1]*grad2[1], 2) / pow(grad1[0]*grad1[0] + grad1[1]*grad1[1], 2)) * s[0],
            -grad2[1] + (pow(grad2[0]*grad2[0] + grad2[1]*grad2[1], 2) / pow(grad1[0]*grad1[0] + grad1[1]*grad1[1], 2)) * s[1]};
}

double lambda(vector<double> grad, vector<double> s, vector<vector<double>> hess) {
    double up = grad[0] * grad[0] + grad[1] * grad[1];
    double down = (hess[0][0] * s[0] + hess[0][1] * s[1]) * s[0] + (hess[1][0] * s[0] + hess[1][1] * s[1]) * s[1];
    return up/down;
}

vector<double> fletcher_reeves(vector<double> x, double eps) {
    vector<double> grad = gradient(x);
    vector<double> prev_grad;
    vector<vector<double>> hess = hessian();
    vector<double> s = {-grad[0], -grad[1]};
    int iter = 1;

    while (sqrt(grad[0] * grad[0] + grad[1] * grad[1]) >= eps) {
        double lam = lambda(grad, s, hess);
        x = {x[0] + lam * s[0], x[1] + lam * s[1]};

        cout << "Крок №" << iter++ << ":" << endl <<
            "delf = [" << grad[0] << ", " << grad[1] << "]" << endl <<
            "||delf|| = " << sqrt(grad[0]*grad[0] + grad[1]*grad[1]) << endl <<
            "s = [" << s[0] << ", " << s[1] << "]"  << endl <<
            "lam = " << lam << endl <<
            "x" << iter << " = [" << x[0] << ", " << x[1] << "]" << endl << endl;

        prev_grad = grad;
        grad = gradient(x);
        s = next_s(prev_grad, grad, s);
    }

    return x;
}

int main() {
    vector<double> x0 = {5, 2};
    double eps = 1e-2;

    vector<double> solution = fletcher_reeves(x0, eps);

    cout << "Minimum point: x = (" << fixed << setprecision(3) << solution[0] << ", " << solution[1] << "), Minimum value: " << f(solution) << endl;
}
