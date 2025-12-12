package lab1;

import java.util.function.UnaryOperator;

public class Runner {
    int n;
    double rL;
    double rR;
    double xi;
    double nu1;
    double nu2;
    UnaryOperator<Double> k;
    UnaryOperator<Double> u;
    UnaryOperator<Double> q;
    double maxError;
    double maxN;
    double innerMaxN;

    public Runner(
            int n,
            double rL,
            double rR,
            double xi,
            UnaryOperator<Double> k,
            UnaryOperator<Double> u,
            UnaryOperator<Double> q) {
        this.n = n; this.rL = rL; this.rR = rR; this.xi = xi;
        this.k = k; this.q = q;
        this.u = u;
        this.nu1 = xi * u.apply(rL) - k.apply(rL) * Differentiation.derivative(u, rL);
        this.nu2 = u.apply(rR);
    }

    public void clear() {
        maxError = 0;
        maxN = 0;
        innerMaxN = 0;
    }
    public void run() {
        System.out.printf("nu1 = %.1f, nu2 = %.1f\n", nu1, nu2);
        countValues();
        nevyazka();
        System.out.printf("Макс. погрешность %.3e\n", maxError);
        System.out.printf("Макс. невязка %.3e\n", maxN);
        System.out.printf("Макс. невязка без учета граничной точки %.3e\n", innerMaxN);
        System.out.print("\n\n");
    }

    public void nevyazka() {
        System.out.println("Невязки:");
        double[] x = new double[n];
        double current = rL;
        double h = rR - rL;
        h /= n;
        for (int i = 0; i < n; i++) {
            x[i] = u.apply(current);
            current += h;
        }
        MatrixSystem matrix = DiffSeqSolver.differenceScheme(rL, rR, n, nu1, nu2, xi,
                k, u, q);
        double[] Av = TridiagonalMatrixCalculator.multiplyTridiagonal(matrix.a, matrix.b, matrix.c, x);
        for (int i = 0; i < n; i++) {
            double res = matrix.f[i] - Av[i];
            if (Math.abs(res) > Math.abs(maxN)) maxN = res;
            if (i != 0) {
                if (Math.abs(res) > Math.abs(innerMaxN)) innerMaxN = res;
            }
            System.out.printf("%.3e\n", res);
        }
    }

    public void countValues() {
        MatrixSystem matrix = DiffSeqSolver.differenceScheme(rL, rR, n, nu1, nu2, xi,
                k, u, q);
        matrix.print();
        double[] x = TridiagonalMatrixCalculator.solveTridiagonal(matrix.a, matrix.b, matrix.c, matrix.f);
        double current = rL;
        double h = rR - rL;
        h /= n;
        for (int i = 0; i < n; i++) {
            System.out.printf("x[%d] = %6.2e", i, x[i]);
            double error = Math.abs(x[i] - u.apply(current));
            if (Math.abs(error) > Math.abs(maxError)) maxError = error;
            System.out.printf(" (погрешность: %.3e)\n", error);
            current += h;
        }
        System.out.println();
    }
}