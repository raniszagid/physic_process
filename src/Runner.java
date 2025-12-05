import java.util.function.UnaryOperator;

public class Runner {
    int n;
    double rL;
    double rR;
    double xi1;
    double xi2;
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
            double xi1,
            double xi2,
            UnaryOperator<Double> k,
            UnaryOperator<Double> u,
            UnaryOperator<Double> q) {
        this.n = n; this.rL = rL; this.rR = rR; this.xi1 = xi1; this.xi2 = xi2;
        this.k = k; this.q = q;
        this.u = u;
        this.nu1 = xi1 * u.apply(rL) - k.apply(rL) * Differentiation.derivative(u, rL, 1e-10);
        this.nu2 = xi2 * u.apply(rR) + k.apply(rR) * Differentiation.derivative(u, rR, 1e-10);
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
        double[] x = new double[n+1];
        double current = rL;
        double h = rR - rL;
        h /= n;
        for (int i = 0; i <= n; i++) {
            x[i] = u.apply(current);
            current += h;
        }
        MatrixSystem matrix = DiffSeqSolver.differenceScheme(rL, rR, n, nu1, nu2, xi1, xi2,
                k, u, q);
        double[] Av = TridiagonalMatrixCalculator.multiplyTridiagonal(matrix.a, matrix.b, matrix.c, x);
        for (int i = 0; i <= n; i++) {
            double res = matrix.f[i] - Av[i];
            if (Math.abs(res) > Math.abs(maxN)) maxN = res;
            if (i > 0 && i < n-1) {
                if (Math.abs(res) > Math.abs(innerMaxN)) innerMaxN = res;
            }
            System.out.printf("%.3e\n", res);
        }
    }

    public void countValues() {
        MatrixSystem matrix = DiffSeqSolver.differenceScheme(rL, rR, n, nu1, nu2, xi1, xi2,
                k, u, q);
        double[] x = TridiagonalMatrixCalculator.solveTridiagonal(matrix.a, matrix.b, matrix.c, matrix.f);
        double current = rL;
        double h = rR - rL;
        h /= n;
        for (int i = 0; i <= n; i++) {
            System.out.printf("x[%d] = %6.2e", i, x[i]);
            double error = Math.abs(x[i] - u.apply(current));
            if (Math.abs(error) > Math.abs(maxError)) maxError = error;
            System.out.printf(" (погрешность: %.3e)\n", error);
            current += h;
        }
        System.out.println();
    }
}