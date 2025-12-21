package lab__2;

import lab1.Differentiation;
import lab1.MatrixSystem;
import lab1.TridiagonalMatrixCalculator;

import java.util.function.BinaryOperator;
import java.util.function.UnaryOperator;

public class EulerSolver {
    int n;
    int timeSteps;
    double a;
    double b;
    double tend;
    BinaryOperator<Double> u;
    BinaryOperator<Double> k;
    BinaryOperator<Double> q;
    double xi1;
    UnaryOperator<Double> phi;
    double maxError;
    double maxNev;
    public double[] countSourceSingleMatrix(double t) {
        double[] f = new double[n];
        double r = a;
        double h = (b - a) / n;
        for (int i = 0; i < n; i++) {
            double dk = Differentiation.derR(k, r, t);
            double du = Differentiation.derR(u, r, t);
            double d2u = Differentiation.der2R(u, r, t);
            double firstTerm = dk * du + k.apply(r, t) * du / r + k.apply(r, t) * d2u;
            double initialRightSide = firstTerm - q.apply(r, t) * u.apply(r, t) // + f.apply(r, t)
            ;
            f[i] = Differentiation.derivativeByT(u, r, t) -
                    initialRightSide;
            r += h;
        }
        return f;
    }

    public MatrixSystem differenceScheme(double t, double[] f) {
        double h = (b - a) / n;
        double hh = h / 2.0;
        double r = a;
        double[] aa = new double[n];
        double[] bb = new double[n];
        double[] cc = new double[n];
        double[] ff = new double[n];

        // для i = 0
        aa[0] = 0;
        bb[0] =-(
                (r + hh) * kmid(r, r+h, t) / (h * hh * r)
                        + q.apply(r, t)
                        + xi1 / hh
        );
        cc[0] = (r + hh) * kmid(r, r+h, t) / (h * hh * r);
        ff[0] = f[0] + v1(t) / hh;

        // для i от 1 до N-2
        for (int i = 1; i < n - 1; i++)
        {
            r += h;
            aa[i] = (r - hh) * kmid(r-h, r, t) / (h * h * r);
            bb[i] = -(
                    ((r + hh) * kmid(r, r+h, t)) / (h * h *r)
                            + ((r - hh) * kmid(r, r-h, t)) / (h * h *r)
                            + q.apply(r, t)
            );
            cc[i] = (r + hh) * kmid(r, r+h, t) / (h * h *r) ;
            ff[i] = f[i];
        }

        // для i = N-1
        r += h;
        aa[n - 1] = (r - hh) * kmid(r, r-h, t) / (h*h*r);
        bb[n - 1] =-(
                (r - hh) * kmid(r, r-h, t) / (h * h * r)
                        + (r + hh) * kmid(r, r+h, t) / (h * h * r)
                        + q.apply(r, t)
        );
        cc[n - 1] = 0.0;
        ff[n - 1] = f[n - 1] + (r + hh) * kmid(r, r+h, t) / (h*h*r) * v2(t);

        return new MatrixSystem(t, aa, bb, cc, ff);
    }

    public void implicit() {
        double currentTime = 0;
        double tau = tend / timeSteps;
        double[][] uu = new double[timeSteps+1][n];
        double hr = (b - a) / n;
        for (int i = 0; i < n; i++) {
            uu[0][i] = phi.apply(hr*i + a);
        }
        for (int j = 0; j < timeSteps; j++) {
            currentTime = tau * (j+1);
            double[] g = countSourceSingleMatrix(currentTime);
            MatrixSystem source = differenceScheme(currentTime, g);
            double[] rightSide = new double[n];
            for (int i = 0; i < n; i++) {
                rightSide[i] = uu[j][i] + tau * source.f[i];
            }
            MatrixSystem next = new MatrixSystem(rightSide);
            for (int i = 0; i < n; i++) {
                next.a[i] = - tau * source.a[i];
                next.b[i] = - tau * source.b[i] + 1;
                next.c[i] = - tau * source.c[i];
            }
            uu[j+1] = TridiagonalMatrixCalculator.solveTridiagonal(next.a, next.b, next.c, next.f);
        }
        double[][] nevyzaki = new double[timeSteps+1][n];
        double[][] err = new double[timeSteps+1][n];
        for (int j = 0; j <= timeSteps; j++) {
            double[] x = new double[n];
            for (int i = 0; i < n; i++) {
                x[i] = u.apply(hr*i+a, j*tau);
            }
            MatrixSystem matrix = differenceScheme(j*tau, countSourceSingleMatrix(j*tau));
            double[] Av = TridiagonalMatrixCalculator.multiplyTridiagonal(matrix.a, matrix.b, matrix.c, x);
            for (int i = 0; i < n; i++) {
                err[j][i] = u.apply(hr*i+a, j*tau) - uu[j][i];
                nevyzaki[j][i] = Differentiation.derivativeByT(u, hr*i+a, j*tau)-Av[i]-matrix.f[i];
            }
        }
        finish(err, nevyzaki);
    }

    public void explicit() {
        double currentTime = 0;
        double tau = tend / timeSteps;
        double[][] uu = new double[timeSteps+1][n];
        double hr = (b - a) / n;
        for (int i = 0; i < n; i++) {
            uu[0][i] = phi.apply(hr*i + a);
        }
        for (int j = 0; j < timeSteps; j++) {
            currentTime = tau * j;
            double[] g = countSourceSingleMatrix(currentTime);
            MatrixSystem matrix = differenceScheme(currentTime, g);
            for (int i = 0; i < n; i++) {
                matrix.a[i] = tau * matrix.a[i];
                matrix.b[i] = matrix.b[i] * tau + 1;
                matrix.c[i] = tau * matrix.c[i];
            }
            double[] values = TridiagonalMatrixCalculator.multiplyTridiagonal(matrix.a, matrix.b, matrix.c, uu[j]);
            for (int i = 0; i < n; i++) {
                values[i] = values[i] + tau * matrix.f[i];
            }
            uu[j+1] = values;
        }
        double[][] nevyzaki = new double[timeSteps+1][n];
        double[][] err = new double[timeSteps+1][n];
        for (int j = 0; j <= timeSteps; j++) {
            double[] x = new double[n];
            for (int i = 0; i < n; i++) {
                x[i] = u.apply(hr*i+a, j*tau);
            }
            MatrixSystem matrix = differenceScheme(j*tau, countSourceSingleMatrix(j*tau));
            double[] Av = TridiagonalMatrixCalculator.multiplyTridiagonal(matrix.a, matrix.b, matrix.c, x);
            for (int i = 0; i < n; i++) {
                err[j][i] = Math.abs(u.apply(hr*i+a, j*tau) - uu[j][i]);
                nevyzaki[j][i] = Math.abs(Differentiation.derivativeByT(u, hr*i+a, j*tau)-Av[i]-matrix.f[i]);
            }
        }
        finish(err, nevyzaki);
    }


    private double kmid(double i, double i_plus, double t) {
        return (k.apply(i_plus, t) + k.apply(i, t)) / 2;
    }

    public EulerSolver(int n, int timeSteps, double a, double b, double tend, BinaryOperator<Double> u, BinaryOperator<Double> k, BinaryOperator<Double> q, double xi1) {
        this.n = n;
        this.timeSteps = timeSteps;
        this.a = a;
        this.b = b;
        this.tend = tend;
        this.u = u;
        this.k = k;
        this.q = q;
        this.xi1 = xi1;
        this.phi = r -> u.apply(r, 0d);
    }

    private double v2(double t) {
        return u.apply(b, t);
    }

    private double v1(double t) {
        return xi1 * u.apply(a, t) - k.apply(a, t) * Differentiation.derR(u, a, t);
    }

    public void printMatrix(double[][] matrix) {
        for (int j = 0; j < matrix.length; j++) {
            for (int i = 0; i < matrix[0].length; i++) {
                System.out.printf("%.3e\t", matrix[j][i]);
            }
            System.out.println();
        }
        System.out.println();
    }

    public double findMax(double[][] matrix) {
        double max = 0;
        for (int j = 0; j < matrix.length; j++) {
            for (int i = 0; i < matrix[0].length; i++) {
                if (matrix[j][i] > max) max = matrix[j][i];
            }
        }
        return max;
    }

    private void finish(double[][] err, double[][] nevyzaki) {
        //printMatrix(err);
        //printMatrix(nevyzaki);
        maxError = findMax(err);
        maxNev = findMax(nevyzaki);
        //System.out.printf("Макс. по модулю погрешность: %.3e\n", maxError);
        //System.out.printf("Макс. по модулю невязка: %.3e\n", maxNev);
    }
}
