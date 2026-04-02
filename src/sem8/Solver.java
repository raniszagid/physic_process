package sem8;

import lab1.Differentiation;

import java.util.function.BinaryOperator;
import java.util.function.UnaryOperator;

public class Solver {
    /*
    i+1, j = a
    i-1, j = b
    i, j   = c
    i, j+1 = d
    i, j-1 = e
    */
    BinaryOperator<Double> k1;
    BinaryOperator<Double> k2;
    BinaryOperator<Double> u;
    UnaryOperator<Double> g1;
    UnaryOperator<Double> g2;
    UnaryOperator<Double> g3;
    UnaryOperator<Double> g4;
    double xi;
    double A;
    double B;
    double C;
    double D;
    int N;
    int M;
    double hx;
    double hy;

    public Solver(BinaryOperator<Double> k1, BinaryOperator<Double> k2,
                  double A, double B, double C, double D, int N, int M,
                  BinaryOperator<Double> u,
                  UnaryOperator<Double> g1, UnaryOperator<Double> g2,
                  UnaryOperator<Double> g3, UnaryOperator<Double> g4, double xi) {
        this.k1 = k1;
        this.k2 = k2;
        this.A = A;
        this.B = B;
        this.C = C;
        this.D = D;
        this.N = N;
        this.M = M;
        this.u = u;
        this.g1 = g1;
        this.g2 = g2;
        this.g3 = g3;
        this.g4 = g4;
        this.xi = xi;
        this.hx = (B - A) / N;
        this.hy = (D - C) / M;
    }

    public void run() {
        MatrixManager calc = new MatrixManager();
        //Main.printMatrix(countSourceMatrix());
        //System.out.println();
        DiffScheme diffScheme = differenceScheme();
        diffScheme.print();
        MatrixManager.TridiagonalSystem system = calc.convertToTridiagonal(
                diffScheme.a, diffScheme.b, diffScheme.c,
                diffScheme.d, diffScheme.e, diffScheme.f, N, M
        );

        // Решение методом матричной прогонки
        double[][] xyCounted = calc.matrixProgonka(
                system.a, system.b, system.c, system.f,
                system.sizeY, system.sizeX
        );
        double x = A; double y = C;
        for (int j = 0; j < M; j++) {
            y = C + j * hy;
            for (int i = 0; i < N; i++) {
                x = A + i * hx;
                double expected = u.apply(x, y);
                double actual = xyCounted[j][i];
                System.out.printf("(%d, %d)\t%.3f\t%.3f\n", i, j, expected, actual);
            }
        }
    }

    public DiffScheme differenceScheme() {
        double x = A; double y = C;
        double[][] a = new double[M][N];
        double[][] b = new double[M][N];
        double[][] c = new double[M][N];
        double[][] d = new double[M][N];
        double[][] e = new double[M][N];
        double[][] f = new double[M][N];
        double[][] source = countSourceMatrix();

        for (int j = 0; j < M; j++) {
            y = C + j * hy;
            if (j == 0) {
                for (int i = 0; i < N - 1; i++) {
                    x = A + i * hx;
                    if (i == 0) {
                        c[j][i] = 1;
                        a[j][i] = 0;
                        b[j][i] = 0;
                        d[j][i] = 0;
                        e[j][i] = 0;
                        f[j][i] = g1.apply(y);
                    }
                    else {
                        c[j][i] = hy / hx / 2 * kmidx(k1, x, x + hx, y)
                                + hy / hx / 2 * kmidx(k1, x, x - hx, y)
                                + hx / hy * kmidy(k2, x, y, y + hy)
                                + hx * xi;
                        d[j][i] = -hx / hy * kmidy(k2, x, y, y + hy);
                        e[j][i] = 0;
                        if (i == 1) {
                            // для i = 1, j = 0
                            a[j][i] = -hy / hx / 2 * kmidx(k1, x, x + hx, y);
                            b[j][i] = 0;
                            f[j][i] = hx * hy / 2 * source[j][i]
                                    + hx * g3.apply(x)
                                    + hy / hx / 2 * kmidx(k1, x, x - hx, y) * g1.apply(y);
                        } else if (i == N - 1) {
                            // для i = N - 1, j = 0
                            a[j][i] = 0;
                            b[j][i] = -hy / hx * kmidx(k1, x, x - hx, y);
                            f[j][i] = hx * hy / 2 * source[j][i]
                                    + hx * g3.apply(x)
                                    + hy / hx / 2 * kmidx(k1, x, x + hx, y) * g2.apply(y);
                        } else {
                            // для i = 2 ... N - 2, j = 0
                            a[j][i] = -hy / hx / 2 * kmidx(k1, x, x + hx, y);
                            b[j][i] = -hy / hx / 2 * kmidx(k1, x, x - hx, y);
                            f[j][i] = hx * hy / 2 * source[j][i]
                                    + hx * g3.apply(x);
                        }
                    }
                }
            }
            else if (j == M - 1) {
                for (int i = 0; i < N - 1; i++) {
                    x = A + i * hx;
                    if (i == 0) {
                        c[j][i] = 1;
                        a[j][i] = 0;
                        b[j][i] = 0;
                        d[j][i] = 0;
                        e[j][i] = 0;
                        f[j][i] = g1.apply(y);
                    }
                    else {
                        c[j][i] = hy / hx * kmidx(k1, x, x + hx, y)
                                + hy / hx * kmidx(k1, x, x - hx, y)
                                + hx / hy * kmidy(k2, x, y, y + hy)
                                + hx / hy * kmidy(k2, x, y, y - hy);
                        d[j][i] = 0;
                        e[j][i] = -hx / hy * kmidy(k2, x, y, y - hy);
                        if (i == 1) {
                            // для i = 1, j = M - 1
                            a[j][i] = -hy / hx * kmidx(k1, x, x + hx, y);
                            b[j][i] = 0;
                            f[j][i] = hx * hy * source[j][i]
                                    + hy / hx * kmidx(k1, x, x - hx, y) * g1.apply(y)
                                    + hx / hy * kmidy(k2, x, y, y + hy) * g4.apply(x);
                        } else if (i == N - 1) {
                            // для i = N - 1, j = M - 1
                            a[j][i] = 0;
                            b[j][i] = -hy / hx * kmidx(k1, x, x - hx, y);
                            f[j][i] = hx * hy * source[j][i]
                                    + hy / hx * kmidx(k1, x, x + hx, y) * g2.apply(y)
                                    + hx / hy * kmidy(k2, x, y, y + hy) * g4.apply(x);
                        } else {
                            // для i = 2 ... N - 2, j = M - 1
                            a[j][i] = -hy / hx * kmidx(k1, x, x + hx, y);
                            b[j][i] = -hy / hx * kmidx(k1, x, x - hx, y);
                            f[j][i] = hx * hy * source[j][i]
                                    + hx / hy * kmidy(k2, x, y, y + hy) * g4.apply(x);
                        }
                    }
                }
            }
            else {
                for (int i = 0; i < N - 1; i++) {
                    x = A + i * hx;
                    if (i == 0) {
                        c[j][i] = 1;
                        a[j][i] = 0;
                        b[j][i] = 0;
                        d[j][i] = 0;
                        e[j][i] = 0;
                        f[j][i] = g1.apply(y);
                    }
                    else {
                        c[j][i] = hy / hx * kmidx(k1, x, x + hx, y)
                                + hy / hx * kmidx(k1, x, x - hx, y)
                                + hx / hy * kmidy(k2, x, y, y + hy)
                                + hx / hy * kmidy(k2, x, y, y - hy);
                        d[j][i] = -hx / hy * kmidy(k2, x, y, y + hy);
                        e[j][i] = -hx / hy * kmidy(k2, x, y, y - hy);
                        if (i == 1) {
                            // для i = 1, j = 1 ... M - 2
                            a[j][i] = -hy / hx * kmidx(k1, x, x + hx, y);
                            b[j][i] = 0;
                            f[j][i] = hx * hy * source[j][i]
                                    + hy / hx * kmidx(k1, x, x - hx, y) * g1.apply(y);
                        } else if (i == N - 1) {
                            // для i = N - 1, j = 1 ... M - 2
                            a[j][i] = 0;
                            b[j][i] = -hy / hx * kmidx(k1, x, x - hx, y);
                            f[j][i] = hx * hy * source[j][i]
                                    + hy / hx * kmidx(k1, x, x + hx, y) * g2.apply(y);
                        } else {
                            // для i = 2 ... N - 2, j = 1 ... M - 2
                            a[j][i] = -hy / hx * kmidx(k1, x, x + hx, y);
                            b[j][i] = -hy / hx * kmidx(k1, x, x - hx, y);
                            f[j][i] = hx * hy * source[j][i];
                        }
                    }
                }
            }
        }
        return new DiffScheme(a, b, c, d, e, f);
    }

    double[][] countSourceMatrix() {
        double x = A; double y = C;
        double[][] f = new double[N][M];
        for (int j = 0; j < M; j++) {
            y = C + j * hy;
            for (int i = 0; i < N; i++) {
                x = A + i * hx;
                double dk1dx = Differentiation.derR(k1, x, y) * Differentiation.derR(u, x, y);
                double d2udx2 = k1.apply(x, y) * Differentiation.der2R(u, x, y);
                double dk2dy = Differentiation.derivativeByT(k2, x, y) * Differentiation.derivativeByT(u, x, y);
                double d2udy2 = k2.apply(x, y) * Differentiation.der2T(u, x, y);
                //System.out.printf("%f\t%f\t%f\t%f\n", dk1dx,d2udx2,dk2dy ,d2udy2);
                f[j][i] = - (dk1dx + d2udx2 + dk2dy + d2udy2);
            }
        }
        return f;
    }

    private double kmidx(BinaryOperator<Double> k, double x, double x_plus_minus, double y) {
        return (k.apply(x, y) + k.apply(x_plus_minus, y)) / 2;
    }

    private double kmidy(BinaryOperator<Double> k, double x, double y, double y_plus_minus) {
        return (k.apply(x, y) + k.apply(x, y_plus_minus)) / 2;
    }
}
