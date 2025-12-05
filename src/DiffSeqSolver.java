import java.util.function.UnaryOperator;

public class DiffSeqSolver {
    public static double[] countSourceMatrix(int N,
                                             double A,
                                             double B,
                                             UnaryOperator<Double> k,
                                             UnaryOperator<Double> u,
                                             UnaryOperator<Double> q) {
        double[] f = new double[N+1];
        double r = A;
        double h = (B - A) / N;
        for (int i = 0; i <= N; i++) {
            double dk = Differentiation.derivative(k, r, 1e-10);
            double du = Differentiation.derivative(u, r, 1e-10);
            double d2u = Differentiation.derivative2(u, r, 1e-10);
            double firstTerm = dk * du + k.apply(r) * du / r + k.apply(r) * d2u;
            f[i] = - firstTerm + q.apply(r) * u.apply(r);
            r += h;
        }
        return f;
    }

    public static MatrixSystem differenceScheme(double A, double B, int N,
                                                double nu1, double nu2, double xi1, double xi2,
                                                UnaryOperator<Double> k,
                                                UnaryOperator<Double> u,
                                                UnaryOperator<Double> q) {
        MatrixSystem sys = new MatrixSystem(countSourceMatrix(N, A, B, k, u, q));
        double h = (B - A) / N;
        double hh = h / 2.0;
        double r = A;

        // для i = 0
        sys.a[0] = 0;
        sys.b[0] = (r + hh) * k.apply(r + hh) / h
                + r * q.apply(r) * hh
                + xi1 * r;
        sys.c[0] = - (r + hh) * k.apply(r + hh) / h;
        sys.f[0] = r * (hh * sys.f[0] + nu1);

        // для i от 1 до N-1
        for (int i = 1; i <= N - 1; i++)
        {
            r += h;
            sys.a[i] = -((r - hh) * k.apply(r - hh)) / h;
            sys.b[i] = ((r + hh) * k.apply(r + hh)) / h
                    + ((r - hh) * k.apply(r - hh)) / h
                    + r * q.apply(r) * h;
            sys.c[i] = -((r + hh) * k.apply(r + hh)) / h;
            sys.f[i] = r * h * sys.f[i];
        }

        // для i = N
        r += h;
        sys.a[N] = -((r - hh) * k.apply(r - hh)) / h;
        sys.b[N] = (r - hh) * k.apply(r - hh) / h
                + r * xi2
                + r * q.apply(r) * hh;
        sys.c[N] = 0.0;
        sys.f[N] = r * hh * sys.f[N]
                + r * nu2;
        return sys;
    }
}
