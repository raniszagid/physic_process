package lab__2;

import java.util.ArrayList;
import java.util.List;
import java.util.function.BiFunction;
import java.util.function.Function;

public class DifferentialEquationSolver {

    public static List<Double> linspace(double start, double end, int num) {
        List<Double> result = new ArrayList<>(num);
        double step = (end - start) / (num - 1);
        for (int i = 0; i < num; i++) {
            result.add(start + i * step);
        }
        return result;
    }

    public static double normInf(List<Double> vec) {
        double maxVal = 0.0;
        for (double value : vec) {
            maxVal = Math.max(maxVal, Math.abs(value));
        }
        return maxVal;
    }

    public static double normInf2d(List<List<Double>> mat) {
        double maxVal = 0.0;
        for (List<Double> row : mat) {
            for (double value : row) {
                maxVal = Math.max(maxVal, Math.abs(value));
            }
        }
        return maxVal;
    }

    public static List<Double> methodProgonki(List<Double> a, List<Double> b, List<Double> c, List<Double> g) {

        int N = g.size();
        List<Double> x = new ArrayList<>(N);
        for (int i = 0; i < N; i++) {
            x.add(0.0);
        }

        List<Double> alpha = new ArrayList<>(N + 1);
        List<Double> beta = new ArrayList<>(N + 1);
        for (int i = 0; i <= N; i++) {
            alpha.add(0.0);
            beta.add(0.0);
        }

        alpha.set(1, -b.get(0) / c.get(0));
        beta.set(1, g.get(0) / c.get(0));

        for (int i = 1; i < N - 1; i++) {
            double denominator = a.get(i) * alpha.get(i) + c.get(i);
            alpha.set(i + 1, -b.get(i) / denominator);
            beta.set(i + 1, (g.get(i) - a.get(i) * beta.get(i)) / denominator);
        }

        x.set(N - 1, (g.get(N - 1) - a.get(N - 1) * beta.get(N - 1)) /
                (a.get(N - 1) * alpha.get(N - 1) + c.get(N - 1)));

        for (int i = N - 2; i >= 0; i--) {
            x.set(i, alpha.get(i + 1) * x.get(i + 1) + beta.get(i + 1));
        }

        return x;
    }

    public static void calculateCoefficients(
            BiFunction<Double, Double, Double> k,
            BiFunction<Double, Double, Double> q,
            BiFunction<Double, Double, Double> f,
            double L,
            Function<Double, Double> nu1,
            double chi2,
            Function<Double, Double> nu2,
            int N,
            double tj,
            List<Double> a,
            List<Double> b,
            List<Double> c,
            List<Double> g) {

        double h = L / N;
        List<Double> x = linspace(0.0, L, N + 1);

        List<Double> xHalf = new ArrayList<>(N);
        for (int i = 0; i < N; i++) {
            xHalf.add(x.get(i) + h / 2);
        }

        a.clear();
        b.clear();
        c.clear();
        g.clear();

        for (int i = 0; i <= N; i++) {
            a.add(0.0);
            b.add(0.0);
            c.add(0.0);
            g.add(0.0);
        }

        for (int i = 1; i < N; i++) {
            a.set(i, k.apply(xHalf.get(i - 1), tj) / (h * h));
            b.set(i, k.apply(xHalf.get(i), tj) / (h * h));
            c.set(i, -a.get(i) - b.get(i) - q.apply(x.get(i), tj));
            g.set(i, f.apply(x.get(i), tj));
        }

        b.set(0, 2.0 * k.apply(xHalf.get(0), tj) / (h * h));
        c.set(0, -b.get(0) - q.apply(x.get(0), tj));
        g.set(0, f.apply(x.get(0), tj) + 2.0 * nu1.apply(tj) / h);

        a.set(N, 2.0 * k.apply(xHalf.get(N - 1), tj) / (h * h));
        c.set(N, -a.get(N) - q.apply(x.get(N), tj) - 2.0 * chi2 / h);
        g.set(N, f.apply(x.get(N), tj) + 2.0 * nu2.apply(tj) / h);
    }

    public static List<Double> tridiagMatvec(List<Double> mainDiag, List<Double> upperDiag, List<Double> lowerDiag,
                                             List<Double> vector) {

        int n = vector.size();
        List<Double> result = new ArrayList<>(n);
        for (int i = 0; i < n; i++) {
            result.add(0.0);
        }

        for (int i = 0; i < n; i++) {
            double sum = mainDiag.get(i) * vector.get(i);
            if (i > 0) {
                sum += lowerDiag.get(i) * vector.get(i - 1);
            }
            if (i < n - 1) {
                sum += upperDiag.get(i) * vector.get(i + 1);
            }
            result.set(i, sum);
        }

        return result;
    }

    public static List<List<Double>> solveDifferentialEquation(
            BiFunction<Double, Double, Double> k,
            BiFunction<Double, Double, Double> q,
            BiFunction<Double, Double, Double> f,
            double L,
            Function<Double, Double> nu1,
            double chi2,
            Function<Double, Double> nu2,
            int N,
            double T,
            Function<Double, Double> phi,
            int M,
            boolean explicitMethod) {

        List<List<Double>> u = new ArrayList<>(M + 1);
        for (int i = 0; i <= M; i++) {
            List<Double> row = new ArrayList<>(N + 1);
            for (int j = 0; j <= N; j++) {
                row.add(0.0);
            }
            u.add(row);
        }

        List<Double> xPoints = linspace(0.0, L, N + 1);
        for (int i = 0; i <= N; i++) {
            u.get(0).set(i, phi.apply(xPoints.get(i)));
        }

        List<Double> t = linspace(0.0, T, M + 1);
        double tau = T / M;

        for (int j = 0; j < M; j++) {
            if (explicitMethod) {
                List<Double> aj = new ArrayList<>();
                List<Double> bj = new ArrayList<>();
                List<Double> cj = new ArrayList<>();
                List<Double> gj = new ArrayList<>();

                calculateCoefficients(k, q, f, L, nu1, chi2, nu2, N, t.get(j),
                        aj, bj, cj, gj);

                List<Double> mainDiag = new ArrayList<>(N + 1);
                List<Double> upperDiag = new ArrayList<>(N + 1);
                List<Double> lowerDiag = new ArrayList<>(N + 1);

                for (int i = 0; i <= N; i++) {
                    mainDiag.add(tau * cj.get(i) + 1.0);
                    upperDiag.add(i < N ? tau * bj.get(i) : 0.0);
                    lowerDiag.add(i > 0 ? tau * aj.get(i) : 0.0);
                }

                List<Double> nextU = tridiagMatvec(mainDiag, upperDiag, lowerDiag, u.get(j));
                for (int i = 0; i <= N; i++) {
                    nextU.set(i, nextU.get(i) + tau * gj.get(i));
                }
                u.set(j + 1, nextU);
            } else {
                List<Double> ajp1 = new ArrayList<>();
                List<Double> bjp1 = new ArrayList<>();
                List<Double> cjp1 = new ArrayList<>();
                List<Double> gjp1 = new ArrayList<>();

                calculateCoefficients(k, q, f, L, nu1, chi2, nu2, N, t.get(j + 1),
                        ajp1, bjp1, cjp1, gjp1);

                List<Double> rightSide = new ArrayList<>(N + 1);
                for (int i = 0; i <= N; i++) {
                    rightSide.add(u.get(j).get(i) + tau * gjp1.get(i));
                }

                List<Double> aMod = new ArrayList<>(N + 1);
                List<Double> bMod = new ArrayList<>(N + 1);
                List<Double> cMod = new ArrayList<>(N + 1);

                for (int i = 0; i <= N; i++) {
                    aMod.add(-tau * ajp1.get(i));
                    bMod.add(-tau * bjp1.get(i));
                    cMod.add(-tau * cjp1.get(i) + 1.0);
                }

                u.set(j + 1, methodProgonki(aMod, bMod, cMod, rightSide));
            }
        }

        return u;
    }

    public static void testModel(
            BiFunction<Double, Double, Double> uExact,
            BiFunction<Double, Double, Double> k,
            BiFunction<Double, Double, Double> q,
            BiFunction<Double, Double, Double> f,
            double L,
            Function<Double, Double> nu1,
            double chi2,
            Function<Double, Double> nu2,
            List<Integer> NArray,
            double T,
            Function<Double, Double> phi,
            List<Integer> MArray,
            boolean explicitMethod) {

        List<List<Double>> epsArray = new ArrayList<>(NArray.size());
        for (int i = 0; i < NArray.size(); i++) {
            List<Double> row = new ArrayList<>(MArray.size());
            for (int j = 0; j < MArray.size(); j++) {
                row.add(0.0);
            }
            epsArray.add(row);
        }

        for (int i = 0; i < NArray.size(); i++) {
            int N = NArray.get(i);
            for (int j = 0; j < MArray.size(); j++) {
                int M = MArray.get(j);

                List<List<Double>> uCalc = solveDifferentialEquation(
                        k, q, f, L, nu1, chi2, nu2, N, T, phi, M, explicitMethod);

                List<Double> xPoints = linspace(0.0, L, N + 1);
                List<Double> tPoints = linspace(0.0, T, M + 1);

                List<List<Double>> uOrig = new ArrayList<>(M + 1);
                for (int ti = 0; ti <= M; ti++) {
                    List<Double> row = new ArrayList<>(N + 1);
                    for (int xi = 0; xi <= N; xi++) {
                        row.add(uExact.apply(xPoints.get(xi), tPoints.get(ti)));
                    }
                    uOrig.add(row);
                }

                List<List<Double>> eps = new ArrayList<>(M + 1);
                for (int ti = 0; ti <= M; ti++) {
                    List<Double> row = new ArrayList<>(N + 1);
                    for (int xi = 0; xi <= N; xi++) {
                        row.add(uCalc.get(ti).get(xi) - uOrig.get(ti).get(xi));
                    }
                    eps.add(row);
                }

                double epsNorm = normInf2d(eps);
                epsArray.get(i).set(j, epsNorm);
            }
        }

        System.out.printf("%-8s", "M/N");
        for (Integer n : NArray) {
            System.out.printf("%-16d", n);
        }
        System.out.println();

        for (int j = 0; j < MArray.size(); j++) {
            System.out.printf("%-8d", MArray.get(j));
            for (int i = 0; i < NArray.size(); i++) {
                System.out.printf("%-16.2e", epsArray.get(i).get(j));
            }
            System.out.println();
        }
    }

    public static void tests() {
        List<Integer> NArray = new ArrayList<>();
        List<Integer> MArray = new ArrayList<>();

        for (int i = 2; i <= 7; i++) {
            NArray.add(1 << i);
        }
        for (int i = 2; i <= 15; i++) {
            MArray.add(1 << i);
        }

        System.out.println("\nTest 1 (explicit Euler):");
        testModel(
                (x, t) -> 1.0,
                (x, t) -> 1.0,
                (x, t) -> 2.0,
                (x, t) -> 2.0,
                1.0,
                t -> 0.0,
                2.0,
                t -> 2.0,
                NArray,
                1.0,
                x -> 1.0,
                MArray,
                true
        );

        System.out.println("\nTest 1.1 (explicit Euler):");
        testModel(
                (x, t) -> 1.0,
                (x, t) -> 1.0,
                (x, t) -> 2.0,
                (x, t) -> 2.0,
                1.0,
                t -> 0.0,
                2.0,
                t -> 2.0,
                NArray,
                1.0,
                x -> 1.0 + 100 * Double.MIN_VALUE,
                MArray,
                true
        );

        System.out.println("\nTest 1 (implicit Euler):");
        testModel(
                (x, t) -> 1.0,
                (x, t) -> 1.0,
                (x, t) -> 2.0,
                (x, t) -> 2.0,
                1.0,
                t -> 0.0,
                2.0,
                t -> 2.0,
                NArray,
                1.0,
                x -> 1.0,
                MArray,
                false
        );

        System.out.println("\nTest 2 (explicit Euler):");
        testModel(
                (x, t) -> Math.sin(t),
                (x, t) -> 1.0 + x,
                (x, t) -> x,
                (x, t) -> x * Math.sin(t) + Math.cos(t),
                1.0,
                t -> 0.0,
                1.0,
                t -> Math.sin(t),
                NArray,
                1.0,
                x -> 0.0,
                MArray,
                true
        );

        System.out.println("\nTest 2 (implicit Euler):");
        testModel(
                (x, t) -> Math.sin(t),
                (x, t) -> 1.0 + x,
                (x, t) -> x,
                (x, t) -> x * Math.sin(t) + Math.cos(t),
                1.0,
                t -> 0.0,
                1.0,
                t -> Math.sin(t),
                NArray,
                1.0,
                x -> 0.0,
                MArray,
                false
        );

        System.out.println("\nTest 3 (explicit Euler):");
        testModel(
                (x, t) -> x * x * x * Math.sin(t) * Math.sin(t),
                (x, t) -> Math.exp(x),
                (x, t) -> 1.0,
                (x, t) -> x * x * x * Math.sin(t) * Math.sin(t) +
                        2 * x * x * x * Math.sin(t) * Math.cos(t) -
                        3 * x * x * Math.exp(x) * Math.sin(t) * Math.sin(t) -
                        6 * x * Math.exp(x) * Math.sin(t) * Math.sin(t),
                1.0,
                t -> 0.0,
                2.0,
                t -> (2 + 3 * Math.exp(1)) * Math.sin(t) * Math.sin(t),
                NArray,
                1.0,
                x -> 0.0,
                MArray,
                true
        );

        System.out.println("\nTest 3 (implicit Euler):");
        testModel(
                (x, t) -> x * x * x * Math.sin(t) * Math.sin(t),
                (x, t) -> Math.exp(x),
                (x, t) -> 1.0,
                (x, t) -> x * x * x * Math.sin(t) * Math.sin(t) +
                        2 * x * x * x * Math.sin(t) * Math.cos(t) -
                        3 * x * x * Math.exp(x) * Math.sin(t) * Math.sin(t) -
                        6 * x * Math.exp(x) * Math.sin(t) * Math.sin(t),
                1.0,
                t -> 0.0,
                2.0,
                t -> (2 + 3 * Math.exp(1)) * Math.sin(t) * Math.sin(t),
                NArray,
                1.0,
                x -> 0.0,
                MArray,
                false
        );
    }

    public static void main(String[] args) {
        tests();
    }
}