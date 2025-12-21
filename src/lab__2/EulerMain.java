package lab__2;

import java.util.List;
import java.util.function.BinaryOperator;

public class EulerMain {
    static int testNumber;
    public static void main(String[] args) {
        test(
                (r, t) -> 1d, // u
                (r, t) -> 1d, // k
                (r, t) -> 1d); // q

        test(
                (r, t) -> r*t, // u
                (r, t) -> 1d, // k
                (r, t) -> 1d); // q

        test(
                (r, t) -> r*t, // u
                (r, t) -> r*t, // k
                (r, t) -> 1d); // q
    }

    public static void test(BinaryOperator<Double> u, BinaryOperator<Double> k, BinaryOperator<Double> q) {
        List<Integer> n = List.of(4, 8, 16, 32, 64, 128, 256);
        List<Integer> m = List.of(10, 20, 50, 100);
        double[][] explicitNev = new double[m.size()][n.size()];
        double[][] explicitError = new double[m.size()][n.size()];
        double[][] implicitNev = new double[m.size()][n.size()];
        double[][] implicitError = new double[m.size()][n.size()];
        for (int j = 0; j < m.size(); j++) {
            for (int i = 0; i < n.size(); i++) {
                EulerSolver eulerSolver = new EulerSolver(n.get(i), m.get(j),
                        1, 9,
                        50, u, k, q, 2);
                eulerSolver.explicit();
                explicitError[j][i] = eulerSolver.maxError;
                explicitNev[j][i] = eulerSolver.maxNev;
                eulerSolver.implicit();
                implicitError[j][i] = eulerSolver.maxError;
                implicitNev[j][i] = eulerSolver.maxNev;
            }
        }
        System.out.println("Тест " + (++testNumber));
        System.out.println("Ошибка явного метода");
        System.out.print("  \t");
        n.forEach(x -> System.out.printf("%d\t", x));
        System.out.println();
        for (int j = 0; j < m.size(); j++) {
            System.out.printf("%d\t", m.get(j));
            for (int i = 0; i < n.size(); i++) {
                System.out.printf("%.3e\t", explicitError[j][i]);
            }
            System.out.println();
        }
        System.out.println();

        System.out.println("Невязка явного метода");
        System.out.print("  \t");
        n.forEach(x -> System.out.printf("%d\t", x));
        System.out.println();
        for (int j = 0; j < m.size(); j++) {
            System.out.printf("%d\t", m.get(j));
            for (int i = 0; i < n.size(); i++) {
                System.out.printf("%.3e\t", explicitNev[j][i]);
            }
            System.out.println();
        }
        System.out.println();


        System.out.println("Ошибка неявного метода");
        System.out.print("  \t");
        n.forEach(x -> System.out.printf("%d\t", x));
        System.out.println();
        for (int j = 0; j < m.size(); j++) {
            System.out.printf("%d\t", m.get(j));
            for (int i = 0; i < n.size(); i++) {
                System.out.printf("%.3e\t", implicitError[j][i]);
            }
            System.out.println();
        }
        System.out.println();

        System.out.println("Невязка неявного метода");
        System.out.print("  \t");
        n.forEach(x -> System.out.printf("%d\t", x));
        System.out.println();
        for (int j = 0; j < m.size(); j++) {
            System.out.printf("%d\t", m.get(j));
            for (int i = 0; i < n.size(); i++) {
                System.out.printf("%.3e\t", implicitNev[j][i]);
            }
            System.out.println();
        }
        System.out.println();
    }
}
