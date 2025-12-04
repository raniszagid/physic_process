public class TridiagonalMatrixCalculator {
    public static double[] multiplyTridiagonal(double[] a, double[] b, double[] c, double[] x) {
        int n = x.length;
        double[] result = new double[n];

        // Первая строка: b[0]*x[0] + c[0]*x[1]
        result[0] = b[0] * x[0] + c[0] * x[1];

        // Средние строки: a[i]*x[i-1] + b[i]*x[i] + c[i]*x[i+1]
        for (int i = 1; i < n - 1; i++) {
            result[i] = a[i] * x[i - 1] + b[i] * x[i] + c[i] * x[i + 1];
        }

        // Последняя строка: a[n-1]*x[n-2] + b[n-1]*x[n-1]
        result[n - 1] = a[n - 1] * x[n - 2] + b[n - 1] * x[n - 1];

        return result;
    }

    /**
     * Метод прогонки для решения трехдиагональной системы
     */
    public static double[] solveTridiagonal(double[] a, double[] b, double[] c, double[] f) {
        int n = f.length;
        double[] x = new double[n];

        // Прямой ход метода прогонки
        double[] alpha = new double[n];
        double[] beta = new double[n];

        // Первый шаг
        alpha[0] = -c[0] / b[0];
        beta[0] = f[0] / b[0];

        // Промежуточные шаги
        for (int i = 1; i < n - 1; i++) {
            double denominator = b[i] + a[i] * alpha[i - 1];
            alpha[i] = -c[i] / denominator;
            beta[i] = (f[i] - a[i] * beta[i - 1]) / denominator;
        }

        // Последний шаг прямого хода
        beta[n - 1] = (f[n - 1] - a[n - 1] * beta[n - 2]) /
                (b[n - 1] + a[n - 1] * alpha[n - 2]);

        // Обратный ход
        x[n - 1] = beta[n - 1];
        for (int i = n - 2; i >= 0; i--) {
            x[i] = alpha[i] * x[i + 1] + beta[i];
        }

        return x;
    }
}
