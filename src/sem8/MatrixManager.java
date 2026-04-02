package sem8;

import java.util.Arrays;

public class MatrixManager {
    public static class TridiagonalSystem {
        public double[][][] a;  // поддиагональные блоки (нижняя диагональ)
        public double[][][] b;  // диагональные блоки
        public double[][][] c;  // наддиагональные блоки (верхняя диагональ)
        public double[][] f;    // правая часть
        public int sizeX;       // количество узлов по x
        public int sizeY;       // количество узлов по y

        public TridiagonalSystem(int sizeY, int sizeX) {
            this.sizeY = sizeY;
            this.sizeX = sizeX;
            this.a = new double[sizeY][sizeX][sizeX];
            this.b = new double[sizeY][sizeX][sizeX];
            this.c = new double[sizeY][sizeX][sizeX];
            this.f = new double[sizeY][sizeX];
        }
    }

    /**
     * Преобразование 5-диагональной системы в блочно-трехдиагональную
     */
    public TridiagonalSystem convertToTridiagonal(double[][] a, double[][] b, double[][] c,
                                                  double[][] d, double[][] e, double[][] f,
                                                  int N, int M) {
        TridiagonalSystem system = new TridiagonalSystem(M, N);

        for (int j = 0; j < M; j++) {
            for (int i = 0; i < N; i++) {
                // Диагональный блок B_j
                system.b[j][i][i] = c[j][i];

                // Связи внутри слоя (по x)
                if (i > 0) {
                    system.b[j][i][i - 1] = b[j][i];
                }
                if (i < N - 1) {
                    system.b[j][i][i + 1] = a[j][i];
                }

                // Связь с предыдущим слоем по y (A_j)
                if (j > 0) {
                    system.a[j][i][i] = e[j][i];
                }

                // Связь со следующим слоем по y (C_j)
                if (j < M - 1) {
                    system.c[j][i][i] = d[j][i];
                }

                // Правая часть
                system.f[j][i] = f[j][i];
            }
        }

        return system;
    }

    /**
     * Матричная прогонка для решения блочно-трехдиагональной системы
     */
    public double[][] matrixProgonka(double[][][] a, double[][][] b, double[][][] c,
                                     double[][] f, int sizeY, int sizeX) {

        System.out.println("=== Начало матричной прогонки ===");
        System.out.printf("Размеры: sizeY=%d, sizeX=%d\n", sizeY, sizeX);

        // Прогоночные коэффициенты
        double[][][] alpha = new double[sizeY][sizeX][sizeX];
        double[][] beta = new double[sizeY][sizeX];

        // ===== Прямой ход =====

        // j = 0: B0 * U0 + C0 * U1 = F0
        System.out.println("\nШаг 1: Решение для j=0");
        double[][] invB0 = solveMatrixSystem(b[0], null);
        if (invB0 == null) {
            throw new RuntimeException("Не удалось найти обратную матрицу для B0");
        }

        // alpha[0] = -B0^(-1) * C0
        alpha[0] = multiplyMatrices(invB0, c[0]);
        for (int i = 0; i < sizeX; i++) {
            for (int j = 0; j < sizeX; j++) {
                alpha[0][i][j] = -alpha[0][i][j];
            }
        }

        // beta[0] = B0^(-1) * F0
        beta[0] = multiplyMatrixVector(invB0, f[0]);

        // Прямой ход для j = 1, 2, ..., sizeY-2
        for (int j = 1; j < sizeY - 1; j++) {
            System.out.printf("\nШаг %d: Прямой ход для j=%d\n", j+1, j);

            // B' = B_j + A_j * alpha[j-1]
            double[][] B_prime = new double[sizeX][sizeX];
            double[][] A_alpha = multiplyMatrices(a[j], alpha[j - 1]);

            for (int i = 0; i < sizeX; i++) {
                for (int k = 0; k < sizeX; k++) {
                    B_prime[i][k] = b[j][i][k] + A_alpha[i][k];
                }
            }

            // F' = F_j - A_j * beta[j-1]
            double[] F_prime = new double[sizeX];
            double[] A_beta = multiplyMatrixVector(a[j], beta[j - 1]);
            for (int i = 0; i < sizeX; i++) {
                F_prime[i] = f[j][i] - A_beta[i];
            }

            // Находим обратную матрицу для B'
            double[][] invB_prime = solveMatrixSystem(B_prime, null);
            if (invB_prime == null) {
                throw new RuntimeException(String.format("Не удалось найти обратную матрицу для B' на шаге j=%d", j));
            }

            // alpha[j] = -inv(B') * C_j
            alpha[j] = multiplyMatrices(invB_prime, c[j]);
            for (int i = 0; i < sizeX; i++) {
                for (int k = 0; k < sizeX; k++) {
                    alpha[j][i][k] = -alpha[j][i][k];
                }
            }

            // beta[j] = inv(B') * F'
            beta[j] = multiplyMatrixVector(invB_prime, F_prime);
        }

        // ===== Обратный ход =====

        // Последний слой j = sizeY-1
        System.out.println("\nШаг: Обратный ход для последнего слоя");
        double[][] B_last_prime = new double[sizeX][sizeX];
        double[][] A_last_alpha = multiplyMatrices(a[sizeY - 1], alpha[sizeY - 2]);

        for (int i = 0; i < sizeX; i++) {
            for (int k = 0; k < sizeX; k++) {
                B_last_prime[i][k] = b[sizeY - 1][i][k] + A_last_alpha[i][k];
            }
        }

        double[] F_last_prime = new double[sizeX];
        double[] A_last_beta = multiplyMatrixVector(a[sizeY - 1], beta[sizeY - 2]);
        for (int i = 0; i < sizeX; i++) {
            F_last_prime[i] = f[sizeY - 1][i] - A_last_beta[i];
        }

        // Решаем B_last_prime * U_last = F_last_prime
        double[][] result = new double[sizeY][sizeX];
        result[sizeY - 1] = solveLinearSystem(B_last_prime, F_last_prime);
        if (result[sizeY - 1] == null) {
            throw new RuntimeException("Не удалось решить систему для последнего слоя");
        }

        // Обратный ход для j = sizeY-2, ..., 0
        System.out.println("\nШаг: Обратный ход для остальных слоев");
        for (int j = sizeY - 2; j >= 0; j--) {
            // U_j = alpha_j * U_{j+1} + beta_j
            result[j] = addVectors(multiplyMatrixVector(alpha[j], result[j + 1]), beta[j]);
        }

        System.out.println("\n=== Матричная прогонка завершена ===");
        return result;
    }

    /**
     * Решение матричной системы или нахождение обратной матрицы
     * @param A матрица
     * @param b правая часть (если null, то находится обратная матрица)
     * @return если b != null, возвращается матрица 1xN с решением,
     *         если b == null, возвращается обратная матрица
     */
    private double[][] solveMatrixSystem(double[][] A, double[] b) {
        int n = A.length;

        // Проверка на вырожденность перед решением
        if (isSingular(A)) {
            System.err.println("Предупреждение: матрица вырождена, добавляем регуляризацию");
            A = regularizeMatrix(A);
        }

        if (b == null) {
            // Находим обратную матрицу
            double[][] inverse = new double[n][n];
            for (int i = 0; i < n; i++) {
                double[] e = new double[n];
                e[i] = 1.0;
                double[] column = solveLinearSystem(A, e);
                if (column == null) {
                    System.err.printf("Не удалось найти %d-й столбец обратной матрицы\n", i);
                    return null;
                }
                for (int j = 0; j < n; j++) {
                    inverse[j][i] = column[j];
                }
            }
            return inverse;
        } else {
            // Решаем систему A * x = b
            double[] solution = solveLinearSystem(A, b);
            if (solution == null) return null;
            double[][] result = new double[1][n];
            result[0] = solution;
            return result;
        }
    }

    /**
     * Решение линейной системы методом Гаусса с выбором главного элемента
     */
    private double[] solveLinearSystem(double[][] A, double[] b) {
        int n = A.length;
        double[][] augmented = new double[n][n + 1];

        // Копируем матрицу и правую часть
        for (int i = 0; i < n; i++) {
            System.arraycopy(A[i], 0, augmented[i], 0, n);
            augmented[i][n] = b[i];
        }

        // Прямой ход метода Гаусса
        for (int i = 0; i < n; i++) {
            // Поиск главного элемента
            int maxRow = i;
            double maxVal = Math.abs(augmented[i][i]);
            for (int k = i + 1; k < n; k++) {
                if (Math.abs(augmented[k][i]) > maxVal) {
                    maxVal = Math.abs(augmented[k][i]);
                    maxRow = k;
                }
            }

            // Проверка на вырожденность
            if (maxVal < 1e-12) {
                System.err.printf("Ошибка: нулевой столбец %d, maxVal = %e\n", i, maxVal);
                printMatrix(A, b);
                return null;
            }

            // Перестановка строк
            if (maxRow != i) {
                double[] temp = augmented[i];
                augmented[i] = augmented[maxRow];
                augmented[maxRow] = temp;
            }

            // Нормировка строки
            double pivot = augmented[i][i];
            for (int j = i; j <= n; j++) {
                augmented[i][j] /= pivot;
            }

            // Исключение в остальных строках
            for (int k = 0; k < n; k++) {
                if (k != i && Math.abs(augmented[k][i]) > 1e-12) {
                    double factor = augmented[k][i];
                    for (int j = i; j <= n; j++) {
                        augmented[k][j] -= factor * augmented[i][j];
                    }
                }
            }
        }

        // Формируем решение
        double[] x = new double[n];
        for (int i = 0; i < n; i++) {
            x[i] = augmented[i][n];
        }

        return x;
    }

    /**
     * Проверка матрицы на вырожденность
     */
    private boolean isSingular(double[][] A) {
        int n = A.length;
        double det = 1.0;
        double[][] temp = new double[n][n];

        for (int i = 0; i < n; i++) {
            System.arraycopy(A[i], 0, temp[i], 0, n);
        }

        // Приводим к верхнетреугольному виду
        for (int i = 0; i < n; i++) {
            int maxRow = i;
            double maxVal = Math.abs(temp[i][i]);
            for (int k = i + 1; k < n; k++) {
                if (Math.abs(temp[k][i]) > maxVal) {
                    maxVal = Math.abs(temp[k][i]);
                    maxRow = k;
                }
            }

            if (maxVal < 1e-12) {
                return true; // Матрица вырождена
            }

            if (maxRow != i) {
                double[] row = temp[i];
                temp[i] = temp[maxRow];
                temp[maxRow] = row;
                det *= -1;
            }

            det *= temp[i][i];

            for (int k = i + 1; k < n; k++) {
                double factor = temp[k][i] / temp[i][i];
                for (int j = i; j < n; j++) {
                    temp[k][j] -= factor * temp[i][j];
                }
            }
        }

        return Math.abs(det) < 1e-12;
    }

    /**
     * Регуляризация матрицы (добавление малой диагонали)
     */
    private double[][] regularizeMatrix(double[][] A) {
        int n = A.length;
        double[][] result = new double[n][n];
        double lambda = 1e-8;

        for (int i = 0; i < n; i++) {
            System.arraycopy(A[i], 0, result[i], 0, n);
            result[i][i] += lambda;
        }

        System.out.printf("Добавлена регуляризация с lambda = %e\n", lambda);
        return result;
    }

    /**
     * Умножение матрицы на вектор
     */
    private double[] multiplyMatrixVector(double[][] matrix, double[] vector) {
        int n = matrix.length;
        double[] result = new double[n];
        for (int i = 0; i < n; i++) {
            result[i] = 0;
            for (int j = 0; j < n; j++) {
                result[i] += matrix[i][j] * vector[j];
            }
        }
        return result;
    }

    /**
     * Умножение двух матриц
     */
    private double[][] multiplyMatrices(double[][] A, double[][] B) {
        int n = A.length;
        int m = B[0].length;
        int p = A[0].length;
        double[][] result = new double[n][m];

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                result[i][j] = 0;
                for (int k = 0; k < p; k++) {
                    result[i][j] += A[i][k] * B[k][j];
                }
            }
        }
        return result;
    }

    /**
     * Сложение двух векторов
     */
    private double[] addVectors(double[] v1, double[] v2) {
        int n = v1.length;
        double[] result = new double[n];
        for (int i = 0; i < n; i++) {
            result[i] = v1[i] + v2[i];
        }
        return result;
    }

    /**
     * Вывод матрицы для отладки
     */
    private void printMatrix(double[][] A, double[] b) {
        System.out.println("Матрица системы:");
        for (int i = 0; i < A.length; i++) {
            for (int j = 0; j < A[0].length; j++) {
                System.out.printf("%12.6f ", A[i][j]);
            }
            if (b != null) {
                System.out.printf("| %12.6f\n", b[i]);
            } else {
                System.out.println();
            }
        }
    }

    /**
     * Проверка диагонального преобладания (для отладки)
     */
    public void checkDiagonalDominance(double[][] A) {
        int n = A.length;
        boolean hasDominance = true;

        for (int i = 0; i < n; i++) {
            double diag = Math.abs(A[i][i]);
            double sum = 0.0;
            for (int j = 0; j < n; j++) {
                if (j != i) sum += Math.abs(A[i][j]);
            }

            if (diag < sum) {
                System.out.printf("Строка %d: |%f| < %f - нет диагонального преобладания\n", i, diag, sum);
                hasDominance = false;
            } else {
                System.out.printf("Строка %d: |%f| >= %f - OK\n", i, diag, sum);
            }
        }

        if (hasDominance) {
            System.out.println("Матрица имеет диагональное преобладание");
        } else {
            System.out.println("Матрица НЕ имеет диагонального преобладания");
        }
    }
    /*public MatrixSystem2D convertToTridiagonal(double[][] a, double[][] b, double[][] c, double[][] d,
                                               double[][] e, double[][] f, int sizeX, int sizeY) {

        // Для матричной прогонки создаем блочную систему:
        // A_j * U_{j-1} + B_j * U_j + C_j * U_{j+1} = F_j

        // Каждый U_j — это вектор [u_{j,0}, u_{j,1}, ..., u_{j,N-1}]
        // A_j, B_j, C_j — матрицы размера N×N

        double[][][] aa = new double[sizeY][sizeX][sizeX]; // Нижняя диагональ (связь с j-1)
        double[][][] bb = new double[sizeY][sizeX][sizeX]; // Центральная диагональ
        double[][][] cc = new double[sizeY][sizeX][sizeX]; // Верхняя диагональ (связь с j+1)
        double[][] ff = new double[sizeY][sizeX];          // Правая часть

        // Заполняем блоки
        for (int j = 0; j < sizeY; j++) {
            // Обнуляем матрицы
            for (int i = 0; i < sizeX; i++) {
                for (int k = 0; k < sizeX; k++) {
                    aa[j][i][k] = 0.0;
                    bb[j][i][k] = 0.0;
                    cc[j][i][k] = 0.0;
                }
            }

            // Заполняем B_j (центральная матрица для строки j)
            for (int i = 0; i < sizeX; i++) {
                // Диагональный элемент B[j][i][i]
                bb[j][i][i] = c[j][i];

                // Внедиагональные элементы (связи по x)
                if (i > 0) {
                    bb[j][i][i - 1] = b[j][i]; // связь с i-1
                }
                if (i < sizeX - 1) {
                    bb[j][i][i + 1] = a[j][i]; // связь с i+1
                }
            }

            // Заполняем A_j (связь с нижней строкой j-1)
            if (j > 0) {
                for (int i = 0; i < sizeX; i++) {
                    aa[j][i][i] = e[j][i]; // e[j][i] — это коэффициент при u_{j-1,i}
                }
            }

            // Заполняем C_j (связь с верхней строкой j+1)
            if (j < sizeY - 1) {
                for (int i = 0; i < sizeX; i++) {
                    cc[j][i][i] = d[j][i]; // d[j][i] — это коэффициент при u_{j+1,i}
                }
            }

            // Заполняем правую часть F_j
            for (int i = 0; i < sizeX; i++) {
                ff[j][i] = f[j][i];
            }
        }
        return new MatrixSystem2D(aa, bb, cc, ff, sizeX, sizeY);
    }
    public double[][] matrixProgonka(double[][][] A, double[][][] B,
                                     double[][][] C, double[][] F,
                                     int sizeY, int sizeX) {
        // Прямой ход
        double[][][] alpha = new double[sizeY][sizeX][sizeX];
        double[][] beta = new double[sizeY][sizeX];

        // Инициализация для j = 0
        double[][] B0_inv = invertMatrix(B[0]);
        alpha[0] = multiplyMatrices(B0_inv, C[0]);
        beta[0] = multiplyMatrixVector(B0_inv, F[0]); // F[0] - это double[], возвращает double[]

        // Прямой ход для j = 1,...,sizeY-1
        for (int j = 1; j < sizeY; j++) {
            // Вычисляем B_j - A_j * alpha_{j-1}
            double[][] temp = subtractMatrices(B[j], multiplyMatrices(A[j], alpha[j-1]));
            double[][] temp_inv = invertMatrix(temp);

            // alpha_j = temp_inv * C_j
            alpha[j] = multiplyMatrices(temp_inv, C[j]);

            // beta_j = temp_inv * (F_j - A_j * beta_{j-1})
            double[] tempVec = subtractVectors(F[j], multiplyMatrixVector(A[j], beta[j-1]));
            beta[j] = multiplyMatrixVector(temp_inv, tempVec); // temp_inv - double[][], tempVec - double[]
        }

        // Обратный ход
        double[][] U = new double[sizeY][sizeX];

        // Последняя строка
        System.arraycopy(beta[sizeY-1], 0, U[sizeY-1], 0, sizeX);

        // Обратный ход для j = sizeY-2,...,0
        for (int j = sizeY-2; j >= 0; j--) {
            // U_j = beta_j - alpha_j * U_{j+1}
            double[] alphaU = multiplyMatrixVector(alpha[j], U[j+1]);
            for (int i = 0; i < sizeX; i++) {
                U[j][i] = beta[j][i] - alphaU[i];
            }
        }

        return U;
    }

    // Умножение матрицы на вектор (возвращает вектор)
    private double[] multiplyMatrixVector(double[][] A, double[] v) {
        int n = A.length;
        int m = A[0].length;
        double[] result = new double[n];

        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < m; j++) {
                sum += A[i][j] * v[j];
            }
            result[i] = sum;
        }
        return result;
    }

    // Умножение матрицы на матрицу (возвращает матрицу)
    private double[][] multiplyMatrices(double[][] A, double[][] B) {
        int n = A.length;
        int m = B[0].length;
        int p = B.length;
        double[][] result = new double[n][m];

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                double sum = 0.0;
                for (int k = 0; k < p; k++) {
                    sum += A[i][k] * B[k][j];
                }
                result[i][j] = sum;
            }
        }
        return result;
    }

    // Вычитание векторов (возвращает вектор)
    private double[] subtractVectors(double[] v1, double[] v2) {
        int n = v1.length;
        double[] result = new double[n];

        for (int i = 0; i < n; i++) {
            result[i] = v1[i] - v2[i];
        }
        return result;
    }

    private double[][] invertMatrix(double[][] matrix) {
        int n = matrix.length;
        double[][] augmented = new double[n][2*n];

        // Создаем расширенную матрицу [A | I]
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                augmented[i][j] = matrix[i][j];
                augmented[i][j + n] = (i == j) ? 1.0 : 0.0;
            }
        }

        // Метод Гаусса для нахождения обратной матрицы
        for (int i = 0; i < n; i++) {
            // Поиск максимального элемента в столбце
            int maxRow = i;
            for (int k = i + 1; k < n; k++) {
                if (Math.abs(augmented[k][i]) > Math.abs(augmented[maxRow][i])) {
                    maxRow = k;
                }
            }

            // Перестановка строк
            double[] temp = augmented[i];
            augmented[i] = augmented[maxRow];
            augmented[maxRow] = temp;

            // Нормализация строки
            double pivot = augmented[i][i];
            for (int j = 0; j < 2*n; j++) {
                augmented[i][j] /= pivot;
            }

            // Исключение остальных строк
            for (int k = 0; k < n; k++) {
                if (k != i) {
                    double factor = augmented[k][i];
                    for (int j = 0; j < 2*n; j++) {
                        augmented[k][j] -= factor * augmented[i][j];
                    }
                }
            }
        }

        // Извлекаем обратную матрицу
        double[][] inverse = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                inverse[i][j] = augmented[i][j + n];
            }
        }

        return inverse;
    }

    private double[][] subtractMatrices(double[][] A, double[][] B) {
        int n = A.length;
        double[][] result = new double[n][n];

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                result[i][j] = A[i][j] - B[i][j];
            }
        }
        return result;
    }*/
}
