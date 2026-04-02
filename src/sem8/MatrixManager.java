package sem8;

import java.util.Arrays;

public class MatrixManager {

    /**
     * Структура для хранения трехдиагональной блочной системы
     */
    public static class TridiagonalSystem {
        public double[][][] a;  // нижняя диагональ (связь с j-1)
        public double[][][] b;  // главная диагональ
        public double[][][] c;  // верхняя диагональ (связь с j+1)
        public double[][] f;    // правая часть
        public int sizeY;       // количество блоков по y
        public int sizeX;       // размер каждого блока (по x)

        public TridiagonalSystem(double[][][] A, double[][][] B, double[][][] C,
                                 double[][] F, int sizeY, int sizeX) {
            this.a = A;
            this.b = B;
            this.c = C;
            this.f = F;
            this.sizeY = sizeY;
            this.sizeX = sizeX;
        }
    }

    /**
     * Преобразование пятидиагональной системы в трехдиагональную блочную
     * для решения методом матричной прогонки
     *
     * @param a коэффициент при u[i+1][j]
     * @param b коэффициент при u[i-1][j]
     * @param c коэффициент при u[i][j]
     * @param d коэффициент при u[i][j+1]
     * @param e коэффициент при u[i][j-1]
     * @param f правая часть
     * @param N количество узлов по x
     * @param M количество узлов по y
     * @return трехдиагональная блочная система
     */
    public TridiagonalSystem convertToTridiagonal(double[][] a, double[][] b, double[][] c,
                                                  double[][] d, double[][] e, double[][] f,
                                                  int N, int M) {
        int sizeY = M;  // количество строк (по y)
        int sizeX = N;  // количество столбцов (по x)

        System.out.println("=== convertToTridiagonal ===");
        System.out.printf("sizeY=%d, sizeX=%d\n", sizeY, sizeX);

        // Инициализация блочных матриц
        double[][][] A = new double[sizeY][sizeX][sizeX];  // нижняя диагональ
        double[][][] B = new double[sizeY][sizeX][sizeX];  // главная диагональ
        double[][][] C = new double[sizeY][sizeX][sizeX];  // верхняя диагональ
        double[][] F = new double[sizeY][sizeX];           // правая часть

        for (int j = 0; j < sizeY; j++) {
            for (int i = 0; i < sizeX; i++) {
                // Обнуляем все элементы
                for (int k = 0; k < sizeX; k++) {
                    if (j > 0) A[j][i][k] = 0.0;
                    B[j][i][k] = 0.0;
                    if (j < sizeY - 1) C[j][i][k] = 0.0;
                }

                // Заполняем главную диагональ B[j]
                // Коэффициент при u[i][j] (центр)
                B[j][i][i] = c[j][i];

                // Связь с левым соседом u[i-1][j]
                if (i > 0) {
                    B[j][i][i-1] = b[j][i];
                }

                // Связь с правым соседом u[i+1][j]
                if (i < sizeX - 1) {
                    B[j][i][i+1] = a[j][i];
                }

                // Связь с нижним слоем u[i][j-1] (матрица A)
                if (j > 0) {
                    A[j][i][i] = e[j][i];  // e - коэффициент при u[i][j-1]
                }

                // Связь с верхним слоем u[i][j+1] (матрица C)
                if (j < sizeY - 1) {
                    C[j][i][i] = d[j][i];  // d - коэффициент при u[i][j+1]
                }

                // Правая часть
                F[j][i] = f[j][i];
            }
        }

        // Диагностический вывод
        System.out.println("\n=== Проверка B[0] ===");
        for (int i = 0; i < Math.min(3, sizeX); i++) {
            System.out.printf("B[0][%d]: ", i);
            for (int k = 0; k < Math.min(3, sizeX); k++) {
                System.out.printf("%.2f ", B[0][i][k]);
            }
            System.out.println();
        }

        System.out.println("\n=== Проверка F[0] ===");
        for (int i = 0; i < Math.min(3, sizeX); i++) {
            System.out.printf("F[0][%d] = %.2f\n", i, F[0][i]);
        }

        return new TridiagonalSystem(A, B, C, F, sizeY, sizeX);
    }
    /**
     * Модифицированная матричная прогонка для решения блочной трехдиагональной системы
     *
     * Система имеет вид:
     * B[0]*U0 + C[0]*U1 = F0
     * A[j]*U_{j-1} + B[j]*Uj + C[j]*U_{j+1} = Fj, j=1..sizeY-2
     * A[sizeY-1]*U_{sizeY-2} + B[sizeY-1]*U_{sizeY-1} = F_{sizeY-1}
     *
     * @param A нижняя диагональ (матрицы sizeX x sizeX)
     * @param B главная диагональ (матрицы sizeX x sizeX)
     * @param C верхняя диагональ (матрицы sizeX x sizeX)
     * @param F правая часть (векторы sizeX)
     * @param sizeY количество слоев по y
     * @param sizeX размер каждого слоя (количество узлов по x)
     * @return решение U (матрица sizeY x sizeX)
     */
    public double[][] matrixProgonka(double[][][] A, double[][][] B, double[][][] C,
                                     double[][] F, int sizeY, int sizeX) {
        System.out.println("=== Начало матричной прогонки ===");
        System.out.printf("Размеры: sizeY=%d, sizeX=%d\n\n", sizeY, sizeX);

        // ДИАГНОСТИКА: выводим B[0] и F[0]
        System.out.println("=== B[0] (первая матрица) ===");
        for (int i = 0; i < sizeX; i++) {
            for (int j = 0; j < sizeX; j++) {
                System.out.printf("%8.3f ", B[0][i][j]);
            }
            System.out.println();
        }

        System.out.println("\n=== F[0] (первая правая часть) ===");
        for (int i = 0; i < sizeX; i++) {
            System.out.printf("F[0][%d] = %.3f\n", i, F[0][i]);
        }

        // Проверка на вырожденность B[0]
        double det = determinant(B[0]);
        System.out.printf("\nОпределитель B[0] = %.6f\n", det);
        if (Math.abs(det) < 1e-10) {
            System.err.println("ВНИМАНИЕ: B[0] вырождена!");
        }

        // Прогоночные коэффициенты: Uj = alpha[j] * U_{j+1} + beta[j]
        double[][][] alpha = new double[sizeY][sizeX][sizeX];
        double[][] beta = new double[sizeY][sizeX];

        // Прямой ход

        // Шаг 1: Вычисляем alpha[0] и beta[0] из первого уравнения
        // B[0]*U0 + C[0]*U1 = F0
        // U0 = -inv(B[0])*C[0]*U1 + inv(B[0])*F0
        // alpha[0] = -inv(B[0])*C[0]
        // beta[0] = inv(B[0])*F0

        System.out.println("Шаг 1: Решение для j=0");

        // Копируем B[0] для обращения
        double[][] B0 = new double[sizeX][sizeX];
        for (int i = 0; i < sizeX; i++) {
            System.arraycopy(B[0][i], 0, B0[i], 0, sizeX);
        }

        // Вычисляем inv(B0)
        double[][] invB0 = invertMatrix(B0);
        if (invB0 == null) {
            System.err.println("Ошибка: матрица B[0] вырождена");
            return null;
        }

        // alpha[0] = -inv(B0) * C[0]
        for (int i = 0; i < sizeX; i++) {
            for (int j = 0; j < sizeX; j++) {
                alpha[0][i][j] = 0;
                for (int k = 0; k < sizeX; k++) {
                    alpha[0][i][j] -= invB0[i][k] * C[0][k][j];
                }
            }
        }

        // beta[0] = inv(B0) * F[0]
        for (int i = 0; i < sizeX; i++) {
            beta[0][i] = 0;
            for (int k = 0; k < sizeX; k++) {
                beta[0][i] += invB0[i][k] * F[0][k];
            }
        }

        // Шаг 2: Прямой ход для j = 1..sizeY-2
        // A[j]*U_{j-1} + B[j]*Uj + C[j]*U_{j+1} = Fj
        // Подставляем U_{j-1} = alpha[j-1]*Uj + beta[j-1]
        // Получаем: (A[j]*alpha[j-1] + B[j])*Uj + C[j]*U_{j+1} = Fj - A[j]*beta[j-1]
        // Обозначим: B_tilda = A[j]*alpha[j-1] + B[j]
        // Тогда: Uj = -inv(B_tilda)*C[j]*U_{j+1} + inv(B_tilda)*(Fj - A[j]*beta[j-1])

        System.out.println("\nШаг 2: Прямой ход для j=1.." + (sizeY-2));

        for (int j = 1; j < sizeY - 1; j++) {
            // Вычисляем A_j * alpha_{j-1}
            double[][] A_alpha = new double[sizeX][sizeX];
            for (int i = 0; i < sizeX; i++) {
                for (int k = 0; k < sizeX; k++) {
                    A_alpha[i][k] = 0;
                    for (int l = 0; l < sizeX; l++) {
                        A_alpha[i][k] += A[j][i][l] * alpha[j-1][l][k];
                    }
                }
            }

            // B_tilda = A_alpha + B[j]
            double[][] B_tilda = new double[sizeX][sizeX];
            for (int i = 0; i < sizeX; i++) {
                for (int k = 0; k < sizeX; k++) {
                    B_tilda[i][k] = A_alpha[i][k] + B[j][i][k];
                }
            }

            // Вычисляем правую часть: F_tilda = F[j] - A[j]*beta[j-1]
            double[] F_tilda = new double[sizeX];
            for (int i = 0; i < sizeX; i++) {
                F_tilda[i] = F[j][i];
                for (int k = 0; k < sizeX; k++) {
                    F_tilda[i] -= A[j][i][k] * beta[j-1][k];
                }
            }

            // Вычисляем inv(B_tilda)
            double[][] invB_tilda = invertMatrix(B_tilda);
            if (invB_tilda == null) {
                System.err.printf("Ошибка: матрица B_tilda[%d] вырождена\n", j);
                return null;
            }

            // alpha[j] = -inv(B_tilda) * C[j]
            for (int i = 0; i < sizeX; i++) {
                for (int k = 0; k < sizeX; k++) {
                    alpha[j][i][k] = 0;
                    for (int l = 0; l < sizeX; l++) {
                        alpha[j][i][k] -= invB_tilda[i][l] * C[j][l][k];
                    }
                }
            }

            // beta[j] = inv(B_tilda) * F_tilda
            for (int i = 0; i < sizeX; i++) {
                beta[j][i] = 0;
                for (int k = 0; k < sizeX; k++) {
                    beta[j][i] += invB_tilda[i][k] * F_tilda[k];
                }
            }
        }

        // Шаг 3: Обратный ход для последнего слоя j = sizeY-1
        // A[sizeY-1]*U_{sizeY-2} + B[sizeY-1]*U_{sizeY-1} = F[sizeY-1]
        // Подставляем U_{sizeY-2} = alpha[sizeY-2]*U_{sizeY-1} + beta[sizeY-2]
        // Получаем: (A[sizeY-1]*alpha[sizeY-2] + B[sizeY-1]) * U_{sizeY-1} = F[sizeY-1] - A[sizeY-1]*beta[sizeY-2]

        System.out.println("\nШаг 3: Обратный ход для последнего слоя j=" + (sizeY-1));

        int last = sizeY - 1;

        // Вычисляем A_last * alpha_{last-1}
        double[][] A_alpha_last = new double[sizeX][sizeX];
        for (int i = 0; i < sizeX; i++) {
            for (int k = 0; k < sizeX; k++) {
                A_alpha_last[i][k] = 0;
                for (int l = 0; l < sizeX; l++) {
                    A_alpha_last[i][k] += A[last][i][l] * alpha[last-1][l][k];
                }
            }
        }

        // B_last_tilda = A_alpha_last + B[last]
        double[][] B_last_tilda = new double[sizeX][sizeX];
        for (int i = 0; i < sizeX; i++) {
            for (int k = 0; k < sizeX; k++) {
                B_last_tilda[i][k] = A_alpha_last[i][k] + B[last][i][k];
            }
        }

        // F_last_tilda = F[last] - A[last]*beta[last-1]
        double[] F_last_tilda = new double[sizeX];
        for (int i = 0; i < sizeX; i++) {
            F_last_tilda[i] = F[last][i];
            for (int k = 0; k < sizeX; k++) {
                F_last_tilda[i] -= A[last][i][k] * beta[last-1][k];
            }
        }

        // Вычисляем inv(B_last_tilda)
        double[][] invB_last = invertMatrix(B_last_tilda);
        if (invB_last == null) {
            System.err.println("Ошибка: матрица B_last_tilda вырождена");
            return null;
        }

        // Решение для последнего слоя
        double[][] U = new double[sizeY][sizeX];
        for (int i = 0; i < sizeX; i++) {
            U[last][i] = 0;
            for (int k = 0; k < sizeX; k++) {
                U[last][i] += invB_last[i][k] * F_last_tilda[k];
            }
        }

        // Шаг 4: Обратный ход для остальных слоев
        // Uj = alpha[j]*U_{j+1} + beta[j]

        System.out.println("\nШаг 4: Обратный ход для j=" + (sizeY-2) + "..0");

        for (int j = sizeY - 2; j >= 0; j--) {
            for (int i = 0; i < sizeX; i++) {
                U[j][i] = beta[j][i];
                for (int k = 0; k < sizeX; k++) {
                    U[j][i] += alpha[j][i][k] * U[j+1][k];
                }
            }
        }

        System.out.println("\n=== Матричная прогонка завершена ===\n");

        return U;
    }

    /**
     * Обращение матрицы методом Гаусса
     * @param matrix исходная матрица
     * @return обратная матрица или null, если матрица вырождена
     */
    private double[][] invertMatrix(double[][] matrix) {
        int n = matrix.length;
        double[][] augmented = new double[n][2*n];

        // Создаем расширенную матрицу [A | I]
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                augmented[i][j] = matrix[i][j];
            }
            augmented[i][n + i] = 1.0;
        }

        // Прямой ход метода Гаусса
        for (int i = 0; i < n; i++) {
            // Поиск главного элемента
            int maxRow = i;
            for (int k = i + 1; k < n; k++) {
                if (Math.abs(augmented[k][i]) > Math.abs(augmented[maxRow][i])) {
                    maxRow = k;
                }
            }

            // Проверка на вырожденность
            if (Math.abs(augmented[maxRow][i]) < 1e-12) {
                System.err.printf("Матрица вырождена: нулевой элемент на позиции [%d][%d]\n", maxRow, i);
                return null;
            }

            // Перестановка строк
            double[] temp = augmented[i];
            augmented[i] = augmented[maxRow];
            augmented[maxRow] = temp;

            // Нормализация строки i
            double pivot = augmented[i][i];
            for (int j = 0; j < 2 * n; j++) {
                augmented[i][j] /= pivot;
            }

            // Вычитание из остальных строк
            for (int k = 0; k < n; k++) {
                if (k != i) {
                    double factor = augmented[k][i];
                    for (int j = 0; j < 2 * n; j++) {
                        augmented[k][j] -= factor * augmented[i][j];
                    }
                }
            }
        }

        // Извлекаем обратную матрицу
        double[][] inverse = new double[n][n];
        for (int i = 0; i < n; i++) {
            System.arraycopy(augmented[i], n, inverse[i], 0, n);
        }

        return inverse;
    }

    /**
     * Упрощенная версия для тестирования (если матрицы диагональные)
     * Используйте этот метод, если ваша система допускает скалярную прогонку
     */
    public double[][] simpleMatrixProgonka(double[][] A, double[][] B, double[][] C,
                                           double[][] F, int sizeY, int sizeX) {
        // Если матрицы диагональные, можно использовать скалярную прогонку
        double[] alpha = new double[sizeY];
        double[] beta = new double[sizeY];

        // Прямой ход
        alpha[0] = -C[0][0] / B[0][0];
        beta[0] = F[0][0] / B[0][0];

        for (int j = 1; j < sizeY - 1; j++) {
            double denominator = B[j][0] + A[j][0] * alpha[j-1];
            alpha[j] = -C[j][0] / denominator;
            beta[j] = (F[j][0] - A[j][0] * beta[j-1]) / denominator;
        }

        // Обратный ход
        double[][] U = new double[sizeY][sizeX];
        int last = sizeY - 1;
        double denominator = B[last][0] + A[last][0] * alpha[last-1];
        U[last][0] = (F[last][0] - A[last][0] * beta[last-1]) / denominator;

        for (int j = last - 1; j >= 0; j--) {
            U[j][0] = alpha[j] * U[j+1][0] + beta[j];
            // Копируем на все x (если решение не зависит от x)
            for (int i = 1; i < sizeX; i++) {
                U[j][i] = U[j][0];
            }
        }

        return U;
    }
    private double determinant(double[][] matrix) {
        int n = matrix.length;
        if (n == 1) return matrix[0][0];
        if (n == 2) return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
        if (n == 3) {
            return matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
                    - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
                    + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);
        }
        return 1; // для больших матриц используем другой метод
    }
}
