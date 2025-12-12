package lab1;

public class MatrixSystem {
    public double[] a; // поддиагональ
    public double[] b; // главная диагональ
    public double[] c; // наддиагональ
    public double[] f; // правая часть
    public double[] expectedX; // ожидаемое решение
    public double t;



    public MatrixSystem(double[] f) {
        int n = f.length;
        this.a = new double[n];
        this.b = new double[n];
        this.c = new double[n];
        this.expectedX = new double[n];
        this.f = f;
    }
    public MatrixSystem(double t, double[] a, double[] b, double[] c, double[] f) {
        this.t = t;
        this.a = a;
        this.b = b;
        this.c = c;
        this.f = f;
        this.expectedX = new double[f.length];
    }

    public MatrixSystem(double[] a, double[] b, double[] c, double[] f, double[] expectedX) {
        this.a = a;
        this.b = b;
        this.c = c;
        this.f = f;
        this.expectedX = expectedX;
    }

    public void print() {
        int n = b.length;
        System.out.printf("Трехдиагональная матрица (%.2f):\n", t);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (j == i - 1) {
                    System.out.printf("%6.1f ", a[i]);
                } else if (j == i) {
                    System.out.printf("%6.1f ", b[i]);
                } else if (j == i + 1) {
                    System.out.printf("%6.1f ", c[i]);
                } else {
                    System.out.printf("%6.1f ", 0.0);
                }
            }
            System.out.printf(" | %6.1f\n", f[i]);
        }
    }
}
