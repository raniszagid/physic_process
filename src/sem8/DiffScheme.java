package sem8;

public class DiffScheme {
    double[][] a;
    double[][] b;
    double[][] c;
    double[][] d;
    double[][] e;
    double[][] f;

    public DiffScheme(double[][] a, double[][] b, double[][] c, double[][] d, double[][] e, double[][] f) {
        this.a = a;
        this.b = b;
        this.c = c;
        this.d = d;
        this.e = e;
        this.f = f;
    }

    public void print() {
        printArray(a);
        printArray(b);
        printArray(c);
        printArray(d);
        printArray(e);
        printArray(f);
    }

    private void printArray(double[][] arr) {
        for (int j = 0; j < arr.length; j++) {
            for (int i = 1; i < arr[0].length; i++) {
                System.out.printf("%.3f\t", arr[j][i]);
            }
            System.out.println();
        }
        System.out.println();
    }
}
