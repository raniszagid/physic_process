package sem8;

public class MatrixSystem2D {
    double[][][] a;
    double[][][] b;
    double[][][] c;
    double[][] f;
    int sizeX;
    int sizeY;

    public MatrixSystem2D(double[][][] a, double[][][] b, double[][][] c, double[][] f, int sizeX, int sizeY) {
        this.a = a;
        this.b = b;
        this.c = c;
        this.f = f;
        this.sizeX = sizeX;
        this.sizeY = sizeY;
    }
}
