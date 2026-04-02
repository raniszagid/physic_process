package sem8;


import java.util.function.BinaryOperator;

public class Main {
    public static void main(String[] args) {
        Solver example = new Solver(
                (x, y) -> x*x*y,
                (x, y) -> x*y,
                1, 9,
                1, 9,
                32, 32,
                (x, y) -> x*y*y,
                1
        );
        //double[][] sourceMatrix = example.countSourceMatrix();
        //printMatrix(sourceMatrix);
        example.run();
    }

    public static void printMatrix(double[][] matrix) {
        for (int j = 0; j < matrix.length; j++) {
            for (int i = 0; i < matrix[1].length; i++) {
                System.out.printf("%.3f\t", matrix[j][i]);
            }
            System.out.print("\n");
        }
    }
}
