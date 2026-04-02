package sem8;


import java.util.function.BinaryOperator;

public class Main {
    public static void main(String[] args) {
        Solver example = new Solver(
                (x, y) -> 1d,
                (x, y) -> 1d,
                1, 9,
                1, 9,
                3, 3,
                (x, y) -> x*y,
                y -> 1d,
                y -> 1d,
                x -> 1d,
                x -> 1d,
                1
        );
        //double[][] sourceMatrix = example.countSourceMatrix();
        //printMatrix(sourceMatrix);
        example.run();
    }

    public static void printMatrix(double[][] matrix) {
        for (int j = 0; j < matrix.length; j++) {
            for (int i = 1; i < matrix[1].length; i++) {
                System.out.printf("%.3f\t", matrix[j][i]);
            }
            System.out.print("\n");
        }
    }
}
