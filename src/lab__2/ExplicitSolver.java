package lab__2;

import lab1.MatrixSystem;

public class ExplicitSolver {
    public static void main(String[] args) {
        EulerSolver eulerSolver = new EulerSolver(8, 10, 1, 9, 10,
                (r, t) -> 1d, // u
                (r, t) -> 1d, // k
                (r, t) -> 1d, // q
                2);
        eulerSolver.run3();
    }
}
