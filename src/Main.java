import java.util.ArrayList;
import java.util.List;

public class Main {
    static int gridSize = 8;
    public static void main(String[] args) {
        Runner case1 = new Runner(
                gridSize, 1, 9,
                2,1,
                k -> 1d,
                u -> 1d,
                q -> 1d
        );
        Runner case2 = new Runner(gridSize, 1, 9,
                2, 2,
                k -> k,
                u -> u*u,
                q -> 1d);
        case1.run();
        case2.run();
    }
}
