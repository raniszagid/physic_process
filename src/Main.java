import java.util.ArrayList;
import java.util.List;

public class Main {
    static int gridSize = 8;
    public static void main(String[] args) {
        Runner case1 = classic();
        Runner case2 = new Runner(
          gridSize, 1, 14,
          2,
          k -> 1d,
                u -> u, q -> 1d
        );
        Runner case3 = new Runner(
                gridSize, 1, 14,
                2,
                k -> k,
                u -> u,
                q -> 1d
        );
        Runner case4 = new Runner(gridSize, 1, 14,
                2,
                k -> k,
                u -> u*u,
                q -> 1d);
        variousDivisions(case1);
        variousDivisions(case2);
        variousDivisions(case3);
        variousDivisions(case4);
    }
    public static void variousDivisions(Runner runner) {
        int n = 4;
        List<Double> stepSize = new ArrayList<>();
        List<Double> errorsX = new ArrayList<>();
        List<Double> maxNevyazka = new ArrayList<>();
        List<Double> maxInnerNevyazka = new ArrayList<>();
        while (n <= 128) {
            runner.n = n;
            runner.run();
            stepSize.add((runner.rR - runner.rL) / runner.n);
            errorsX.add(runner.maxError);
            maxNevyazka.add(runner.maxN);
            maxInnerNevyazka.add(runner.innerMaxN);
            runner.clear();
            n *= 2;
        }
        printL(stepSize);
        printList(errorsX);
        printList(maxNevyazka);
        printList(maxInnerNevyazka);
    }

    private static Runner classic() {
        return new Runner(
                gridSize, 1, 14,
                2,
                k -> 1d,
                u -> 1d,
                q -> 1d
        );
    }

    private static void printList(List<Double> list) {
        for (int i = 0; i < list.size(); i++) {
            System.out.printf("%.3e\t", list.get(i));
        }
        System.out.println();
    }
    private static void printL(List<Double> list) {
        for (int i = 0; i < list.size(); i++) {
            System.out.printf("%.3f\t", list.get(i));
        }
        System.out.println();
    }
}
