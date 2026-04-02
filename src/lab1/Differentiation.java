package lab1;

import java.util.function.BinaryOperator;
import java.util.function.UnaryOperator;

public class Differentiation {

    public static double derivative(UnaryOperator<Double> f, double x) {
        final double eps = 1e-10;
        final int MAX_ITERATIONS = 40;

        if (f == null) {
            throw new IllegalArgumentException("Callback function is null");
        }

        if (eps <= 0.0) {
            throw new IllegalArgumentException("Epsilon is zero or negative");
        }

        int iteration = 0;

        double delta = 0.0625;
        double derivativePrevious = (f.apply(x + delta) - f.apply(x)) / delta;
        double derivativeCurrent = 0.0;

        delta *= 0.25;
        derivativeCurrent = (f.apply(x + delta) - f.apply(x)) / delta;

        while ((Math.abs(derivativeCurrent - derivativePrevious) > eps)
                && (iteration < MAX_ITERATIONS)) {
            derivativePrevious = derivativeCurrent;
            delta *= 0.25;
            derivativeCurrent = (f.apply(x + delta) - f.apply(x)) / delta;
            ++iteration;
        }

        if (iteration == MAX_ITERATIONS) {
            throw new RuntimeException("Can't find a derivative for MAX_ITERATIONS times");
        }

        return derivativeCurrent;
    }

    public static double derivative2(UnaryOperator<Double> f, double x) {
        final double eps = 1e-10;
        final int MAX_ITERATIONS = 20;

        if (f == null) {
            throw new IllegalArgumentException("Callback function is null");
        }

        if (eps <= 0.0) {
            throw new IllegalArgumentException("Epsilon is zero or negative");
        }

        int iteration = 0;

        double delta = 0.0625;
        double derivativePrevious = (f.apply(x + 2.0 * delta) - 2.0 * f.apply(x + delta) + f.apply(x))
                / (delta * delta);
        double derivativeCurrent = 0.0;

        delta *= 0.25;
        derivativeCurrent = (f.apply(x + 2.0 * delta) - 2.0 * f.apply(x + delta) + f.apply(x))
                / (delta * delta);

        while ((Math.abs(derivativeCurrent - derivativePrevious) > eps)
                && (iteration < MAX_ITERATIONS)) {
            derivativePrevious = derivativeCurrent;
            delta *= 0.25;
            derivativeCurrent = (f.apply(x + 2.0 * delta) - 2.0 * f.apply(x + delta) + f.apply(x))
                    / (delta * delta);
            ++iteration;
        }

        if (iteration == MAX_ITERATIONS) {
            throw new RuntimeException("Can't find a derivative for MAX_ITERATIONS times");
        }

        return derivativeCurrent;
    }

    public static double derivativeByT(BinaryOperator<Double> func, double r, double t) {
        UnaryOperator<Double> onlyT = x -> func.apply(r, x);
        return derivative(onlyT, t);
    }

    public static double derR(BinaryOperator<Double> func, double r, double t) {
        UnaryOperator<Double> onlyR = x -> func.apply(x, t);
        return derivative(onlyR, r);
    }

    public static double der2R(BinaryOperator<Double> func, double r, double t) {
        UnaryOperator<Double> onlyR = x -> func.apply(x, t);
        return derivative2(onlyR, r);
    }
    public static double der2T(BinaryOperator<Double> func, double r, double t) {
        UnaryOperator<Double> onlyT = x -> func.apply(r, x);
        return derivative2(onlyT, t);
    }
}