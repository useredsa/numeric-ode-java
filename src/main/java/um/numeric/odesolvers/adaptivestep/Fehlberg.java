package um.numeric.odesolvers.adaptivestep;

import um.numeric.odesolvers.AdaptiveStepSolver;
import um.numeric.odesolvers.InitialValueProblem;

public class Fehlberg implements AdaptiveStepSolver {
  private InitialValueProblem problem;
  private double tolerance;
  private double currentStep;
  private double minimumStep;
  private double maximumStep;

  public Fehlberg(
      InitialValueProblem problem, double tolerance, double minimumStep, double maximumStep) {
    this.problem = problem;
    this.tolerance = tolerance;
    this.minimumStep = minimumStep;
    this.maximumStep = maximumStep;
  }

  @Override
  public int getOrder() {
    return 4; // REVISE controversial
  }

  @Override
  public double getTolerance() {
    return tolerance;
  }

  @Override
  public void setTolerance(double tolerance) {
    this.tolerance = tolerance;
  }

  @Override
  public double getStepSize() {
    return currentStep;
  }

  @Override
  public void setStepSize(double stepSize) {
    this.currentStep = stepSize;
  }

  @Override
  public double iterate(double time, double[] state) {
    return hintedIterate(time, state, problem.getDerivative(time, state));
  }

  @Override
  public double hintedIterate(double time, double[] state, double[] derivative) {
    int n = problem.getDimension();
    double[] aux = new double[n];
    double[] rk5 = new double[n];

    double[] k1 = derivative;
    for (int i = 0; i < n; i++) {
      aux[i] = state[i] + currentStep * (1.0 / 4.0 * k1[i]);
    }
    double[] k2 = problem.getDerivative(time + currentStep / 4.0, aux);
    for (int i = 0; i < n; i++) {
      aux[i] = state[i] + currentStep * (3.0 / 32.0 * k1[i] + 9.0 / 32.0 * k2[i]);
    }
    double[] k3 = problem.getDerivative(time + 3.0 / 8.0 * currentStep, aux);
    for (int i = 0; i < n; i++) {
      aux[i] = state[i] + currentStep * (1932.0 / 2197.0 * k1[i] - 7200.0 / 2197.0 * k2[i] + 7296.0 / 2197.0 * k3[i]);
    }
    double[] k4 = problem.getDerivative(time + 12.0 / 13.0 * currentStep, aux);
    for (int i = 0; i < n; i++) {
      aux[i] = state[i] + currentStep * (439.0 / 216.0 * k1[i] - 8.0 * k2[i] + 3680.0 / 513.0 * k3[i] - 845.0 / 4104.0 * k4[i]);
    }
    double[] k5 = problem.getDerivative(time + currentStep, aux);
    for (int i = 0; i < n; i++) {
      aux[i] = state[i] + currentStep * (-8.0 / 27.0 * k1[i] + 2.0 * k2[i] - 3544.0 / 2565.0 * k3[i] + 1859.0 / 4104.0 * k4[i] - 11.0 / 40.0 * k5[i]);
    }
    double[] k6 = problem.getDerivative(time + 1.0 / 2.0 * currentStep, aux);

    for (int i = 0; i < n; i++) {
      aux[i] = state[i] + currentStep * (25.0 / 216.0 * k1[i] + 1408.0 / 2565.0 * k3[i] + 2197.0 / 4104.0 * k4[i] - 1.0 / 5.0 * k5[i]);
      rk5[i] = state[i] + currentStep * (16.0 / 135.0 * k1[i] + 6656.0 / 12825.0 * k3[i] + 28561.0 / 56430.0 * k4[i] - 9.0 / 50.0 * k5[i] + 2.0 / 55.0 * k6[i]);
    }

    double error = 0;
    for (int i = 0; i < n; i++) {
      error += (rk5[i] - aux[i]) * (rk5[i] - aux[i]);
    }
    error = Math.sqrt(error);
    double q = Math.pow((tolerance * currentStep) / (2 * error), 0.25);
    q = Math.min(4, Math.max(q, 0.1));
    // Update time first
    time += currentStep;
    currentStep *= q;

    if (error < tolerance * currentStep) {
      for (int i = 0; i < n; i++) {
        state[i] = rk5[i];
      }
      return time;
    }
    return Double.NaN;
  }
}
