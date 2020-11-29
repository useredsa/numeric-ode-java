package um.numeric.odesolvers.adaptivestep;

import um.numeric.odesolvers.AdaptiveStepSolver;
import um.numeric.odesolvers.FixedStepSolver;
import um.numeric.odesolvers.InitialValueProblem;

public class RichardsonExtrapolationOfFixedStep implements AdaptiveStepSolver {
  private InitialValueProblem problem;
  private FixedStepSolver baseMethod;
  private double tolerance;
  private double currentStep;
  private double minimumStep;
  private double maximumStep;

  public RichardsonExtrapolationOfFixedStep(
      InitialValueProblem problem,
      FixedStepSolver baseMethod,
      double tolerance,
      double minimumStep,
      double maximumStep) {
    this.problem = problem;
    this.baseMethod = baseMethod;
    this.tolerance = tolerance;
    this.minimumStep = minimumStep;
    this.maximumStep = maximumStep;
  }

  @Override
  public int getOrder() {
    return baseMethod.getOrder() + 1;
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
    double[] fullStep = state.clone();
    double[] halfStep = state.clone();
    baseMethod.hintedStep(currentStep, time, fullStep, derivative);
    baseMethod.hintedStep(currentStep / 2, time, halfStep, derivative);
    baseMethod.step(currentStep / 2, time + currentStep / 2, halfStep);

    int n = problem.getDimension();
    double error = 0;
    for (int i = 0; i < n; i++) {
      error += (fullStep[i] - halfStep[i]) * (fullStep[i] - halfStep[i]);
    }
    error = Math.sqrt(error);
    double q = Math.pow((tolerance * currentStep) / (2 * error), 1.0 / baseMethod.getOrder());
    q = Math.min(4, Math.max(q, 0.1));
    // Update time first
    time += currentStep;
    currentStep *= q;

    if (error < tolerance * currentStep) {
      for (int i = 0; i < n; i++) {
        state[i] = halfStep[i];
      }
      return time;
    }
    return Double.NaN;
  }
}
