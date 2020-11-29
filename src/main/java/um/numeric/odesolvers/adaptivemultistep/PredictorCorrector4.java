package um.numeric.odesolvers;

import um.numeric.odesolvers.multistep.AdamsBashforth4;

public class PredictorCorrector4 implements AdaptiveMultiStepSolver {
  private InitialValueProblem problem;
  private AdamsBashforth4 predictor;
  private double tolerance;
  private double currentStep;
  private double minimumStep;
  private double maximumStep;
  private boolean peceMode;

  public PredictorCorrector4(
      InitialValueProblem problem, double tolerance, double minimumStep, double maximumStep) {
    this.problem = problem;
    this.predictor = new AdamsBashforth4(problem);
    this.tolerance = tolerance;
    this.minimumStep = minimumStep;
    this.maximumStep = maximumStep;
    this.peceMode = false;
  }

  @Override
  public int getOrder() {
    return 4;
  }

  @Override
  public int getNeededSteps() {
    return 3;
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

  public boolean getPeceMode() {
    return peceMode;
  }

  public void setPeceMode(boolean peceMode) {
    this.peceMode = peceMode;
  }

  @Override
  public double iterate(MultiStepNumericalSolution solution) {
    // Predict before getting m
    predictor.step(solution);

    int n = problem.getDimension();
    int m = solution.size();
    double deltaTime = solution.getStep();
    double[] d0 = solution.getDerivative(m - 4);
    double[] d1 = solution.getDerivative(m - 3);
    double[] d2 = solution.getDerivative(m - 2);
    double[] d3 = solution.getDerivative(m - 1);
    // We change this reference is on purpose
    double[] predictor = solution.getState(m - 1);
    double[] corrector = predictor.clone();

    for (int i = 0; i < n; i++) {
      corrector[i] += deltaTime / 24 * (9 * d3[i] + 19 * d2[i] - 5 * d1[i] + d0[i]);
    }

    double error = 0;
    for (int i = 0; i < n; i++) {
      error += (corrector[i] - predictor[i]) * (corrector[i] - predictor[i]);
    }
    error = 19.0 * Math.sqrt(error) / 270.0;
    double q = 1.5 * Math.pow(tolerance * currentStep / error, 0.25);
    q = Math.min(4, Math.max(q, 0.1));

    if (error < tolerance * currentStep) {
      if (q < 1 / 10 || q > 10) {
        currentStep *= q;
      }
      for (int i = 0; i < n; i++) {
        predictor[i] = corrector[i];
      }
      if (peceMode) {
        double[] eval = problem.getDerivative(solution.getTime(m - 1), corrector);
        for (int i = 0; i < n; i++) {
          d3[i] = eval[i];
        }
      }
      return solution.getTime(m - 1);
    } else {
      currentStep *= q;
      solution.remove(m - 1);
      return Double.NaN;
    }
  }
}
