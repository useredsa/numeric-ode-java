package um.numeric.odesolvers.fixedstep;

import um.numeric.odesolvers.FixedStepSolver;
import um.numeric.odesolvers.InitialValueProblem;

public class RK4 implements FixedStepSolver {
  private InitialValueProblem problem;

  public RK4(InitialValueProblem problem) {
    this.problem = problem;
  }

  @Override
  public int getOrder() {
    return 4;
  }

  @Override
  public void step(double deltaTime, double time, double[] state) {
    hintedStep(deltaTime, time, state, problem.getDerivative(time, state));
  }

  @Override
  public void hintedStep(double deltaTime, double time, double[] state, double[] derivative) {
    int n = problem.getDimension();
    double[] auxState = new double[n];
    double h = deltaTime, h2 = deltaTime / 2;

    double[] k1 = derivative;
    for (int i = 0; i < n; i++) {
      auxState[i] = state[i] + h2 * k1[i];
    }
    double[] k2 = problem.getDerivative(time + h2, auxState);
    for (int i = 0; i < n; i++) {
      auxState[i] = state[i] + h2 * k2[i];
    }
    double[] k3 = problem.getDerivative(time + h2, auxState);
    for (int i = 0; i < n; i++) {
      auxState[i] = state[i] + h * k3[i];
    }
    double[] k4 = problem.getDerivative(time + h, auxState);

    for (int i = 0; i < n; i++) {
      state[i] += h / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
    }
  }
}
