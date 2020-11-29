package um.numeric.odesolvers.fixedstep;

import um.numeric.odesolvers.FixedStepSolver;
import um.numeric.odesolvers.InitialValueProblem;

public class Euler implements FixedStepSolver {
  private InitialValueProblem problem;

  public Euler(InitialValueProblem problem) {
    this.problem = problem;
  }

  @Override
  public int getOrder() {
    return 1;
  }

  @Override
  public void step(double deltaTime, double time, double[] state) {
    hintedStep(deltaTime, time, state, problem.getDerivative(time, state));
  }

  @Override
  public void hintedStep(double deltaTime, double time, double[] state, double[] derivative) {
    int n = problem.getDimension();
    for (int i = 0; i < n; i++) {
      state[i] += deltaTime * derivative[i];
    }
  }
}
