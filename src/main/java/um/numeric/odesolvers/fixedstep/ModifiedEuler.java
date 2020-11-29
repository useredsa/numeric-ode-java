package um.numeric.odesolvers.fixedstep;

import um.numeric.odesolvers.FixedStepSolver;
import um.numeric.odesolvers.InitialValueProblem;

public class ModifiedEuler implements FixedStepSolver {
  private InitialValueProblem problem;

  public ModifiedEuler(InitialValueProblem problem) {
    this.problem = problem;
  }

  @Override
  public int getOrder() {
    return 2;
  }

  @Override
  public void step(double deltaTime, double time, double[] state) {
    hintedStep(deltaTime, time, state, problem.getDerivative(time, state));
  }

  @Override
  public void hintedStep(double deltaTime, double time, double[] state, double[] derivative) {
    int n = problem.getDimension();
    double[] state2 = new double[n];
    for (int i = 0; i < n; i++) {
      state2[i] = state[i] + deltaTime * derivative[i];
    }
    double[] derivative2 = problem.getDerivative(time + deltaTime, state2);
    for (int i = 0; i < n; i++) {
      state[i] += deltaTime / 2 * (derivative[i] + derivative2[i]);
    }
  }
}
