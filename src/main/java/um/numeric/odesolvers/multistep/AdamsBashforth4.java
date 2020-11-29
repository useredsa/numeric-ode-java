package um.numeric.odesolvers.multistep;

import um.numeric.odesolvers.InitialValueProblem;
import um.numeric.odesolvers.LinearMultiStepSolver;
import um.numeric.odesolvers.MultiStepNumericalSolution;

public class AdamsBashforth4 implements LinearMultiStepSolver {
  private InitialValueProblem problem;

  public AdamsBashforth4(InitialValueProblem problem) {
    this.problem = problem;
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
  public void step(MultiStepNumericalSolution solution) {
    int n = problem.getDimension();
    int m = solution.size();
    double deltaTime = solution.getStep();
    double[] d0 = solution.getDerivative(m - 4);
    double[] d1 = solution.getDerivative(m - 3);
    double[] d2 = solution.getDerivative(m - 2);
    double[] d3 = solution.getDerivative(m - 1);
    double[] state = solution.getState(m - 1).clone();

    for (int i = 0; i < n; i++) {
      state[i] += deltaTime / 24 * (55 * d3[i] - 59 * d2[i] + 37 * d1[i] - 9 * d0[i]);
    }
    double time = solution.getTime(m - 1) + deltaTime;
    solution.add(time, state, problem.getDerivative(time, state));
  }
}
