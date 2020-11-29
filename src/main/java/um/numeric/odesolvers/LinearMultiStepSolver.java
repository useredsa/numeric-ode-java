package um.numeric.odesolvers;

public interface LinearMultiStepSolver {
  public int getOrder();

  public int getNeededSteps();

  public void step(MultiStepNumericalSolution solution);
}
