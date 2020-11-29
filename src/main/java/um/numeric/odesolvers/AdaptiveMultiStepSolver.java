package um.numeric.odesolvers;

public interface AdaptiveMultiStepSolver {
  public int getOrder();

  public int getNeededSteps();

  public double getTolerance();

  public void setTolerance(double tolerance);

  public double getStepSize();

  public void setStepSize(double stepSize);

  public double iterate(MultiStepNumericalSolution solution);
}
