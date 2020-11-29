package um.numeric.odesolvers;

public interface AdaptiveStepSolver {
  public int getOrder();

  public double getTolerance();

  public void setTolerance(double tolerance);

  public double getStepSize();

  public void setStepSize(double stepSize);

  public double iterate(double time, double[] state);

  public double hintedIterate(double time, double[] state, double[] derivative);
}
