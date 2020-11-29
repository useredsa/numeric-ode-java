package um.numeric.odesolvers;

public interface FixedStepSolver {
  public int getOrder();

  public void step(double deltaTime, double time, double[] state);

  public void hintedStep(double deltaTime, double time, double[] state, double[] derivative);
}
