package um.numeric.odesolvers;

public interface InitialValueProblem {
  public int getDimension();

  public double getInitialTime();

  public double[] getInitialState();

  public double[] getDerivative(double time, double[] state);
}
