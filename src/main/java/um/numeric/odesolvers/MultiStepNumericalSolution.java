package um.numeric.odesolvers;

import java.util.List;

public interface MultiStepNumericalSolution {
  public int size();

  public int implicitSolutionSize();

  public int derivativesSize();

  public void add(double time, double[] state, double[] derivative);

  public void addDerivative(double time, double[] derivative);

  public void addState(double time, double[] state);

  public void remove(int i);

  public void removeDerivative(int i);

  public void removeState(int i);

  public double getTime(int i);

  public double[] getState(int i);

  public double[] getDerivative(int i);

  public double getStep();

  public List<Double> getTimes();

  public List<double[]> getPoints();
}
