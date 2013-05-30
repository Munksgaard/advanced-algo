import lpsolve.*;

public class Demo1 {

  public static void main(String[] args) {
    try {
      // Create a problem with 4 variables and 0 constraints
      LpSolve solver = LpSolve.makeLp(0, 4);
      solver.setDebug(false);

      // add constraints
      double[] arr = new double[] {0, 3,2,2,1};
      double[] arr2 = new double[] {0, 0,4,3,1};
      solver.addConstraint(arr, LpSolve.LE, 4);
      solver.addConstraint(arr2, LpSolve.GE, 3);
//      solver.strAddConstraint("3 2 2 1", LpSolve.LE, 4);
//      solver.strAddConstraint("0 4 3 1", LpSolve.GE, 3);

      // set objective function
      double[] objFn = new double[] {0, 2,3,-2,3};
      solver.setObjFn(objFn);
//      solver.strSetObjFn("2 3 -2 3");

      // solve the problem
      solver.solve();

      // print solution
      System.out.println("Value of objective function: " + solver.getObjective());
      double[] var = solver.getPtrVariables();
      for (int i = 0; i < var.length; i++) {
        System.out.println("Value of var[" + i + "] = " + var[i]);
      }

      // delete the problem and free memory
      solver.deleteLp();
    }
    catch (LpSolveException e) {
       e.printStackTrace();
    }
  }

}

