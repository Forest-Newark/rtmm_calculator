import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.interpolation.*;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

/**
 * Created by Forest on 6/19/17.
 */

// http://introcs.cs.princeton.edu/java/22library/Matrix.java.html
public class helper {

    public static double linearInterp(double[] x, double[] y, double xi) {
        // return linear interpolation of (x,y) on xi


        LinearInterpolator li = new LinearInterpolator();
        PolynomialSplineFunction psf = li.interpolate(x, y);
        double yi = psf.value(xi);
        return yi;
    }


    public static double univariateInterpolator(double[] x, double[] y, double xi){

        double xPoints[] = x;
        double yPoints[] = y;
        UnivariateInterpolator interpolator = new SplineInterpolator();
        UnivariateFunction function = interpolator.interpolate(xPoints, yPoints);
        double interpolatedY = function.value(xi);
        return interpolatedY;

    }

    // return c = a + b
    public static double[][] add(double[][] a, double[][] b) {
        int m = a.length;
        int n = a[0].length;
        double[][] c = new double[m][n];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                c[i][j] = a[i][j] + b[i][j];
        return c;
    }


    // return c = a - b
    public static double[][] subtract(double[][] a, double[][] b) {
        int m = a.length;
        int n = a[0].length;
        double[][] c = new double[m][n];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                c[i][j] = a[i][j] - b[i][j];
        return c;
    }


    //matrix multiplication
    // return c = a * b
    public static double[][] multiply(double[][] a, double[][] b) {
        int m1 = a.length;
        int n1 = a[0].length;
        int m2 = b.length;
        int n2 = b[0].length;
        if (n1 != m2) throw new RuntimeException("Illegal matrix dimensions.");
        double[][] c = new double[m1][n2];
        for (int i = 0; i < m1; i++)
            for (int j = 0; j < n2; j++)
                for (int k = 0; k < n1; k++)
                    c[i][j] += a[i][k] * b[k][j];
        return c;
    }


    //Matrix Transposition to facilitate division of matrices
    // return B = A^T
    public static double[][] transpose(double[][] a) {
        int m = a.length;
        int n = a[0].length;
        double[][] b = new double[n][m];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                b[j][i] = a[i][j];
        return b;
    }




}
