import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.interpolation.*;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

/**
 * Created by Forest on 6/19/17.
 */

//
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




}
