import org.apache.commons.math3.analysis.interpolation.HermiteInterpolator;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.interpolation.SmoothingPolynomialBicubicSplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

/**
 * Created by Forest on 6/19/17.
 */
public class helper {

    public static double linearInterp(double[] x, double[] y, double xi) {
        // return linear interpolation of (x,y) on xi


        LinearInterpolator li = new LinearInterpolator();
        PolynomialSplineFunction psf = li.interpolate(x, y);
        double yi = psf.value(xi);
        return yi;
    }

    public static double Hermite(double[] x, double[] y, double xi){

        HermiteInterpolator interpolator = new HermiteInterpolator();

        return 0.0;
    }



}
