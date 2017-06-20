import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

import java.util.Arrays;

/**
 * Created by Forest on 6/19/17.
 */

//https://stackoverflow.com/questions/36523396/interpolate-function-using-apache-commons-math/36523685


public class Main {


    public static void main(String[] args) {


    /*
    Constant Variables as defined on line 24 - 51 of Test_Neill_1
     */

        //energy density of carbohydrate
        final double RoC = 4.2;

        //energy density of glucose or glycogen Kcal/g
        final double RoG = 4.187;

        //energy density of fat Kcal/g
        final  double RoF = 9.4;

        //energy density of protien Kcal/g
        final  double RoP =4.7;

        //energy density lean body mass Kcal/g
        final  double RoL =1.8;

        //synthesis cost fat Kcal/g
        final  double nuF = .18;

        //synthesis cost glucose?
        final  double nuG = .21;

        //synthesis cost of protein
        final  double nuP = .86;

        //de novo lipogenesis efficiency
        final  double epsd = .8;

        //gluconeogenesis efficiency
        final  double epsg = .8;

        //protein degradation cost Kcal/g
        final  double epsP = .17;

        //short term thermic effect for fat
        final  double alfaF = .025;

        //short term thermic effect for carbohydrate
        final  double alfaC = .075;

        //short term thermic effect for protein
        final  double alfaP = .25;

        //molecular weight triglyceride g/mol
        final  double MTG =860;

        //molecular weight free fatty acid g/mol
        final  double MFFA = 273;

        //molecular weight of glycepRLe g/mol
        final  double MG = 98;

        //digestibility factors for carbohydrate
        final  double Digest_C= .95;

        //digestibility factors for Fat
        final  double Digest_F = .96;

        //digestibility factors for protein
        final  double Digest_P = .90;

        //molecule weight ratios
        final  double rFFA = (3* MFFA) / MTG;

        //molecule weight ratios
        final  double rG = (RoC * MG) / MTG;

        //molecule weight ratios
        final  double rGF = 1 - rFFA;

        //for reach gram of fat 2.32 Kcal carbohydrate energy is needed
        final  double fs = 2.32;

        //Identity matrix
        int [][] I_ = new int[][] {{1,0,0},{0,1,0},{0,0,1}};




        int days = 371;


        int[] time = new int[371];
        //populate "time" with 1 - 371
        for(int x = 0; x < 371; x++){

            time[x] = x+1;
        }


        //daily measured weight in grams
        double[] ywg = new double[371];

        //daily measured lean body mass
        double[] yleang = new double[371];

        //daily measured lean fat
        double[] yfatg = new double[371];




        //Generation of missing data with interpolation

        double time1[] = new double[308];
        for(int x = 0;x < 308; x++){
            time1[x] = x + 64;
        }

        double[] fatday1 = {64,148,232,316,330,353,371};
        double[] fat1 = {9.05,4.34,3.06,6.86,9.58,11.38,13.73};

        double[] fatday = {0,64,148,232,316,330,353,371};
        double[] fat = {9.05,9.05,4.34,3.06,6.86,9.58,11.38,13.73};


        for (double timepoint: time1){

            System.out.println(helper.linearInterp(fatday1,fat1,timepoint));

        }







//
//        double[] x = { 0, 50, 100 };
//        double[] y = { 0, 50, 200 };
//
//        LinearInterpolator interp = new LinearInterpolator();
//        PolynomialSplineFunction f = interp.interpolate(x, y);
//
////        System.out.println("Piecewise functions:");
////        Arrays.stream(f.getPolynomials()).forEach(System.out::println);
//
//        System.out.println(f.getPolynomials().length);
//
//        double value = f.value(70);
//        System.out.println("y for xi = 70: " + value);
//








    }




}