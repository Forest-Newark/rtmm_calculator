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


        double[] yfat1 = new double[308];

        for (int x = 0; x < time1.length;x++){
            yfat1[x] = helper.linearInterp(fatday1,fat1,time1[x]);

        }

        double[] yfat = new double[371];
        for(int x = 0; x <371 ; x++){
            if(x < 64){
                yfat[x] = yfat1[0];
            }
            else {
                yfat[x] = yfat1[x - 64];
            }
        }


        //Checking yfat
//
//        for(int x = 0; x < yfat.length;x++){
//            System.out.println("Index - " + x + " number - " + yfat[x]);
//        }


        double[] WKGram1 = {67.608,55.783,51.692,58.158,64.467,69.698,70.824};
        double[] WKGday1 = {64,148,232,316,330,353,371};

        double[] yw1 = new double[308];

        for (int x = 0; x < time1.length;x++){
            yw1[x] = helper.linearInterp(WKGday1,WKGram1,time1[x]);

        }

        double[] yw = new double[371];
        for(int x = 0; x <371 ; x++){
            if(x < 64){
                yw[x] = yw1[0];
            }
            else {
                yw[x] = yw1[x - 64];
            }
        }

//        Checking yw
//
//        for(int x = 0; x < yw.length;x++){
//            System.out.println("Index - " + x + " number - " + yw[x]);
//        }


        //Measured fat trajectory grams
        for(int x = 0;x <yfatg.length;x++){
            yfatg[x] = yfat[x] * 1000.;

        }

        //Measured lean mass trajectory grams
        for(int x = 0;x <yleang.length;x++){
            yleang[x] = 1000 * (yw[x] - yfat[x]);

        }

        //Measured weight trajectory grams
        for(int x = 0;x <ywg.length;x++){
            ywg[x] = yw[x] * 1000.;

        }


        //GENERATING RMR DATA//

        double RMR1 = 1591.7;
        double height =172.0;
        double age = 35.0;
        double BMR0 = (10 * yw[0]) + (6.25 * height) - (5 * age);
        double cBMR0 = RMR1 - BMR0;


        //** double check that this should be an array of doubles and not a single value...
        double[] RMR = new double[371];
        for(int x = 0; x < yw.length;x++){
            RMR[x] = BMR0 = (10 * yw[x]) + (6.25 * height) - (5 * age) + cBMR0;
        }


        //GENERATING ECW DATA//

        double ECW1 =  11109.0;
        double gW = 0.191;
        double ECW0 = 1000 * ( (gW * yw[0]) + (.0957 * height) + (.025 * age) - 12.424);
        double cECW0 = ECW1 - ECW0;


        double[] ECW = new double[371];
        for(int x = 0; x < yw.length;x++){
            ECW[0] = 1000 * ( (gW * yw[x]) + (.0957 * height) + (.025 * age) - 12.424) + cECW0;
        }






















    }




}
