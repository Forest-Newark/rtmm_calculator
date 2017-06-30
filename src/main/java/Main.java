
import org.ejml.*;

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
        final  double MG = 92;

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
        double fs = 2.32;

        //Identity matrix
        double[][] I_ = {{1,0,0},{0,1,0},{0,0,1}};


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
        double[] fat1 = {9.05,4.34,3.08,6.86,9.58,11.38,13.73};

        double[] fatday = {0,64,148,232,316,330,353,371};
        double[] fat = {9.05,9.05,4.34,3.08,6.86,9.58,11.38,13.73};


        double[] yfat1 = new double[308];

        for (int x = 0; x < time1.length;x++){
            yfat1[x] = helper.univariateInterpolator(fatday1,fat1,time1[x]);

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
            yw1[x] = helper.univariateInterpolator(WKGday1,WKGram1,time1[x]);

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



        //GENERATING measured protein mass nP DATA

        double BWb = ywg[0];
        double Fb = yfatg[0];
        double BM = 0.04 * BWb;
        double ECWb = ECW[0];
        double ECPb = 0.732 * BM + 0.01087 * ECWb;

        double P1 = 10661.0;
        double gLP = P1 / (yleang[0] - BM - ECPb - ECW1);
        P1 = gLP * (yleang[0] - BM - ECPb - ECW1);

        double[] nP = new double[371];
        for(int x = 0; x<yleang.length;x++){
            nP[x] = gLP * (yleang[x] - BM - ECPb - ECW1);
        }

        //DATA FROM EXCEL FILE Minnesota_Data_Hall2.xlsx

        double [] pae = {1750.54, 1750.54, 1750.54, 1750.54, 1750.54, 1750.54, 1750.54, 1750.54, 1750.54, 1750.54, 1750.54, 1750.54, 1750.54, 1750.54, 1750.54, 1750.54, 1750.54, 1750.54, 1750.54, 1750.54, 1750.54, 1750.54, 1750.54, 1750.54, 1750.54, 1750.54, 1750.54, 1750.54, 1750.54, 1750.54, 1750.54, 1750.54, 1750.54, 1750.54, 1750.54, 1750.54, 1750.54, 1750.54, 1750.54, 1750.54, 1750.54, 1750.54, 1750.47, 1750.31, 1750.07, 1749.75, 1749.39, 1748.97, 1748.52, 1748.08, 1747.7, 1747.36, 1747.06, 1746.79, 1746.53, 1746.3, 1746.29, 1746.62, 1747.19, 1747.78, 1748.37, 1748.96, 1749.54, 1707.1, 1639.94, 1585.56, 1539.62, 1499.85, 1464.93, 1433.98, 1406.33, 1381.48, 1359, 1338.55, 1319.85, 1302.65, 1286.74, 1271.93, 1258.06, 1245.02, 1232.69, 1220.97, 1209.8, 1199.1, 1188.79, 1178.8, 1169.12, 1159.69, 1150.49, 1141.5, 1132.69, 1124.17, 1115.99, 1108.1, 1100.45, 1093, 1085.74, 1078.63, 1071.55, 1064.42, 1057.25, 1050.08, 1042.92, 1035.77, 1028.63, 1021.57, 1014.6, 1007.73, 1000.92, 994.189, 987.517, 980.902, 974.37, 967.933, 961.577, 955.294, 949.078, 942.924, 936.826, 930.758, 924.7, 918.66, 912.642, 906.649, 900.681, 894.741, 888.745, 882.644, 876.472, 870.251, 863.994, 857.712, 851.414, 845.085, 838.716, 832.322, 825.912, 819.493, 813.072, 806.652, 800.121, 793.412, 786.576, 779.643, 772.637, 765.574, 758.468, 751.56, 744.995, 738.69, 732.592, 726.667, 720.892, 715.245, 709.594, 703.856, 698.067, 692.247, 686.409, 680.56, 674.707, 668.8, 662.811, 656.764, 650.678, 644.562, 638.427, 632.28, 626.254, 620.43, 614.761, 609.219, 603.784, 598.442, 593.181, 588.071, 583.15, 578.38, 573.737, 569.202, 564.76, 560.401, 555.961, 551.341, 546.595, 541.757, 536.847, 531.883, 526.878, 521.967, 517.238, 512.643, 508.158, 503.763, 499.445, 495.194, 490.987, 486.807, 482.655, 478.529, 474.428, 470.35, 466.293, 462.255, 458.233, 454.226, 450.235, 446.258, 442.294, 438.343, 434.337, 430.235, 426.064, 421.843, 417.581, 413.29, 408.974, 404.726, 400.6, 396.564, 392.6, 388.694, 384.836, 381.019, 377.264, 373.582, 369.958, 366.379, 362.839, 359.33, 355.847, 366.167, 379.801, 391.639, 402.364, 412.331, 421.739, 430.711, 439.326, 447.645, 455.711, 463.56, 471.222, 478.722, 486.081, 493.317, 500.447, 507.484, 514.44, 521.326, 528.152, 534.924, 541.651, 548.339, 554.992, 561.616, 568.216, 574.793, 581.352, 587.896, 594.426, 600.946, 607.455, 613.957, 620.453, 626.943, 633.428, 639.91, 646.389, 652.865, 659.34, 665.813, 672.285, 685.416, 703.04, 718.47, 732.546, 745.69, 758.137, 770.026, 781.454, 792.49, 803.189, 813.597, 823.751, 833.684, 843.424, 852.996, 862.422, 871.719, 880.904, 889.992, 898.995, 907.924, 916.787, 925.594, 934.352, 943.067, 951.744, 960.387, 969.002, 982.086, 998.057, 1012.44, 1025.87, 1038.66, 1050.96, 1062.88, 1074.49, 1085.82, 1096.93, 1107.83, 1118.55, 1129.12, 1139.55, 1151.5, 1164.35, 1177.76, 1191.68, 1206.09, 1220.98, 1236.37, 1252.16, 1268.15, 1284.3, 1300.59, 1317.06, 1333.71, 1350.56, 1366.48, 1380.99, 1394.64, 1407.68, 1420.22, 1432.33, 1444.01, 1456.73, 1470.99, 1485.98, 1501.36, 1516.94, 1532.68, 1548.53, 1563.44, 1577.14, 1590.26, 1603.04, 1615.6, 1627.97, 1640.18, 1651.78, 1662.63, 1672.93, 1682.78, 1692.22, 1701.26, 1709.93, 1718.88, 1728.31, 1737.83, 1747.28, 1756.6, 1765.75, 1774.69, 1783.51, 1792.31, 1801.11, 1809.89, 1818.66, 1827.4, 1836.12};

        double[] time_ = {0, 7, 14, 21, 28, 35, 42, 49, 56, 63, 64, 70, 77, 84, 91, 98, 105, 112, 119, 126, 133, 140, 147, 154, 161, 168, 175, 182, 189, 196, 203, 210, 217, 224, 231, 232, 238, 245, 252, 259, 266, 273, 274, 280, 287, 294, 301, 302, 308, 315, 322, 329, 336, 343, 350, 357, 364, 371};
        double[] CI_= {1811.51, 1811.51, 1811.51, 1811.51, 1811.51, 1811.51, 1811.51, 1800.87, 1799.82, 1833.85, 1117.49, 1117.49, 1117.49, 1112.78, 1102.59, 1129.3, 1120.14, 1126.91, 1142.58, 1149.61, 1126.84, 1095.69, 1018.68, 1036.09, 1003.07, 945.148, 948.308, 990.685, 951.34, 984.493, 1008.49, 1031.45, 1007.97, 1047.44, 1108.66, 1671.76, 1671.76, 1671.76, 1671.76, 1671.76, 1671.76, 1671.76, 2050.53, 2050.53, 2050.53, 2050.53, 2050.53, 2257.02, 2257.02, 2257.02, 2442.57, 2715.02, 2337.5, 2846.34, 2668.86, 2257.47, 2246.94, 2217.49} ;
        double[] FI_ = {1332.52, 1332.52, 1332.52, 1332.52, 1332.52, 1332.52, 1332.52, 1324.69, 1323.92, 1348.96, 300.098, 300.098, 300.098, 299.581, 298.462, 301.394, 300.388, 301.132, 302.853, 303.625, 301.125, 297.704, 289.246, 291.159, 287.532, 281.171, 281.518, 286.172, 281.851, 285.492, 288.127, 290.649, 288.07, 292.405, 299.128, 480.488, 480.488, 480.488, 480.488, 480.488, 480.488, 480.488, 721.747, 721.747, 721.747, 721.747, 721.747, 811.948, 811.948, 811.948, 2019.82, 2355.34, 2067.59, 1728.89, 1621.09, 1371.2, 1067.48, 1053.49} ;
        double[] PI_ = {457.38, 457.38, 457.38, 457.38, 457.38, 457.38, 457.38, 454.692, 454.427, 463.021, 209.461, 209.461, 207.25, 206.392, 204.536, 209.4, 207.732, 208.965, 211.821, 213.102, 208.954, 203.279, 189.249, 192.421, 186.406, 175.854, 176.43, 184.15, 176.982, 183.021, 187.393, 191.576, 187.298, 194.49, 205.64, 306.499, 306.499, 306.499, 306.499, 306.499, 306.499, 306.499, 403.69, 403.69, 403.69, 403.69, 403.69, 458.198, 458.198, 458.198, 756.779, 841.193, 687.499, 695.772, 652.388, 551.826, 418.035, 412.557};

        //Generate missing data of above

        double[] CI = new double[time.length];

        for (int x = 0; x < time.length;x++){
            CI[x] = helper.univariateInterpolator(time_,CI_,time[x]);

        }

        double[] FI = new double[time.length];

        for (int x = 0; x < time.length;x++){
            FI[x] = helper.univariateInterpolator(time_,FI_,time[x]);

        }

        double[] PI = new double[time.length];

        for (int x = 0; x < time.length;x++){
            PI[x] = helper.univariateInterpolator(time_,PI_,time[x]);

        }


        //Adjustments for Thermic Effect of Feeding TEF

        double[] mCI = new double[CI.length];

        for(int x = 0; x < CI.length;x++){
            mCI[x] = Digest_C*(1 - alfaC) * CI[x];
        }

        double[] mFI = new double[FI.length];

        for(int x = 0; x < FI.length;x++){
            mFI[x] = Digest_F*(1 - alfaF) * FI[x];
        }


        double[] mPI = new double[PI.length];

        for(int x = 0; x < PI.length;x++){
            mPI[x] = Digest_P*(1 - alfaP) * PI[x];
        }


        //Initial/Baseline value calculations

        double Lb = BWb - Fb;
       //Two calculations done above... "Extraceullar water" & "extracellular protein mass

        //baseline protein mass
        double Pb = P1;

        //baseline fat intake
        double mFIb = mFI[0];

        //baseline carbohydrate intake
        double mCIb = mCI[0];

        //baseline protein inake
        double mPIb = mPI[0];

        //adaptive thermogeneises dimensionless
        double T = 0;

        //metabolizible energy intake in Kcal
        double MEIb = mFIb + mCIb + mPIb;

        //physical activity expenditure at steady state
        double PAEb = pae[0];

        //De novo lipogenesis at steady state
        double sigmab = mCIb / (1 + Math.pow(2,4));

        //Net gluconeogenesis from amino acids at steady state
        double GNGPb = 100;

        double day = 1;

        //total energy expenditure at baseline
        double TEEb = MEIb;

        //daily PAE
        double[] PAE = pae;

        //new corrected PAE
        PAEb = PAE[0];

        for (int x=0; x < 64;x++){
            PAE[x] = PAEb;
        }

        //calculation of the bias value dEc
        double dEc = RMR[0] + PAE[0] - TEEb;

        //calculation of total energy expenditure
        double[] TEE = new double[RMR.length];
        for (int x = 0;x < RMR.length;x++){
            TEE[x] = RMR[x] -dEc + PAE[x];
        }

        //Protein Oxidation at baseline
        double ProtOxb = mPIb - GNGPb;

        //Non-protein energy expenditure
        double Enpb = TEEb - ProtOxb;

        //Protein burning energy expenditure
        double Epb = ProtOxb;

        double FatOxb = rFFA * mFIb + sigmab;

        double CarbOxb = Enpb - FatOxb;

        double Abb = 10400;

        double Bbb = Lb - Abb * Math.log(Fb);

        //Initial R ratio adjusted
        double Ro = (RoF * TEEb - FatOxb * RoF) / (RoL * FatOxb);

        //Initial R ratio unadjusted
        double Ab = Ro * Fb;

        double rr = Ro;

        double fFb = FatOxb/Enpb;

        double[][] IBC = {{1,rGF,1},{0,rFFA,0},{0,0,1}};

        //double[][] OBC = {{-(1-fFb),0,-1},{0,-fFb,0},{0,0.-1}};
        double[][] OBC = new double[3][3];
        OBC[0][0] = -(1-fFb);
        OBC[0][1] = 0;
        OBC[0][2] = -1;

        OBC[1][0] = 0;
        OBC[1][1] = -(1-fFb);
        OBC[1][2] = 0;

        OBC[2][0] = 0;
        OBC[2][1] = 0;
        OBC[2][2] = -1;






        //TODO: double check this matrix setup...
        //double[][] EI = {{mCIb},{mFIb},{mPIb}};
        double[][] EI = new double[3][371];
        EI[0][0] = mCIb;
        EI[1][0] = mFIb;
        EI[2][0] = mPIb;


        //double[][] EE = {{Enpb},{Enpb},{Epb}};
        double[][] EE = new double[3][371];
        EE[0][0] = Enpb;
        EE[1][0] = Enpb;
        EE[2][0] = Epb;


        //double[][] ZBC = {{-sigmab},{sigmab},{-GNGPb}};
        double[][] ZBC = new double[3][371];
        ZBC[0][0] = -sigmab;
        ZBC[1][0] = sigmab;
        ZBC[2][0] = -GNGPb;



        double[][] RoB = {{RoL,rG-RoF,0},{0,RoF-rG,0},{0,0,RoP}};

        //Total energy balance
        double EBb = mCIb + mFIb + mPIb + TEEb;

        //Vector variables initial values

        double[] CarbOx = new double[371];
        for (int x = 0; x < 371; x++){
            CarbOx[x] = CarbOxb;
        }

        double[] FatOx = new double[371];
        for (int x = 0; x < 371; x++){
            FatOx[x] = FatOxb;
        }

        double[] ProtOx = new double[371];
        for (int x = 0; x < 371; x++){
            ProtOx[x] = ProtOxb;
        }

        double[] EB = new double[371];
        for (int x = 0; x < 371; x++){
            EB[x] = EBb;
        }

        double[] AAk = new double[371];
        for (int x = 0; x < 371; x++){
            AAk[x] = RoL;
        }

        double[] BBk = new double[371];
        for (int x = 0; x < 371; x++){
            BBk[x] = RoF;
        }

        double[] RLek = new double[371];
        for (int x = 0; x < 371; x++){
            RLek[x] = RoL;
        }

        double[] OK = new double[371];
        for (int x = 0; x < 371; x++){
            OK[x] = Bbb;
        }

        double[] ACK = new double[371];
        for (int x = 0; x < 371; x++){
            ACK[x] = Abb;
        }

        double[] F = new double[371];
        for (int x = 0; x < 371; x++){
            F[x] = Fb;
        }

        double[] zF = new double[371];
        for (int x = 0; x < 371; x++){
            zF[x] = Fb;
        }

        double[] L = new double[371];
        for (int x = 0; x < 371; x++){
            L[x] = Lb;
        }

        double[] zL = new double[371];
        for (int x = 0; x < 371; x++){
            zL[x] = Lb;
        }

        double[] P = new double[371];
        for (int x = 0; x < 371; x++){
            P[x] = Pb;
        }

        double[] zP = new double[371];
        for (int x = 0; x < 371; x++){
            zP[x] = Pb;
        }

        double[] RR = new double[371];
        for (int x = 0; x < 371; x++){
            RR[x] = Ro;
        }

        double[] mR = new double[371];
        for (int x = 0; x < 371; x++){
            mR[x] = Ro;
        }

        double[] DDL = new double[371];

        double[] DDF = new double[371];

        double[] DDP = new double[371];

        double[] DDE = new double[371];

        double[] zDDL = new double[371];

        double[] zDDF = new double[371];

        double[] zDDP = new double[371];

        double[] zDDE = new double[371];

        double[] zCI = new double[371];
        for (int x = 0; x < 371; x++){
            zCI[x] = mCIb;
        }

        double[] zFI = new double[371];
        for (int x = 0; x < 371; x++){
            zFI[x] = mFIb;
        }

        double[] zPI = new double[371];
        for (int x = 0; x < 371; x++){
            zPI[x] = mPIb;
        }


        //Filter initiations

        double[][] QBC = {{100,0,0},{0,100,0},{0,0,1}};
        double[][] RBC = {{100,0,0},{0,100,0},{0,0,1}};

        double[][] aPBC = QBC;
        double[][] pPBC = aPBC;

        double[][] S2 = {{0},{0},{0}};
        double[][] zdeltaBCC = {{0},{0},{0}};

         T = 0.0;

         double AT = 0.0;

         double zdeltaFat = 0.0;

         double pRLe = RoL;
         double aRLe = RoL;
         double rFs = RoF;
         double aPRLe = 1;
         double pPRLe = 1;
         double QRLe = 1;
         double RRLe = 1;
         double sigmak = sigmab;
         double omegak = 0;
         double Ack = Abb;
         double Ok = Bbb;
         double [][] pTAck = {{Ack},{Ok}};
         double [][] aPAck = {{100,0},{0,1}};
         double [][] pABk = {{RoL},{RoF}};
         double [][] aPABk = {{100,0},{0,1}};
         double [][] pPABk = {{100,0},{0,1}};

        //Forward Calculations

        // fill with days 2 - 371


        //Calculations of gluconeogenesis from protien

        double GNGCIo = (-0.5 * 100.0) / mCIb;
        double GNGPIo = (0.3 * 100.0) / mPIb;
        double GNGPo = 100.0/Pb;
        double GNGO = 20;


        //ANABOLIC vs KATABOLIC process

        //Fuck it we are doing it live...

        //INITIALIZE THE SHIT HERE

        //TODO: I changed them from 370 to 371 to match everything else.. lets see what that does
        double[] iFI = new double[371];
        double[] iCI = new double[371];
        double[] iPI = new double[371];

        //TODO: EI_ could be initiaalized after the loop...

        double[] PP = new double[371];
        double[] GNGP_ = new double[371];


        double[] TOTEE = new double[371];
        double[] eb = new double [371];

        double[] aDL = new double[371];
        double[] aDF = new double[371];
        double[] aDP = new double[371];
        double[] Ep = new double[371];

        double[] Enp = new double [371];
        double[] wfatOx = new double [371];
        double[] wcarbOx = new double[371];

        double[][] wOx_ = new double[3][371];
        double[] fF = new double[371];

        double[][] S1 = new double[3][371];


        //TODO: I NEED TO RETHINK X= 2 I DONT KNOW ITS CORRECT CUZ I AM DUMB>>>>

        for(int x = 0; x < days;x++){
           iFI[x] = mFI[x];
           iCI[x] = mCI[x];
           iPI[x] = mPI[x];
           //this line of code is strange.. indexes are probably not correct
           PP[x] = zP[x];
           GNGP_[x] = GNGPo * PP[x] + GNGCIo * iCI[x] + GNGPIo * iPI[x] + GNGO;


            fs = 2.32;

            //TODO: THis might need to move down the loop...
            if(zdeltaFat < 1){
                fs = 0.0;
            }

            TOTEE[x] = TEE[x];

            eb[x] = iCI[x] + iFI[x] + iPI[x] - TOTEE[x];

            EB[x] = eb[x];

            aDL[x] = (iCI[x] + rGF * iFI[x] + iPI[x] - sigmak + omegak - (pRLe *rr + rFs)) / pRLe;

            aDF[x] = (rFFA * iFI[x] + sigmak - (rFs * TOTEE[x])/(pRLe * rr + rFs)) / rFs;

            aDP[x] = gLP * (1 - gW) * aDL[x] - gLP * gW * aDF[x];

            Ep[x] = iPI[x]- GNGP_[x] - RoP * aDP[x];

            Enp[x] = TOTEE[x] - Ep[x];

            wfatOx[x] = (rFs * TOTEE[x]) /(pRLe * rr + rFs);

            wcarbOx[x] = Enp[x] - wfatOx[x];

            //TODO: check this shit

            wOx_[0][x] = wcarbOx[x];
            wOx_[1][x] = wfatOx[x];
            wOx_[2][x] = Ep[x];

            fF[x] = wOx_[1][x]/Enp[x];

            //IBC already initialized and declared

            //CHECK OBC ...

            //CHECK RoB...

            //THIS MIGHT NEED TO be +1 on the right half for all of these
            if(x > 0){
                EI[0][x] = iCI[x];
                EI[1][x] = iFI[x];
                EI[2][x] = iPI[x];

                EE[0][x] = Enp[x];
                EE[1][x] = Enp[x];
                EE[2][x] = Ep[x];

                ZBC[0][x] = -sigmak + omegak;
                ZBC[1][x] = sigmak;
                ZBC[2][x] = -GNGP_[x];

            }

        } //END OF DAY LOOPS

        double [][] IBCxEI = helper.multiply(IBC,EI);
        double[][] OBCxEE = helper.multiply(OBC,EE);

       //S1 = helper.add(ZBC,helper.add(helper.multiply(IBC,EI),helper.multiply(OBC,EE)));
        S1 = helper.add(helper.add(IBCxEI,OBCxEE),ZBC);

        S2 = helper.multiply(helper.transpose(RoB),S1);


        double[] zp = new double[371];
        double[] zf = new double[371];


        double[][] mBC = new double[3][371];
        double[] ddF = new double[371];
        double[] ddL = new double[371];
        double[] yF = new double[371];
        double[] mRR = new double[371];

        double[][] XLF = new double[2][371];

        for(int x = 0; x <days;x++){

            DDL[x] = S2[0][0];
            DDF[x] = S2[1][0];
            DDP[x] = S2[2][0];

            //TODO: These two are a little strange
            zf[x] = S2[1][0] + F[x];
            F[x] = zf[x];

            zp[x]=S2[2][0] + P[x];
            P[x] = zp[x];



            //MEASURED BODY COMPOSITION (Looping through days)

            mBC[0][x] = yleang[x];
            mBC[1][x] = yfatg[x];
            mBC[2][x] = P[x];

            if(x == 0){
                ddF[x] = 0;
                ddL[x] = 0;
                yF[x] = 0;
            }

            if(x > 0) {
                ddF[0] = yfatg[x] - yfatg[x-1];
                ddF[x] = yleang[x] - yleang[x-1];
                yF[x] = yfatg[x-1];

            }


            mRR[x] = Abb * ((2 * yF[x] - ddF[x]) / (2 * yF[x]*yF[x]));

            mR[x] = mRR[x];



            //CANONIC REPRESENTATION

            double alfa = 1;


            if(Math.abs(ddL[x]) > 2) {
                XLF[0][x] = ddL[x];
                XLF[1][x] = ddF[x];

            }



        }




    }




}
