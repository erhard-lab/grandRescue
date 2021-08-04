package gedi.statistics;

import jdistlib.Uniform;

import static jdistlib.Beta.*;

public class Rtbeta {

    public static void main(String[] args) {
       /* double[] a = qtbeta(new double[]{0.1,0.329131,0.3847324,0.263721,0.3432420,0.97459,0,1},0.0001,0.04,2,2, false);
        for(int i = 0; i < a.length; i++){
            System.out.println(a[i]);
        }*/

        System.out.println(rtbeta(1,0.0001,0.04,2,2)[0]);

        double upper = 0.04;
        double lower = 0.0001;
        double param1 = 0;
        double param2 = 0;

    }

    public static double[] rtbeta(int number, double lower, double upper, double param1, double param2) {
        return qtbeta(new Uniform(0, 1).random(number), lower, upper, param1, param2, false);
    }

    private static double[] qtbeta(double[] p, double lower, double upper, double param1, double param2, boolean log_p) {

        if (!log_p) {
            for (int i = 0; i < p.length; i++) {
                p[i] = Math.log(p[i]);
            }
        }

        double logdiff = logdiff(cumulative_raw(upper, param1, param2, true, true), cumulative_raw(lower, param1, param2, true, true));

        for (int i = 0; i < p.length; i++) {
            p[i] = p[i] + logdiff;
        }


        double[] lf = lse(p, cumulative_raw(lower, param1, param2, true, true));

        for (int i = 0; i < lf.length; i++) {
            lf[i] = quantile(lf[i], param1, param2, true, true);
        }

        return lf;

    }

    public static double[] lse(double[] u, double v) {

        //pmax(u,v)
        double[] m = new double[u.length];
        double[] out = new double[u.length];

        for (int i = 0; i < u.length; i++) {
            m[i] = Math.max(u[i], v);

        }


        for (int i = 0; i < m.length; i++) {
            double uVal;
            double vVal;
            double mVal = m[i];

            if (i >= u.length) {
                uVal = u[i % u.length];
                //System.out.println("lse: u.length < v.length! using u[i%u.length]");
            } else {
                uVal = u[i];
            }
            out[i] = Math.log(Math.exp(uVal - mVal) + Math.exp(v - mVal)) + mVal;

        }

        return out;
    }

    private static double logdiff(double l1, double l2) {
        return (l1 + Math.log1p(-Math.exp(-(l1 - l2))));
    }
}
