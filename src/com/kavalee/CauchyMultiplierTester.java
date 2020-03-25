package com.kavalee;


import com.sun.xml.internal.xsom.impl.scd.Iterators;
import org.junit.Test;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.Arrays;

public class CauchyMultiplierTester {
    public double[] generateRandomS(int n) {
        double[] s = new double[n];
        for (int i = 0; i < n; i++) {
            s[i] = (double) Math.random();
        }
        Arrays.sort(s);
        return s;
    }
    @Test
    public void test () {
        double[] s = generateRandomS(4);
        double[] t = generateRandomS(4);
        double[] g = generateRandomS(4);
        CauchyMultiplier cm = new CauchyMultiplier(s,t,g,10);
        System.out.println(Arrays.toString(new double[10]));
        System.out.println(Arrays.toString(cm.sortIntoCells(s,2)));
    }
    @Test
    public void testMultiplySlow () {
        double[] s = generateRandomS(1);
        double[] t = generateRandomS(1);
        double[] g = generateRandomS(1);
        System.out.println(Arrays.toString(s));
        System.out.println(Arrays.toString(t));
        System.out.println(Arrays.toString(g));
        CauchyMultiplier cm = new CauchyMultiplier(s,t,g,10);
        System.out.println(Arrays.toString(cm.multiplySlow()));
    }
    @Test
    public void testMultiply2Fast() {
        int n = 8;
        double[] s = generateRandomS(n);
        double[] t = generateRandomS(n);
        double[] g = generateRandomS(n);
        CauchyMultiplier cm = new CauchyMultiplier(s,t,g,10);
        double[] ff = cm.multiplyFast2();
        double[] ss = cm.multiplySlow();
        System.out.println(Arrays.toString(ff));
        System.out.println(Arrays.toString(ss));
        double[] e = new double[ss.length];
        for (int i = 0; i < ff.length; i++) {
            e[i] = ff[i] - ss[i];
        }
        System.out.println(Arrays.toString(e));
    }
    @Test
    public void testMultiply3Fast() {
        int n = 100000;
        double[] s = generateRandomS(n);
        double[] t = generateRandomS(n);
        double[] g = generateRandomS(n);
        CauchyMultiplier cm = new CauchyMultiplier(s,t,g,20);
        double[] ff = cm.multiplyFast3();
        double[] ss = cm.multiplySlow();
        System.out.println(Arrays.toString(ff));
        System.out.println(Arrays.toString(ss));
        double[] e = new double[ss.length];
        for (int i = 0; i < ff.length; i++) {
            e[i] = ff[i] - ss[i];
        }
        System.out.println(Arrays.toString(e));
    }
    @Test
    public void testMultiply3FastAlone() {
        int n = 2;
        double[] s = generateRandomS(n);
        double[] t = generateRandomS(n);
        double[] g = generateRandomS(n);
        CauchyMultiplier cm = new CauchyMultiplier(s,t,g,10);
        double[] ff = cm.multiplyFast3();
        System.out.println(Arrays.toString(ff));
    }
    @Test
    public void testUtils () {
        int n = 2;
        int q = 4;
        int ti = 3; int si = 2;
        //System.out.println((ti - si)/(double) q + 1f/(double) (2*q));
        System.out.println(CauchyMultiplier.nCk(5,3));
    }
    @Test
    public void generateTime() throws FileNotFoundException, UnsupportedEncodingException {
        int numTrials = 10;
        float deltaN = 1000;
        int maxN = 100000;
        PrintWriter writer = new PrintWriter("fast.csv", "UTF-8");

        for(int currN = 10; currN < maxN; currN += deltaN ) {
            long totalTime = 0;
            for(int i = 0; i < numTrials; i++) {
                double[] s = generateRandomS(currN);
                double[] t = generateRandomS(currN);
                double[] g = generateRandomS(currN);
                long ti = System.currentTimeMillis();
                CauchyMultiplier cm = new CauchyMultiplier(s,t,g,10);
                cm.multiplyFast3();
                long time = System.currentTimeMillis() - ti;
                totalTime += time;


            }
            writer.println(currN + "," + (float) totalTime/ (float) numTrials);

        }
        writer.close();
    }
    @Test
    public void generateError() throws FileNotFoundException, UnsupportedEncodingException {
        int numTrials = 10;
        float deltaN = 1.1f;
        int maxN = 10000;
        PrintWriter writer = new PrintWriter("error.csv", "UTF-8");

        for(int currN = 10; currN < maxN; currN *= deltaN ) {
            double totalMaxError = 0;
            for(int i = 0; i < numTrials; i++) {
                double[] s = generateRandomS(currN);
                double[] t = generateRandomS(currN);
                double[] g = generateRandomS(currN);
                CauchyMultiplier cm = new CauchyMultiplier(s,t,g,10);
                double[] ff = cm.multiplyFast3();
                double[] ss = cm.multiplySlow();
                double[] e = new double[ss.length];
                for (int j = 0; j < ff.length; j++) {
                    e[j] = Math.abs(ff[j] - ss[j])/Math.abs(ss[j]);
                }
                Arrays.sort(e);
                totalMaxError += e[0];

            }
            writer.println(currN + "," + totalMaxError/ (double) numTrials);

        }
        writer.close();
    }

}
