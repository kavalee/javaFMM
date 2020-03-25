package com.kavalee;

import java.util.*;

public class CauchyMultiplier {
    double[] s;
    double[] t;
    double[] g;
    int p;
    int N;
    int L;
    NavigableMap<Double,Integer> reverseMapT;
    NavigableMap<Double,Integer> reverseMapS;

    CauchyMultiplier(double[] s, double[] t, double[] f, int p) {
        this.s = s;
        this.t = t;
        this.g = f;
        this.p = p;
        N = g.length;
        L = (int) (Math.log(N)/Math.log(2.0));
        reverseMapT = new TreeMap<>();
        reverseMapS = new TreeMap<>();
        for (int i = 0; i < t.length; i++) {
            reverseMapT.put(t[i],i);
            reverseMapS.put(s[i],i);
        }
    }
    public double[] multiplySlow() {
        int n = g.length;
        double[] f = new double[n];
        for(int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                f[i] += g[j] * cij(i,j);
            }
        }
        return f;
    }
    public double[] multiplyFast3() {
        double[] f = new double[N];
        Stack<Node> toBeProc = new Stack<>();
        Stack<Node> toBeProcNext = new Stack<>();
        toBeProc.add(new Node(0,0));
        toBeProc.add(new Node(1,0));
        toBeProc.add(new Node(0,1));
        toBeProc.add(new Node(1,1));
        for (int l = 1; l <= L; l++) {
            int q = 2 << (l - 1);
            double[][] tMoments = new double[q][p];
            double[] cen = new double[q];
            for (int i = 0; i < q; i++) {
                cen[i] = (double) i / (double) q + 1.0/ (2.0 * q);
            }
            //compute sMoments
            List<Double>[] sCells = sortIntoCells(s,q);
            List<Double>[] tCells = sortIntoCells(t,q);
            double[][] sMoments = computeSMoments(sCells);
            //compute tMoments
            while (!toBeProc.isEmpty()) {
                Node n = toBeProc.pop();
                if(n.isWellSeperated()) {
                    double dc = cen[n.t] - cen[n.s];
                    for (int k = 0; k < p; k++) {
                        double sum = 0;
                        for (int m = 0; m < p; m++) {
                            sum += (Math.pow(dc, -k) *
                                    (Math.pow(dc, -m - 1) * sMoments[n.s][m]))
                                    * nCk(m + k, k);
                        }
                        tMoments[n.t][k] += sum;
                    }
                } else {
                    if (l == L){
                        for (double tv : tCells[n.t]) {
                            int i = reverseMapT.get(tv);
                            double sum = 0;
                            for (double sv : sCells[n.s]) {
                                int j = reverseMapS.get(sv);
                                sum += g[j] * 1.0 / (tv - sv);
                            }
                            f[i] += sum;
                        }
                    } else {
                        toBeProcNext.push(new Node(2 * n.s, 2 * n.t));
                        toBeProcNext.push(new Node(2 * n.s + 1, 2 * n.t));
                        toBeProcNext.push(new Node(2 * n.s, 2 * n.t + 1));
                        toBeProcNext.push(new Node(2 * n.s + 1, 2 * n.t + 1));
                    }
                }
            }

            for (int ti = 0; ti < q; ti++) {
                for (double tv : tCells[ti]) {
                    double sum = 0;
                    int i = reverseMapT.get(tv);
                    for (int k = 0; k < p; k++) {
                        sum += tMoments[ti][k] * Math.pow(cen[ti] - tv, k);
                    }
                    f[i] += sum;
                }
            }
            toBeProc.addAll(toBeProcNext);
            toBeProcNext.clear();
        }
        return f;
    }

    public double[] multiplyFast2() {
        double[] f = new double[N];
        boolean[][] doneWork = new boolean[2][2];
        for (int l = 2; l <= L; l++) {
            int q = 2 << (l - 1);
            //int ps = Math.min(p, q);
            int ps = p;
            boolean[][] doneWorkTemp = new boolean[q][q];
            List<Double>[] sCells = sortIntoCells(s,q);
            List<Double>[] tCells = sortIntoCells(t,q);
            double[][] sMoments = computeSMoments(sCells);
            double[][] tMoments = new double[q][ps];
            for (int si = 0; si < q; si++) {
                for (int ti = 0; ti < q; ti++) {
                    if(doneWork[si/2][ti/2]) {
                        doneWorkTemp[si][ti] = true;
                    } else {
                        if (Math.abs(ti-si) > 1) {
                            double dc = (ti - si)/(double) q;
                            for (int k = 0; k < ps; k++) {
                                double sum = 0;
                                for (int m = 0; m < ps; m++) {
                                    sum += (Math.pow(dc,-k) *  (Math.pow(dc, -m - 1) * sMoments[si][m])) * nCk(m + k, k);
                                }
                                tMoments[ti][k] += sum;
                            }
                            doneWorkTemp[si][ti] = true;
                        } else {
                            doneWorkTemp[si][ti] = false;
                        }
                    }
                }
            }
            doneWork = doneWorkTemp;
            if (l == L) {
                for (int si = 0; si < q; si++) {
                    for (int ti = 0; ti < q; ti++) {
                        if (Math.abs(ti - si) <= 1) {
                            for (double tv : tCells[ti]) {
                                int i = reverseMapT.get(tv);
                                double sum = 0;
                                for (double sv : sCells[si]) {
                                    int j = reverseMapS.get(sv);
                                    sum += g[j] * 1.0/(tv-sv);
                                }
                                f[i] += sum;
                            }
                        }
                    }
                }
            }
            for (int ti = 0; ti < q; ti++) {
                double tc = (double) ti / (double) q + 1.0 / (2.0 * q);
                for (double tv : tCells[ti]) {
                    double sum = 0;
                    int i = reverseMapT.get(tv);
                    for (int k = 0; k < ps; k++) {
                        sum += tMoments[ti][k] * Math.pow(tc - tv, k);
                    }
                    f[i] += sum;
                }
            }
        }
        return f;
    }

    public double[][] computeSMoments(List<Double>[] sCells) {
        int n = sCells.length;
        //int ps = Math.min(n,p);
        int ps = p;
        double[][] sMoments = new double[n][ps];
        for (int si = 0; si < n; si++) {
            double center = (double) si / (double) n + 1f/ (double) (2*n);
            for (int k = 0; k < ps; k++) {
                float sum = 0;
                for (double sv : sCells[si]) {
                    int i = reverseMapS.get(sv);
                    sum += g[i] * Math.pow(sv - center, k);
                }
                sMoments[si][k] = sum;
            }
        }
        return sMoments;
    }
        class Node {
            int s;
            int t;
            Node(int s, int t) {
                this.s = s;
                this.t = t;
            }

            public boolean isWellSeperated() {
                return  Math.abs(s - t) > 1;
            }
            @Override
            public String toString(){
                return "(" + s + ", " + t + ")";
            }
            @Override
            public boolean equals(Object o) {
                Node n = (Node) o;
                return (n.s == s && n.t == t);
            }
        }
    public ArrayList[] sortIntoCells(double[] s, int n) {
        ArrayList<Double>[] cells = new ArrayList[n];

        for(int i = 0; i < n; i++) {
            cells[i] = new ArrayList();
        }
        for (double x : s) {
            cells[(int)(x * n)].add(x);
        }
        return cells;
    }
    public void printMatrix(boolean[][] m) {
        int q = m.length;
        System.out.println("=====" + q + "=====");
        for (int s = 0; s < q; s++) {
            for (int t = 0; t < q; t++) {
                String c = (m[s][t]) ? "w" :"n";
                System.out.print(c);
            }
            System.out.println();
        }
        System.out.println("=============");
    }
    public void printNodes(List<Node> nodes,int q) {
        System.out.println("=====" + q + "=====");
        for (int si = 0; si < q; si++) {
            for (int ti = 0; ti < q; ti++) {
                String s = nodes.contains(new Node(si,ti))? "." : " ";
                System.out.print(s);
            }
            System.out.println();
        }
        System.out.println("=============");

    }


    public double cij(int i, int j){
        return 1/(t[i] - s[j]);
    }
    public static int factorial(int n) {
        if (n < 0) {
            throw new Error("factorial not defined");
        }
        if(n == 0) {
            return 1;
        }
        int prod = 1;
        for (int i = 2; i <= n; i++) {
            prod *= i;
        }
        return prod;
    }
    public static double nCk(int n, int k) { ;
        double prod = 1;
        for (int i = 0; i < k; i++) {
            prod *= (double) (n- i) / (double) (k - i);
        }
        return prod;
    }
}
