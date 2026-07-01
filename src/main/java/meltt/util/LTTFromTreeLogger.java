/*
 * Copyright (c) 2026 ETH Zürich
 *
 * This file is part of meLTT.
 *
 * meLTT is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or (at
 * your option) any later version.
 *
 * meLTT is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with meLTT. If not, see <https://www.gnu.org/licenses/>.
 */

package meltt.util;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.CalculationNode;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

@Description("Logs LTT corresponding to given tree.")
public class LTTFromTreeLogger extends CalculationNode implements Loggable {

    public Input<Tree> treeInput = new Input<>("tree",
            "Tree from which to log LTT.", Input.Validate.REQUIRED);

    int[] k;
    double[] t;

    List<Double> tipAges;
    List<Integer> tipAgeMultiplicities;

    Tree tree;

    @Override
    public void initAndValidate() {
        tree = treeInput.get();

        List<Double> allTipAges = tree.getExternalNodes().stream().map(Node::getHeight).sorted().toList();

        tipAges = new ArrayList<>();
        tipAgeMultiplicities = new ArrayList<>();
        double prevAge = Double.NEGATIVE_INFINITY;
        for (double age : allTipAges) {
            if (age>prevAge) {
                tipAges.add(age);
                tipAgeMultiplicities.add(1);
                prevAge = age;
            } else {
                tipAgeMultiplicities.set(tipAgeMultiplicities.size()-1,
                        tipAgeMultiplicities.getLast()+1);
            }

        }

        k = new int[tipAges.size() + tree.getInternalNodeCount()];
        t = new double[tipAges.size() + tree.getInternalNodeCount()];
    }

    public void updateLTT() {
        List<Double> intAges = tree.getInternalNodes().stream().map(Node::getHeight).sorted().toList();

        int thisk = 0;
        int sampleIdx = 0;
        int intNodeIdx = 0;
        int idx=0;

        for (double internalAge : intAges) {
            while (sampleIdx<tipAges.size() && tipAges.get(sampleIdx)<intAges.getLast()) {
                thisk += tipAgeMultiplicities.get(sampleIdx);
                k[idx] = thisk;
                t[idx] = tipAges.get(sampleIdx);
                idx += 1;
                sampleIdx += 1;
            }

            thisk -= 1;
            k[idx] = thisk;
            t[idx] = internalAge;
            idx += 1;
        }

    }

    @Override
    public void init(PrintStream out) {

        for (int i=0; i<k.length; i++) {
            out.format("k%d\t", i);
        }
        for (int i=0; i<k.length; i++) {
            out.format("t%d\t", i);
        }
    }

    @Override
    public void log(long sample, PrintStream out) {
        updateLTT();
        for (int i=0; i<k.length; i++) {
            out.format("%d\t", k[i]);
        }
        for (int i=0; i<k.length; i++) {
            out.format("%g\t", t[i]);
        }
    }

    @Override
    public void close(PrintStream out) {
    }
}
