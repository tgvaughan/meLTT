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

package meltt;

import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.tree.TreeParser;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.evolution.likelihood.TreeLikelihood;
import beast.base.spec.evolution.sitemodel.SiteModel;
import beast.base.spec.evolution.substitutionmodel.JukesCantor;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.inference.parameter.SimplexParam;
import beast.base.util.Randomizer;
import beastfx.app.seqgen.SimulatedAlignment;
import org.junit.jupiter.api.Test;

import javax.swing.*;
import java.util.ArrayList;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class LTTLikelihoodTest {

    @Test
    public void test2taxon() {

        // Define substitution model and rate

        SiteModel siteModel = new SiteModel();
        siteModel.initByName("mutationRate", new RealScalarParam<>(1.0, PositiveReal.INSTANCE),
                "substModel", new JukesCantor());

        // Define taxa

        List<Sequence> sequenceList = new ArrayList<>();
        sequenceList.add(new Sequence("t1", "?"));
        sequenceList.add(new Sequence("t2", "?"));
        Alignment dummyAlignment = new Alignment(sequenceList, "nucleotide");

        // Define tree

        TreeParser tree = new TreeParser(dummyAlignment, "(t1:0.05,t2:0.05):0.0;");

        // Simulate sequences

        Randomizer.setSeed(1);
        SimulatedAlignment simulatedAlignment = new SimulatedAlignment();
        simulatedAlignment.initByName("data", dummyAlignment,
                "tree", tree,
                "siteModel", siteModel,
                "sequencelength", 100);

        // Estimate likelihood

        Randomizer.setSeed(42);
        LTT ltt = new LTT();
        ltt.initByName(
                "relCoalIntervals", new SimplexParam(new double[] {1.0}),
                "tMRCA", new RealScalarParam<>(0.05, NonNegativeReal.INSTANCE),
                "alignment", simulatedAlignment);

        LTTLikelihood lttLikelihood = new LTTLikelihood();
        lttLikelihood.initByName(
                "ltt", ltt,
                "nParticles", 100,
                "siteModel", siteModel);
        double lttLikelihoodValue = lttLikelihood.calculateLogP();

        TreeLikelihood treeLikelihood = new TreeLikelihood();
        treeLikelihood.initByName(
                "data", simulatedAlignment,
                "tree", tree,
                "siteModel", siteModel);
        double treeLikelihoodValue = treeLikelihood.calculateLogP();

        System.out.println(" LTTLikelihood: " + lttLikelihoodValue);
        System.out.println("TreeLikelihood: " + treeLikelihoodValue);

        assertEquals(treeLikelihoodValue, lttLikelihoodValue, 1e-10);
    }

    @Test
    public void test3taxon() {

        // Define substitution model and rate

        SiteModel siteModel = new SiteModel();
        siteModel.initByName("mutationRate", new RealScalarParam<>(1.0, PositiveReal.INSTANCE),
                "substModel", new JukesCantor());

        // Define taxa

        List<Sequence> sequenceList = new ArrayList<>();
        sequenceList.add(new Sequence("t1", "?"));
        sequenceList.add(new Sequence("t2", "?"));
        sequenceList.add(new Sequence("t3", "?"));
        Alignment dummyAlignment = new Alignment(sequenceList, "nucleotide");

        // Define tree

        TreeParser tree1 = new TreeParser(dummyAlignment, "(t3:0.05,(t1:0.025,t2:0.025):0.025):0.0;");
        TreeParser tree2 = new TreeParser(dummyAlignment, "(t2:0.05,(t1:0.025,t3:0.025):0.025):0.0;");
        TreeParser tree3 = new TreeParser(dummyAlignment, "(t1:0.05,(t3:0.025,t2:0.025):0.025):0.0;");

        // Simulate sequences

        Randomizer.setSeed(1);
        SimulatedAlignment simulatedAlignment = new SimulatedAlignment();
        simulatedAlignment.initByName("data", dummyAlignment,
                "tree", tree1,
                "siteModel", siteModel,
                "sequencelength", 100);

        // Estimate likelihood

        Randomizer.setSeed(System.currentTimeMillis());
        LTT ltt = new LTT();
        ltt.initByName(
                "relCoalIntervals", new SimplexParam(new double[] {0.5, 0.5}),
                "tMRCA", new RealScalarParam<>(0.05, NonNegativeReal.INSTANCE),
                "alignment", simulatedAlignment);

        LTTLikelihood lttLikelihood = new LTTLikelihood();
        lttLikelihood.initByName(
                "ltt", ltt,
                "nParticles", 100000,
                "siteModel", siteModel);
        double lttLikelihoodValue = lttLikelihood.calculateLogP();

        TreeLikelihood treeLikelihood1 = new TreeLikelihood();
        treeLikelihood1.initByName(
                "data", simulatedAlignment,
                "tree", tree1,
                "siteModel", siteModel);
        TreeLikelihood treeLikelihood2 = new TreeLikelihood();
        treeLikelihood2.initByName(
                "data", simulatedAlignment,
                "tree", tree2,
                "siteModel", siteModel);
        TreeLikelihood treeLikelihood3 = new TreeLikelihood();
        treeLikelihood3.initByName(
                "data", simulatedAlignment,
                "tree", tree3,
                "siteModel", siteModel);

        double[] treeLogLiks = new double[3];
        treeLogLiks[0] = treeLikelihood1.calculateLogP();
        treeLogLiks[1] = treeLikelihood2.calculateLogP();
        treeLogLiks[2] = treeLikelihood3.calculateLogP();

        double maxTreeLogLik = Math.max(treeLogLiks[0],
                Math.max(treeLogLiks[1], treeLogLiks[2]));
        double logAvTreeLik = Math.log(
                (Math.exp(treeLogLiks[0]-maxTreeLogLik)+
                Math.exp(treeLogLiks[1]-maxTreeLogLik) +
                Math.exp(treeLogLiks[2]-maxTreeLogLik))/3.0) + maxTreeLogLik;


        System.out.println("         LTTLikelihood: " + lttLikelihoodValue);
        System.out.println("Average TreeLikelihood: " + logAvTreeLik);

        assertEquals(logAvTreeLik, lttLikelihoodValue, 0.05);
    }
}
