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

import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.core.Loggable;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.tree.TraitSet;
import beast.base.inference.CalculationNode;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.inference.parameter.SimplexParam;
import beast.base.spec.type.RealScalar;
import beast.base.spec.type.Simplex;
import beast.base.util.Randomizer;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Instances of this class represent an LTT function.
 *
 * The LTT function itself is computed from the tMRCA, the set of
 * inter-coalescent interval durations relative to the tMRCA, and
 * (optionally) a traitset specifying tip sampling times.
 *
 * On initialization the relCoalInterval input parameter
 * is modified such that its dimension matches the number of
 * intervals requires by the number of taxa.
 *
 * The computed LTT function is exposed via the public fields t[] and k[].
 *
 * Tip age sampling is not (yet?) supported.
 */
public class LTT extends CalculationNode implements Loggable {

    public Input<Simplex> relCoalIntervalsInput = new Input<>("relCoalIntervals",
            "Intervals between coalescent events, relative to tMRCA",
            Input.Validate.REQUIRED);

    public Input<RealScalar<? extends NonNegativeReal>> tMRCAInput = new Input<>("tMRCA",
            "Age of MRCA relative to most recent sample.",
            Input.Validate.REQUIRED);

    public Input<Alignment> alignmentInput = new Input<>("alignment",
            "Taxon set", Input.Validate.REQUIRED);

    public Input<TraitSet> traitSetInput = new Input<>("traitSet",
            "Trait set defining tip dates.");

    Simplex relCoalIntervals;
    RealScalar<?> tMRCA;

    Alignment alignment;
    TraitSet traitSet;

    int[] k;
    double[] t;

    List<Double> tipAges;
    List<Integer> tipAgeMultiplicities;

    boolean validLTT;

    @Override
    public void initAndValidate() {
        relCoalIntervals = relCoalIntervalsInput.get();
        tMRCA = tMRCAInput.get();

        alignment = alignmentInput.get();
        traitSet = traitSetInput.get();

        // For now assume tip ages are fixed
        tipAges = new ArrayList<>();
        tipAgeMultiplicities = new ArrayList<>();
        if (traitSet != null) {
            if (!traitSet.isDateTrait())
                throw new IllegalArgumentException("traitSet provided to LTT must specify a date trait.");

            for (int i=0; i<alignment.getTaxonCount(); i++) {
                tipAges.add(traitSet.getValue(alignment.getTaxaNames().get(i)));
                tipAgeMultiplicities.add(1);
            }
            tipAges.sort(Double::compareTo);
            int idx=0;
            while (idx<tipAges.size()-1) {
                if (tipAges.get(idx + 1).equals(tipAges.get(idx))) {
                    tipAgeMultiplicities.set(idx, tipAgeMultiplicities.get(idx) + tipAgeMultiplicities.get(idx + 1));
                    tipAgeMultiplicities.remove(idx + 1);
                    tipAges.remove(idx + 1);
                } else {
                    idx += 1;
                }
            }
        } else {
            tipAges.add(0.0);
            tipAgeMultiplicities.add(alignment.getTaxonCount());
        }

        if (relCoalIntervals.size() != alignment.getTaxonCount()-1) {

            if (relCoalIntervals instanceof SimplexParam relCoalIntervalsParam &&
                    tMRCA instanceof RealScalarParam<?> tMRCAParam) {

                Log.debug.println("Attempting to initialise relCoalIntervals and tMRCA.");
                initLTTParams(relCoalIntervalsParam, tMRCAParam);
            } else
                throw new IllegalArgumentException("relCoalIntervals has wrong dimension and can't initialise");
        }

        int nEvents = tipAges.size() + alignment.getTaxonCount() - 1;
        t = new double[nEvents];
        k = new int[nEvents];

        updateLTT();
    }

    public void updateLTT() {

        int nextSampleIdx = 0;
        int nextCoalIdx = 0;
        int nextEventIdx = 0;
        int thisk = 0;

        double thisCoalTime = 0.0;

        while (nextCoalIdx < relCoalIntervals.size()) {
            thisCoalTime += tMRCA.get()*relCoalIntervals.get(nextCoalIdx);

            while (nextSampleIdx<tipAges.size() && tipAges.get(nextSampleIdx)<thisCoalTime) {
                t[nextEventIdx] = tipAges.get(nextSampleIdx);
                thisk += tipAgeMultiplicities.get(nextSampleIdx);
                k[nextEventIdx] = thisk;

                nextSampleIdx += 1;
                nextEventIdx += 1;
            }

            t[nextEventIdx] = thisCoalTime;
            thisk -= 1;

            if (thisk<1) {
                validLTT = false;
                return;
            }

            k[nextEventIdx] = thisk;
            nextCoalIdx += 1;
            nextEventIdx += 1;
        }

        validLTT = true;
    }

    /**
     * @return true only when k>0 for all intervals
     */
    public boolean isValid() {
        return validLTT;
    }

    public void initLTTParams(SimplexParam relCoalIntervalsParam, RealScalarParam<?> tMRCAParam) {

        List<Double> sampleTimes = new ArrayList<>(tipAges);
        List<Integer> sampleTimeMultiplicities = new ArrayList<>(tipAgeMultiplicities);

        Collections.reverse(sampleTimes); // Youngest samples at end
        Collections.reverse(sampleTimeMultiplicities);

        List<Double> coalTimes = new ArrayList<>();

        double t = 0;
        int k = 0;
        while (!sampleTimes.isEmpty() || k>1) {

            double coalProp = k*(k-1)/2.0;

            double dt = coalProp>0 ? Randomizer.nextExponential(coalProp) : Double.POSITIVE_INFINITY;

            if (!sampleTimes.isEmpty() && t+dt > sampleTimes.getLast()) {
                t = sampleTimes.removeLast();
                k += sampleTimeMultiplicities.removeLast();
                continue;
            }

            t += dt;
            k -= 1;

            coalTimes.add(t);
        }

        tMRCAParam.set(coalTimes.getLast());

        relCoalIntervalsParam.setDimension(coalTimes.size());
        for (int i=0; i<coalTimes.size(); i++) {
            double delta = coalTimes.get(i);
            if (i>0) delta -= coalTimes.get(i-1);
            relCoalIntervalsParam.set(i, delta/coalTimes.getLast());
        }

    }

    @Override
    protected boolean requiresRecalculation() {
        updateLTT();
        return true;
    }

    @Override
    protected void restore() {
        super.restore();
        updateLTT();
    }

    /*
     * Loggable implementation
     */

    @Override
    public void init(PrintStream out) {

        updateLTT();
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
