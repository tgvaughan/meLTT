/*
 * Copyright (c) 2026 ETH Zürich
 *
 * This file is part of TopoSMC.
 *
 * TopoSMC is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or (at
 * your option) any later version.
 *
 * TopoSMC is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with TopoSMC. If not, see <https://www.gnu.org/licenses/>.
 */

package toposmc;

import beast.base.core.Input;
import beast.base.evolution.tree.coalescent.PopulationFunction;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.util.Binomial;

import java.util.List;
import java.util.Random;

public class CoalescentLTTPrior extends Distribution {

    public Input<LTT> lttInput = new Input<>(
            "ltt",
            "Lineage through time plot.",
            Input.Validate.REQUIRED);
    public Input<PopulationFunction> popFuncInput = new Input<>(
            "populationFunction",
            "Effective population size function.",
            Input.Validate.REQUIRED);

    LTT ltt;
    PopulationFunction popFunc;

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        ltt = lttInput.get();
        popFunc = popFuncInput.get();
    }

    @Override
    public double calculateLogP() {
        logP = 0.0;

        for (int idx=1; idx<ltt.t.length; idx++) {
            logP += -Binomial.choose(ltt.k[idx-1],2)*popFunc.getIntegral(ltt.t[idx-1], ltt.t[idx]);
            if (ltt.k[idx]>ltt.k[idx-1])
                logP +=  Binomial.logChoose(ltt.k[idx-1], 2) - Math.log(popFunc.getPopSize(ltt.t[idx]));
        }

        return logP;
    }

    @Override
    public List<String> getArguments() {
        return List.of();
    }

    @Override
    public List<String> getConditions() {
        return List.of();
    }

    @Override
    public void sample(State state, Random random) {

    }
}
