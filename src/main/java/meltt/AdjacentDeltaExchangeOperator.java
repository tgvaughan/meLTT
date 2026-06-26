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
import beast.base.inference.Operator;
import beast.base.inference.operator.kernel.KernelDistribution;
import beast.base.spec.inference.parameter.SimplexParam;
import beast.base.util.Randomizer;

import java.awt.image.Kernel;

public class AdjacentDeltaExchangeOperator extends Operator {

    public Input<SimplexParam> paramInput = new Input<>("parameter",
            "Simplex parameter to operate on",
            Input.Validate.REQUIRED);

//    public Input<Double> deltaInput = new Input<>("delta",
//            "Maximum offset to use", 0.5);

//    public final Input<Boolean> autoOptimizeiInput =
//            new Input<>("autoOptimize", "Enable auto-tuning of delta input", true);

    public Input<KernelDistribution> kernelDistributionInput = new Input<>("kernel",
            "Kernel distribution from which updates offsets are drawn.",
            KernelDistribution.newDefaultKernelDistribution());

    SimplexParam param;
    double delta;
    KernelDistribution kernelDistribution;
    boolean autoOptimize;

    @Override
    public void initAndValidate() {
        param = paramInput.get();
//        delta = deltaInput.get();
        kernelDistribution = kernelDistributionInput.get();

//        autoOptimize = autoOptimizeiInput.get();
    }

    @Override
    public double proposal() {

        if (param.size() < 2)
            return Double.NEGATIVE_INFINITY;

        int idx = Randomizer.nextInt(param.size() - 1);

        double sum = param.get(idx) + param.get(idx+1);
        double newVal1 = Randomizer.nextDouble()*sum;
        double newVal2 = sum - newVal1;

//        double change = delta * (Randomizer.nextDouble() - 0.5);
//        double change = kernelDistribution.getRandomDelta(0, delta);

//        double newVal1 = param.get(idx) + change;
//        double newVal2 = param.get(idx + 1) - change;

        if (param.isValid(newVal1) && param.isValid(newVal2)) {
            param.set(idx, newVal1);
            param.set(idx + 1, newVal2);
        } else {
            return Double.NEGATIVE_INFINITY;
        }

        return 0;
    }

//    @Override
//    public void optimize(double logAlpha) {
//        if (autoOptimize) {
//            double _delta = calcDelta(logAlpha);
//            _delta += Math.log(delta);
//            delta = Math.exp(_delta);
//        }
//    }
}
