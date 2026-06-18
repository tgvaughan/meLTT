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
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.spec.evolution.sitemodel.SiteModel;

import java.util.List;
import java.util.Random;

public class LTTLikelihood extends Distribution {

    public Input<LTT> lttInput = new Input<>("ltt",
            "Lineage through time plot (relative interval )", Input.Validate.REQUIRED);

    public Input<Integer> nParticlesInput = new Input<>("nParticles",
            "Number of particles to use in the SMC.", 100);

    public Input<SiteModel> siteModelInput = new Input<>("siteModel",
            "Site model describing evolution of sites.",
            Input.Validate.REQUIRED);

    SiteModel siteModel;
    LTT ltt;

    Particle[] particles;
    int nParticles;

    double[] tpMatrix;

    @Override
    public void initAndValidate() {
        siteModel = siteModelInput.get();
        ltt = lttInput.get();

        int nStates = ltt.alignment.getMaxStateCount();
        tpMatrix = new double[nStates*nStates];

        nParticles = nParticlesInput.get();
        particles = new Particle[nParticles];
        for (int i=0; i<nParticles; i++)
            particles[i] = new Particle(ltt.alignment);

        System.out.println("LTTLikelihood initialised.");
    }

    @Override
    public double calculateLogP() {
        logP = 0.0;

        if (!ltt.isValid()) {
            logP = Double.NEGATIVE_INFINITY;
            return logP;
        }

        // TODO Define particle state
        for (Particle p : particles)
            p.reset();

        for (int idx=1; idx<ltt.t.length; idx++) {

            siteModel.getSubstitutionModel().getTransitionProbabilities(null,
                    ltt.t[idx-1], ltt.t[idx], 1.0, tpMatrix);


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

    public static void main(String[] args) {
        System.out.println("Testing.");
    }
}
