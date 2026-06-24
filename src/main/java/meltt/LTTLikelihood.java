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
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.spec.evolution.sitemodel.SiteModel;
import beast.base.util.Randomizer;

import java.util.List;
import java.util.Random;

public class LTTLikelihood extends Distribution {

    public Input<LTT> lttInput = new Input<>("ltt",
            "Lineage through time plot (relative interval )", Input.Validate.REQUIRED);

    public Input<Integer> nParticlesInput = new Input<>("nParticles",
            "Number of particles to use in the SMC.", 50);

    public Input<SiteModel> siteModelInput = new Input<>("siteModel",
            "Site model describing evolution of sites.",
            Input.Validate.REQUIRED);

    SiteModel siteModel;
    LTT ltt;

    Particle[] particles, particlesPrime;
    int nParticles;
    double[] logWeights, weightsWorking;

    double[] tMatrix;

    @Override
    public void initAndValidate() {
        siteModel = siteModelInput.get();
        ltt = lttInput.get();

        int nStates = ltt.alignment.getMaxStateCount();
        tMatrix = new double[nStates*nStates];

        nParticles = nParticlesInput.get();
        particles = new Particle[nParticles];
        particlesPrime = new Particle[nParticles];
        for (int i=0; i<nParticles; i++) {
            particles[i] = new Particle(ltt.alignment);
            particlesPrime[i] = new Particle(ltt.alignment);
        }

        logWeights = new double[nParticles];
        weightsWorking = new double[nParticles];

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

        int idx=0;
        while (true) {
            idx += 1;

            siteModel.getSubstitutionModel().getTransitionProbabilities(null,
                    ltt.t[idx], ltt.t[idx-1],
                    siteModel.getRateForCategory(0,null),
                    tMatrix);

            for (Particle p : particles) {
                // Propagate partial likelihoods
                p.propagateLineages(ltt.k[idx-1], tMatrix);

                // Randomly merge lineages
                p.mergeLineages(ltt.k[idx-1]);

                // Compute weight
                p.computeWeight(ltt.k[idx], siteModel.getSubstitutionModel().getFrequencies());
            }

            // Estimate likelihood contribution

            double maxLogWeight = Double.NEGATIVE_INFINITY;
            for (int pidx=0; pidx<nParticles; pidx++) {
                logWeights[pidx] = particles[pidx].logWeight;
                maxLogWeight = Math.max(maxLogWeight, logWeights[pidx]);
            }

            int nextParticleToAssign = 0;

            double cumSum = 0.0;
            for (int pidx=0; pidx<nParticles; pidx++) {
                weightsWorking[pidx] = Math.exp(logWeights[pidx] - maxLogWeight);
                cumSum += weightsWorking[pidx];
            }

            if (idx==ltt.t.length-1) {
                // We're done.  Compute final likelihood estimate and quit.
                logP = Math.log(cumSum) + maxLogWeight - Math.log(nParticles);
                return logP;
            }

            double cumSumResiduals = 0.0;
            for (int pidx=0; pidx<nParticles; pidx++) {
                weightsWorking[pidx] *= nParticles / cumSum;

                int nForSure = (int) weightsWorking[pidx];
                weightsWorking[pidx] -= nForSure;
                cumSumResiduals += weightsWorking[pidx];

                for (int i=0; i<nForSure; i++)
                    particlesPrime[nextParticleToAssign++].assignFrom(particles[pidx], ltt.k[idx]);
            }
            for (int pidx=0; pidx<nParticles; pidx++) {
                weightsWorking[pidx] /= cumSumResiduals;
            }

            ReplacementSampler replacementSampler = new ReplacementSampler(weightsWorking);
            for (int i=nextParticleToAssign; i<nParticles; i++)
                particlesPrime[i].assignFrom(particles[replacementSampler.next()], ltt.k[idx]);
            Particle[] tmp = particles;
            particles = particlesPrime;
            particlesPrime = tmp;
        }
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

    @Override
    public boolean isStochastic() {
        return true;
    }
}
