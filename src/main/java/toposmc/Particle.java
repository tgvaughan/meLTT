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

import beast.base.evolution.alignment.Alignment;

public class Particle {

    double[][] initialState, state;
    double[] working;
    int nPatterns;
    int[] patternWeights;
    int nChars;

    double logForestLik, logPrevForestLik, logWeight;

    public Particle(Alignment alignment) {
        nChars = alignment.getMaxStateCount();
        nPatterns = alignment.getPatternCount();
        patternWeights = alignment.getWeights();

        initialState = new double[alignment.getTaxonCount()][nChars*nPatterns];
        state = new double[alignment.getTaxonCount()][nChars*nPatterns];
        working = new double[nChars];

        for (int i=0; i<alignment.getTaxonCount(); i++) {
            for (int patIdx=0; patIdx<nPatterns; patIdx++) {
                int code = alignment.getPattern(i, patIdx);
                for (int s : alignment.getDataType().getStatesForCode(code))
                    initialState[i][patIdx*nChars + s] = 1.0;
            }
        }
    }

    public void propagateLineages(int k, double[] tMatrix) {
        for (int i=0; i<k; i++) {
            for (int j=0; j<nPatterns; j++) {
                for (int m=0; m<nChars; m++) {
                    working[m] = 0.0;
                    for (int n=0; n<nChars; n++) {
                        working[m] += tMatrix[m*nChars+n]*state[i][j*nChars + n];
                    }
                }
                System.arraycopy(working, 0, state[i], j*nChars, nChars);
            }
        }
    }

    public void mergeLineages(int k, int lineage1, int lineage2) {

        // Enforce lineage1<lineage2
        int tmp;
        if (lineage1>lineage2) {
            tmp = lineage1;
            lineage1 = lineage2;
            lineage2 = tmp;
        }

        // Compute the merged lineage likelihoods and store in the space occupied
        // by lineage1
        for (int i=0; i<state[lineage1].length; i++) {
            state[lineage1][i] *= state[lineage2][i];
        }

        // Swap pointers such that the lineage to discard is the last position
        if (lineage2<k-1) {
            double[] tmpState = state[k-1];
            state[k-1] = state[lineage2];
            state[lineage2] = tmpState;
        }
    }

    public void computeWeight(int k, double[] frequencies) {
        logForestLik = 0.0;
        for (int lineageIdx=0; lineageIdx<k; lineageIdx++) {
            for (int patIdx=0; patIdx<nPatterns; patIdx++) {
                double siteWeight = 0.0;
                for (int charIdx=0; charIdx<nChars; charIdx++) {
                    siteWeight += state[lineageIdx][patIdx*nChars + charIdx]*frequencies[charIdx];
                }
                logForestLik += Math.log(siteWeight)*patternWeights[patIdx];
            }
        }
        logWeight = logForestLik - logPrevForestLik;
        logPrevForestLik = logForestLik;
    }

    public void reset() {
        for (int i=0; i<initialState.length; i++)
            System.arraycopy(initialState[i], 0, state[i], 0, initialState[i].length);

        logForestLik = 0.0;
        logPrevForestLik = 0.0;
        logWeight = 0.0;
    }

    public void assignFrom(Particle other, int k) {
        for (int lineageIdx=0; lineageIdx<k; lineageIdx++) {
            System.arraycopy(other.state[lineageIdx], 0,
                    state[lineageIdx], 0, state[lineageIdx].length);
        }
        logWeight = other.logWeight;
        logForestLik = other.logForestLik;
        logPrevForestLik = other.logPrevForestLik;
    }
}
