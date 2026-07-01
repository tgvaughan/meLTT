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
import beast.base.evolution.distance.JukesCantorDistance;
import beast.base.util.Binomial;
import beast.base.util.Randomizer;

import java.util.Arrays;

public class Particle {

    double[][] initialState, state;
    double[] working;
    int nPatterns;
    int[] patternWeights;
    int nChars;
    int nTaxa;

    double logPrevForestLik, logWeight;

    double[][] initialDistMatrix, distMatrix, mergeProbs;
    int[] cladeSizes;
    int mergeLineage2, mergeLineage1;

    public Particle(Alignment alignment) {
        nChars = alignment.getMaxStateCount();
        nPatterns = alignment.getPatternCount();
        patternWeights = alignment.getWeights();
        nTaxa = alignment.getTaxonCount();

        initialState = new double[nTaxa][nChars*nPatterns];
        state = new double[nTaxa][nChars*nPatterns];
        working = new double[nChars];

        for (int i=0; i<nTaxa; i++) {
            for (int patIdx=0; patIdx<nPatterns; patIdx++) {
                int code = alignment.getPattern(i, patIdx);
                for (int s : alignment.getDataType().getStatesForCode(code))
                    initialState[i][patIdx*nChars + s] = 1.0;
            }
        }

        JukesCantorDistance distanceFunc = new JukesCantorDistance();
        distanceFunc.setPatterns(alignment);
        distMatrix = new double[nTaxa][nTaxa];
        initialDistMatrix = new double[nTaxa][nTaxa];
        cladeSizes = new int[nTaxa];
        for (int i=1; i<nTaxa; i++) {
            cladeSizes[i]=1;
            for (int j=0; j<i; j++) {
                distMatrix[i][j] = distanceFunc.pairwiseDistance(i,j);
                initialDistMatrix[i][j] = distMatrix[i][j];
            }
        }

        mergeProbs = new double[nTaxa][nTaxa];
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

    public void mergeLineages(int k) {

        computeMergeProbsAndSample(k);

        // Compute the merged lineage likelihoods and store in the space occupied
        // by lineage1
        for (int i = 0; i<state[mergeLineage1].length; i++) {
            state[mergeLineage1][i] *= state[mergeLineage2][i];
        }

        // Compute merged distance matrices and clade sizes
        int n1 = cladeSizes[mergeLineage1];
        int n2 = cladeSizes[mergeLineage2];
        int ntot = n1+n2;
        for (int i=0; i<k; i++) {
            if (i==mergeLineage1 || i==mergeLineage2)
                continue;

            double d1 = i>mergeLineage1 ? distMatrix[i][mergeLineage1] : distMatrix[mergeLineage1][i];
            double d2 = i>mergeLineage2 ? distMatrix[i][mergeLineage2] : distMatrix[mergeLineage2][i];

            if (i>mergeLineage1)
                distMatrix[i][mergeLineage1] = (d1*n1 + d2*n2)/(double)ntot;
            else
                distMatrix[mergeLineage1][i] = (d1*n1 + d2*n2)/(double)ntot;
        }
        cladeSizes[mergeLineage1] += cladeSizes[mergeLineage2];

        if (mergeLineage2<k-1) {
            // Swap pointers such that the lineage to discard is the last position
            double[] tmpState = state[k-1];
            state[k-1] = state[mergeLineage2];
            state[mergeLineage2] = tmpState;

            // Swap shift the distances to/from the original lineage k-1 to its
            // new position at mergeLineage2
            for (int i=0; i<k-1; i++) {
                if (i==mergeLineage2)
                    continue;

                if (i>mergeLineage2)
                    distMatrix[i][mergeLineage2] = distMatrix[k-1][i];
                else
                    distMatrix[mergeLineage2][i] = distMatrix[k-1][i];
            }

            // Update lineage index
            mergeLineage2 = k-1;
        }

    }

    public void computeMergeProbsAndSample(int k) {

        double cumsum = 0.0;
        for (int lineage2=1; lineage2<k; lineage2++) {
            for (int lineage1=0; lineage1<lineage2; lineage1++) {
                mergeProbs[lineage2][lineage1] = Math.exp(-k*distMatrix[lineage2][lineage1]);
                cumsum += mergeProbs[lineage2][lineage1];
            }
        }

        double u = Randomizer.nextDouble()*cumsum;
        for (mergeLineage2=1; mergeLineage2<k; mergeLineage2++) {
            for (mergeLineage1=0; mergeLineage1<mergeLineage2; mergeLineage1++) {
                u -= mergeProbs[mergeLineage2][mergeLineage1];
                if (u < 0)
                    break;
            }
            if (u<0)
                break;
        }

        logWeight = -Math.log(mergeProbs[mergeLineage2][mergeLineage1])
                + Math.log(cumsum) - Binomial.logChoose(k, 2);
    }

    public void computeWeight(int k, double[] frequencies) {
        double logForestLik = 0.0;
        for (int lineageIdx=0; lineageIdx<k; lineageIdx++) {
            for (int patIdx=0; patIdx<nPatterns; patIdx++) {
                double siteWeight = 0.0;
                for (int charIdx=0; charIdx<nChars; charIdx++) {
                    siteWeight += state[lineageIdx][patIdx*nChars + charIdx]*frequencies[charIdx];
                }
                logForestLik += Math.log(siteWeight)*patternWeights[patIdx];
            }
        }
        logWeight += logForestLik - logPrevForestLik;
        logPrevForestLik = logForestLik;
    }

    public void reset() {
        for (int i=0; i<nTaxa; i++) {
            System.arraycopy(initialState[i], 0, state[i],
                    0, initialState[i].length);
            System.arraycopy(initialDistMatrix[i], 0, distMatrix[i],
                    0, initialDistMatrix[i].length);
        }
        Arrays.fill(cladeSizes, 1);

        logPrevForestLik = 0.0;
        logWeight = 0.0;
    }

    public void assignFrom(Particle other, int k) {
        for (int lineageIdx=0; lineageIdx<k; lineageIdx++) {
            System.arraycopy(other.state[lineageIdx], 0,
                    state[lineageIdx], 0, state[lineageIdx].length);
            System.arraycopy(other.distMatrix[lineageIdx], 0,
                    distMatrix[lineageIdx], 0, distMatrix[lineageIdx].length);
        }
        System.arraycopy(other.cladeSizes, 0, cladeSizes, 0, cladeSizes.length);

        logWeight = other.logWeight;
        logPrevForestLik = other.logPrevForestLik;
    }

    public void printMatrix(double[][] mat, int k, String desc) {
        System.out.println(desc + " matrix with k=" + k + ":");
        System.out.print("     ");
        for (int j=0; j<k; j++)
            System.out.format("%-6d ", j);
        System.out.println();
        for (int i=0; i<k; i++) {
            System.out.format("%-2d [ ", i);
            for (int j=0; j<k; j++) {
                System.out.format("%-6.3g ", mat[i][j]);
            }
            System.out.println("]");
        }
    }
}
