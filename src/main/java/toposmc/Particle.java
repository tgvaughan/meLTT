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
    int nPatterns;
    int nChars;

    double logWeight;

    public Particle(Alignment alignment) {
        nChars = alignment.getMaxStateCount();
        nPatterns = alignment.getPatternCount();

        initialState = new double[alignment.getTaxonCount()][nChars*nPatterns];
        state = new double[alignment.getTaxonCount()][nChars*nPatterns];

        for (int i=0; i<alignment.getTaxonCount(); i++) {
            for (int j=0; j<alignment.getPatternCount(); j++) {
                initialState[i][j*nChars + alignment.getPattern(i, j)] = 1.0;
            }
        }
    }

    public void reset() {
        for (int i=0; i<initialState.length; i++)
            System.arraycopy(initialState[i], 0, state[i], 0, initialState[i].length);

        logWeight = 0.0;
    }
}
