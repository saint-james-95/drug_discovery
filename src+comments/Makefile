#===============================================================================
#         ___________.__            .____________.__  __          
#     ____\_   _____/|__| ____    __| _/   _____/|__|/  |_  ____  
#   _/ __ \|    __)  |  |/    \  / __ |\_____  \ |  \   __\/ __ \ 
#   \  ___/|     \   |  |   |  \/ /_/ |/        \|  ||  | \  ___/ 
#    \___  >___  /   |__|___|  /\____ /_______  /|__||__|  \___  >
#        \/    \/            \/      \/       \/               \/ 
#
#                                                  
#   eFindSite - ligand-binding site prediction from meta-threading
#
#   Computational Systems Biology Group
#   Department of Biological Sciences
#   Center for Computation & Technology
#   Louisiana State University
#   407 Choppin Hall, Baton Rouge, LA 70803, USA
#
#   http://www.brylinski.org
#
#   Report bugs to michal@brylinski.org wfeinstein@lsu.edu
#
#   Copyright 2013 Michal Brylinski, Wei Feinstein
#
#   This file is part of eFindSite.
#
#   eFindSite is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   eFindSite is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with eFindSite. If not, see <http://www.gnu.org/licenses/>.
#
#===============================================================================

EXE = partA partB partC parts_phi

OBJ_CMN = cluster.o cmps.o coords.o distance.o data.o frtmalign_omp.o list.o nwalign.o pocket.o refine.o runsvm.o tanimoto.o target.o template.o walltime.o
OBJ_PTA = efindsite_omp_pA.o
OBJ_PTB = efindsite_omp_pB.o
OBJ_PTC = efindsite_omp_pC.o
PARTDIR = parts

SH = sh

CXX = g++

FC = gfortran

CPPFLAGS = -O0 -g -fopenmp -Wall -Wno-write-strings -I../lib -I. -I/usr/include -D_GLIBCXX_USE_CXX11_ABI=0

FFLAGS = -O0 -g -fopenmp 

OMPLDFLAGS = -lgzstream -lsvm -lz -lgfortran -fopenmp -lm -ldl -L../lib 

SCRLDFLAGS = -lgzstream -lsvm -lz -fopenmp -lm -ldl -L../lib -L. -L/usr/lib

default: $(EXE)

all: $(EXE)

partA: $(OBJ_CMN) $(OBJ_PTA)
	$(CXX) -o $@ $(OBJ_CMN) $(OBJ_PTA) $(OMPLDFLAGS)
	@mv partA ../bin/
$(OBJ_PTA): %.o: $(PARTDIR)/%.C
	$(CXX) $(CPPFLAGS) -c -o $@ $<

partB: $(OBJ_CMN) $(OBJ_PTB)
	$(CXX) -o $@ $(OBJ_CMN) $(OBJ_PTB) $(OMPLDFLAGS)
	@mv partB ../bin/
$(OBJ_PTB): %.o: $(PARTDIR)/%.C
	$(CXX) $(CPPFLAGS) -c -o $@ $<

partC: $(OBJ_CMN) $(OBJ_PTC)
	$(CXX) -o $@ $(OBJ_CMN) $(OBJ_PTC) $(OMPLDFLAGS)
	@mv partC ../bin/
$(OBJ_PTC): %.o: $(PARTDIR)/%.C
	$(CXX) $(CPPFLAGS) -c -o $@ $<
parts_phi:
	$(SH) $(PARTDIR)/parts_phi.shar
	@chmod +x parts_phi
	@mv parts_phi ../bin/

cluster.o: cluster.C
	$(CXX) $(CPPFLAGS) -c -o cluster.o cluster.C

cmps.o: cmps.C
	$(CXX) $(CPPFLAGS) -c -o cmps.o cmps.C

coords.o: coords.C
	$(CXX) $(CPPFLAGS) -c -o coords.o coords.C

distance.o: distance.C
	$(CXX) $(CPPFLAGS) -c -o distance.o distance.C

data.o: data.C
	$(CXX) $(CPPFLAGS) -c -o data.o data.C

efindsite_scr.o: efindsite_scr.C
	$(CXX) $(CPPFLAGS) -c -o efindsite_scr.o efindsite_scr.C

list.o: list.C
	$(CXX) $(CPPFLAGS) -c -o list.o list.C

nwalign.o: nwalign.f
	$(FC) $(FFLAGS) -c -o nwalign.o nwalign.f

pocket.o: pocket.C
	$(CXX) $(CPPFLAGS) -c -o pocket.o pocket.C

refine.o: refine.C
	$(CXX) $(CPPFLAGS) -c -o refine.o refine.C

runsvm.o: runsvm.C
	$(CXX) $(CPPFLAGS) -c -o runsvm.o runsvm.C

target.o: target.C
	$(CXX) $(CPPFLAGS) -c -o target.o target.C

tanimoto.o: tanimoto.C
	$(CXX) $(CPPFLAGS) -c -o tanimoto.o tanimoto.C

template.o: template.C
	$(CXX) $(CPPFLAGS) -c -o template.o template.C

walltime.o: walltime.C
	$(CXX) $(CPPFLAGS) -c -o walltime.o walltime.C

efindsite_omp.o: efindsite_omp.C
	$(CXX) $(CPPFLAGS) -c -o efindsite_omp.o efindsite_omp.C

frtmalign_omp.o: frtmalign_omp.f
	$(FC) $(FFLAGS) -c -o frtmalign_omp.o frtmalign_omp.f

clean:
	@(rm -f ${OBJ_CMN} ${OBJ_PTA} ${OBJ_PTB} ${OBJ_PTC} ../bin/partA ../bin/partB ../bin/partC ../bin/parts_phi)
