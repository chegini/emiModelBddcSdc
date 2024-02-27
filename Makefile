include /data/numerik/software/KaskadeDependencies/Kaskade7.5Dependencies-10.2/installed/Makefile.Local
include /data/numerik/people/fchegini/kaskade7/Makefile.Local
include /data/numerik/people/fchegini/kaskade7/Makefile.Rules



INCLUDE = $(DUNEINC) $(UGINC) $(BOOSTINC) $(AMIRAINC) $(KASKADEINC) $(ITSOLINC) $(HYPREINC) $(TAUCSINC) $(TRILINOSINC) $(UMFPACKINC) $(MUMPSINC) $(PSURFACEINC)
FLAGS = $(WFLAGS) $(OPTFLAGS)
 
emiModel: EMI_model.o Makefile $(KASKADE7)/libs/libkaskade.a
	$(CXX) $(FLAGS) -std=c++17 $(HAVE_UG) $< $(KASKADE7)/libs/umfpack_solve.o $(KASKADE7)/libs/mumps_solve.o \
	$(KASKADELIB) \
 $(DUNELIB) \
 $(UGLIB) \
 $(BOOSTLIB) \
 $(UMFPACKLIB) \
 $(MUMPSLIB) \
 $(SUPERLULIB) \
 $(ITSOLLIB) \
 $(HYPRELIB) \
 $(TAUCSLIB) \
 $(TRILINOSLIB) \
 $(PARDISOLIB) \
 $(AMIRALIB) \
 $(BLASLIB) $(FTNLIB) $(NUMALIB) $(LINKFLAGS) -o $@

# depend:
# 	cp Makefile.gen Makefile; ../../tools/gccmakedep -f Makefile $(INCLUDE)  $(MAKEDEPENDFLAGS) EMI_model.cpp / 
# 	   ../../tools/remove_install_deps Makefile

depend:
	cp Makefile.gen Makefile; /data/numerik/people/fchegini/kaskade7/tools/gccmakedep -f Makefile $(INCLUDE)  $(MAKEDEPENDFLAGS) EMI_model.cpp / 
	   /data/numerik/people/fchegini/kaskade7/tools/remove_install_deps Makefile

clean:
	rm -f gccerr.txt *.o emiModel EMI_model-*.vtu sol-*.vtu cor-*.vtu *.txt *.m ./output/*.vtu *.vtu ./matlab_dir/*

# DO NOT DELETE

%.o: %.cpp
	$(CXX) $(FLAGS) -std=c++17 $(HAVE_UG) $(INCLUDE) $(CLANGINC) $< -c -o $@ 2> gccerr.txt
