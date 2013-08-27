FC	= g95
FFLAGS	= -O3 -freal-loops
OBJS	= cosmo.o linearpk.o sigma.o spline.o hunt.o qromb.o chebyshev.o growth.o dverk.o mf_jenkins.o
.SUFFIXES: .f90
.f90.o: 
	$(FC) $(FFLAGS) -c $<
default:  sample compute_mf_jenkins compute_mf_jenkins2
compute_mf_jenkins: $(OBJS) compute_mf_jenkins.o
	$(FC) $(FFLAGS) -o $@ $(OBJS) $@.o $(LDFLAGS)
compute_mf_jenkins2: $(OBJS) compute_mf_jenkins2.o
	$(FC) $(FFLAGS) -o $@ $(OBJS) $@.o $(LDFLAGS)
sample: $(OBJS) sample.o
	$(FC) $(FFLAGS) -o $@ $(OBJS) $@.o $(LDFLAGS)
clean:
	-rm -f *.o *.mod
tidy: clean
	-rm -f sample
	-rm -f compute_mf_jenkins
	-rm -f compute_mf_jenkins2
	-rm -f Rh_dndlnRh.txt
	-rm -f Mh_dndlnMh.txt
	-rm -f Mhom0_dndlnMh.txt
