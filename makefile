COMPILER = g++-9

INCL_DIR =  -I/usr/local/include -I/usr/local/include/gsl/ -I.

HEADERS =
LIBRARY = -L/usr/local/lib/ -L/usr/lib/ -L. -lgsl -lstdc++ -lgslcblas -lm -lgfortran -mfpmath=sse -msse3 -funroll-loops 
SOURCES =
CONCERTLIBDIR =

OBJECTS =   MICMain.o  MICParaMgr.o AlignmentParser.o Seq.o AlignData.o \
			MICData.o JointMat.o PseudoCount.o Inversion.o Aggregation.o SanderAgg.o \
			glasso_psicov.o DiscrData.o CntnsData.o SanderPipeline.o ExpressionParser.o \
			DataPartition.o ExpEqlPipe.o CovMat.o IntegratedAgg.o LatentMap.o 
CFLAGS = -g -m64  -fexceptions -DNDEBUG 


INVERSEOBJECTS= InvMat.o InvMain.cc



COVOBJECTS = ComputeCov.o  MICParaMgr.o AlignmentParser.o Seq.o AlignData.o \
			MICData.o JointMat.o PseudoCount.o Inversion.o Aggregation.o SanderAgg.o \
			glasso_psicov.o DiscrData.o CntnsData.o SanderPipeline.o ExpressionParser.o \
			DataPartition.o ExpEqlPipe.o CovMat.o IntegratedAgg.o LatentMap.o 


CFLAGS += -I/usr/local/include -I. 
#-I/usr/include

all: MICDirect Inverse Covariance

clean:
	rm -f *.o Decompose GraphDiff core

MICDirect: $(OBJECTS) $(SOURCES) $(HEADERS)
	$(COMPILER) -o $@ $(CFLAGS) $(OBJECTS)  $(LIBRARY) -lm

Inverse: $(INVERSEOBJECTS) $(SOURCES) $(HEADERS)
	$(COMPILER) -o $@ $(CFLAGS) $(INVERSEOBJECTS)  $(LIBRARY) -lm

Covariance: $(COVOBJECTS) $(SOURCES) $(HEADERS)
	$(COMPILER) -o $@ $(CFLAGS) $(COVOBJECTS)  $(LIBRARY) -lm


%.o:%.cpp
	$(COMPILER) $(CFLAGS) $(INCL_DIR) -c $<

%.o:%.cc
	$(COMPILER) $(CFLAGS) $(INCL_DIR) -c $<

%.o:../common/%.cc
	$(COMPILER) $(CFLAGS) $(INCL_DIR) -c $<

%.o:../common/%.cpp
	$(COMPILER) $(CFLAGS) $(INCL_DIR) -c $<

%.o:%.f90
	$(COMPILER) $(CFLAGS) $(INCL_DIR) -c $<

