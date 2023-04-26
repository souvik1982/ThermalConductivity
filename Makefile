ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS  = $(shell root-config --libs)

all: CorrelatedDataFitter.o ThroughPlaneAnalysis_ReadRaw ThroughPlaneAnalysis_ProcessDatacard ThroughPlaneAnalysis_Metaslope InPlaneAnalysis_ReadRaw InPlaneAnalysis_ProcessDatacard InPlaneAnalysis_Metaslope
clean:
	rm CorrelatedDataFitter.o ThroughPlaneAnalysis_ReadRaw ThroughPlaneAnalysis_ProcessDatacard ThroughPlaneAnalysis_Metaslope InPlaneAnalysis_ReadRaw InPlaneAnalysis_ProcessDatacard InPlaneAnalysis_Metaslope

CorrelatedDataFitter.o: CorrelatedDataFitter.cc
	g++ -c CorrelatedDataFitter.cc -o CorrelatedDataFitter.o $(ROOTFLAGS) -mmacosx-version-min=12.6

ThroughPlaneAnalysis_ReadRaw: ThroughPlaneAnalysis_ReadRaw.cc
	g++ -O2 ThroughPlaneAnalysis_ReadRaw.cc -o ThroughPlaneAnalysis_ReadRaw $(ROOTFLAGS) $(ROOTLIBS) -mmacosx-version-min=12.6

ThroughPlaneAnalysis_ProcessDatacard: ThroughPlaneAnalysis_ProcessDatacard.cc
	g++ -O2 ThroughPlaneAnalysis_ProcessDatacard.cc -o ThroughPlaneAnalysis_ProcessDatacard $(ROOTFLAGS) $(ROOTLIBS) -mmacosx-version-min=12.6

ThroughPlaneAnalysis_Metaslope: ThroughPlaneAnalysis_Metaslope.cc
	g++ -O2 ThroughPlaneAnalysis_Metaslope.cc -o ThroughPlaneAnalysis_Metaslope $(ROOTFLAGS) $(ROOTLIBS) -mmacosx-version-min=12.6

InPlaneAnalysis_ReadRaw: InPlaneAnalysis_ReadRaw.cc
	g++ -O2 InPlaneAnalysis_ReadRaw.cc -o InPlaneAnalysis_ReadRaw $(ROOTFLAGS) $(ROOTLIBS) -mmacosx-version-min=12.6

InPlaneAnalysis_ProcessDatacard: InPlaneAnalysis_ProcessDatacard.cc
	g++ -O2 InPlaneAnalysis_ProcessDatacard.cc -o InPlaneAnalysis_ProcessDatacard $(ROOTFLAGS) $(ROOTLIBS) -mmacosx-version-min=12.6

InPlaneAnalysis_Metaslope: InPlaneAnalysis_Metaslope.cc
	g++ -O2 InPlaneAnalysis_Metaslope.cc -o InPlaneAnalysis_Metaslope $(ROOTFLAGS) $(ROOTLIBS) -mmacosx-version-min=12.6
