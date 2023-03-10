ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS  = $(shell root-config --libs)

all: ThroughPlaneAnalysis_ReadRaw ThroughPlaneAnalysis_ProcessDatacard ThroughPlaneAnalysis_Metaslope

ThroughPlaneAnalysis_ReadRaw: ThroughPlaneAnalysis_ReadRaw.cc
	g++ -O2 ThroughPlaneAnalysis_ReadRaw.cc -o ThroughPlaneAnalysis_ReadRaw $(ROOTFLAGS) $(ROOTLIBS) -mmacosx-version-min=12.6

ThroughPlaneAnalysis_ProcessDatacard: ThroughPlaneAnalysis_ProcessDatacard.cc
	g++ -O2 ThroughPlaneAnalysis_ProcessDatacard.cc -o ThroughPlaneAnalysis_ProcessDatacard $(ROOTFLAGS) $(ROOTLIBS) -mmacosx-version-min=12.6

ThroughPlaneAnalysis_Metaslope: ThroughPlaneAnalysis_Metaslope.cc
	g++ -O2 ThroughPlaneAnalysis_Metaslope.cc -o ThroughPlaneAnalysis_Metaslope $(ROOTFLAGS) $(ROOTLIBS) -mmacosx-version-min=12.6
