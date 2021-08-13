All:
	#rootcint -f Dict_analysis.cpp -c `root-config --cflags` -p include/DRSWAVE.h include/mcp.h include/energy.h include/timing.h include/testbeam.h LinkDef.h 
	rootcint -f Dict_analysis.cpp include/DRSWAVE.h include/mcp.h include/energy.h include/timing.h include/testbeam.h LinkDef.h 
	g++ analysis.cc src/testbeam.cxx src/mcp.cxx src/energy.cxx src/timing.cxx src/DRSWAVE.cxx Dict_analysis.cpp -Iinclude `root-config --cflags --glibs`
