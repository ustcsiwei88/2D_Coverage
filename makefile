## Needs to differentiate nmake in win and make

SYS=LINUX
GUROBI_INCLUDE=/opt/gurobi910/linux64/include/
GUROBI_LIB=/opt/gurobi910/linux64/lib/
ifeq ($(SYS), WIN)
	CXX=cl
	CXXFLAGS=/permissive- /GS /GL /W3 /Gy /Zc:wchar_t /I"C:\gurobi900\win64\include" \
		/Zi /Gm- /O2 /Zc:inline /fp:precise /D "NDEBUG" \
		/D "_CONSOLE" /D "_MBCS" /errorReport:prompt /WX- /Zc:forScope /Gd /Oi /MT \
		/FC /EHsc /nologo  \
		/diagnostics:classic
	LINKFLAGS= /link "C:\gurobi900\win64\lib\gurobi_c++mt2017.lib" "C:\gurobi900\win64\lib\gurobi90.lib" \
                /LIBPATH:"C:\gurobi900\win64\lib"
else
	CXX=g++
	CXXFLAGS= -I$(GUROBI_INCLUDE) -Iinclude/ -O2 -m64 
	LINKFLAGS= -L$(GUROBI_LIB) -Wl,-rpath=$(GUROBI_LIB) -lgurobi_c++ -lgurobi91 -lm 
endif


polygen: polygen.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LINKFLAGS) 

org2d: ORG2D.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LINKFLAGS) 

opg2d: OPG2D.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LINKFLAGS) 

all: opg2d org2d polygen


#main: main.cpp
#	cl /std:c++17 /O2 /EHsc /openmp main.cpp


#ilp: ILP.cpp
#	cl /permissive- /GS /GL /W3 /Gy /Zc:wchar_t /I"C:\gurobi900\win64\include" \
#		/Zi /Gm- /O2 /Zc:inline /fp:precise /D "NDEBUG" \
#		/D "_CONSOLE" /D "_MBCS" /errorReport:prompt /WX- /Zc:forScope /Gd /Oi /MT \
#		/FC /EHsc /nologo  \
#		/diagnostics:classic \
#		ILP.cpp \
#		/link "C:\gurobi900\win64\lib\gurobi_c++mt2017.lib" "C:\gurobi900\win64\lib\gurobi90.lib" \
#		/LIBPATH:"C:\gurobi900\win64\lib"
#gen_polygen: polygen.cpp
#	cl /permissive- /GS /GL /W3 /Gy /Zc:wchar_t /I"C:\gurobi900\win64\include" \
#		/Zi /Gm- /O2 /Zc:inline /fp:precise /D "NDEBUG" \
#		/D "_CONSOLE" /D "_MBCS" /errorReport:prompt /WX- /Zc:forScope /Gd /Oi /MT \
#		/FC /EHsc /nologo  \
#		/diagnostics:classic \
#		polygen.cpp \
#		/link "C:\gurobi900\win64\lib\gurobi_c++mt2017.lib" "C:\gurobi900\win64\lib\gurobi90.lib" \
#		/LIBPATH:"C:\gurobi900\win64\lib"



# ilplink:
# 	link ILP.obj /OUT:"C:\Users\siwei\2D-geometry\IPsolution\x64\Release\ConsoleApplication1.exe" /MANIFEST /LTCG:incremental /NXCOMPAT /PDB:"C:\Users\siwei\2D-geometry\IPsolution\x64\Release\ConsoleApplication1.pdb" /DYNAMICBASE "gurobi_c++mt2017.lib" "gurobi90.lib" "kernel32.lib" "user32.lib" "gdi32.lib" "winspool.lib" "comdlg32.lib" "advapi32.lib" "shell32.lib" "ole32.lib" "oleaut32.lib" "uuid.lib" "odbc32.lib" "odbccp32.lib" /DEBUG /MACHINE:X64 /OPT:REF /PGD:"C:\Users\siwei\2D-geometry\IPsolution\x64\Release\ConsoleApplication1.pgd" /SUBSYSTEM:CONSOLE /MANIFESTUAC:"level='asInvoker' uiAccess='false'" /ManifestFile:"x64\Release\ConsoleApplication1.exe.intermediate.manifest" /OPT:ICF /ERRORREPORT:PROMPT /NOLOGO /LIBPATH:"C:\gurobi900\win64\lib" /TLBID:1 
# 	link /MANIFEST /LTCG:incremental /NXCOMPAT /DYNAMICBASE "gurobi_c++mt2017.lib" \
# 		"gurobi90.lib" "kernel32.lib" "user32.lib" "gdi32.lib" "winspool.lib" "comdlg32.lib"\
# 		"advapi32.lib" "shell32.lib" "ole32.lib" "oleaut32.lib" "uuid.lib" "odbc32.lib"\
# 		"odbccp32.lib" /DEBUG /MACHINE:X64 /OPT:REF /SUBSYSTEM:CONSOLE\
# 		/MANIFESTUAC:"level='asInvoker' uiAccess='false'" \
# 		/OPT:ICF /ERRORREPORT:PROMPT /NOLOGO /LIBPATH:"C:\gurobi900\win64\lib" /TLBID:1 
#maingcc:
#	g++ -O2 -std=c++11 main.cpp -fopenmp
#minenc:
#	cl /O2 /EHsc minenclosingdisc.cpp
