CXX=g++
SAMTOOLS_ROOT=vendor/samtools-0.1.19
XGBOOST_ROOT=vendor/xgboost

STATIC_LIST=vendor/xgboost/lib/libxgboost.a vendor/xgboost/dmlc-core/libdmlc.a vendor/xgboost/rabit/lib/librabit.a
CXXFLAGS+=    -O2 -fopenmp
LDFLAGS+=    -L${SAMTOOLS_ROOT}
LIBS+=    -lm -lbam -lz -lpthread
INCLUDES+=    -I${SAMTOOLS_ROOT} -I${XGBOOST_ROOT}/include -I${XGBOOST_ROOT}/dmlc-core/include -I${XGBOOST_ROOT}/rabit/include
SOURCE = cmds scan distribution refseq polyscan param utilities homo window bamreader sample chi somatic md5 
OBJS= $(patsubst %,%.o,$(SOURCE))

%.o:%.cpp
	        $(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

all: msisensor-ct

msisensor-ct: 
		$(MAKE) -C ${SAMTOOLS_ROOT}
		cmake vendor/xgboost/
		$(MAKE) -C ${XGBOOST_ROOT}
		$(MAKE) $(OBJS)
		$(CXX) $(OBJS) $(STATIC_LIST) $(CXXFLAGS) $(LDFLAGS) $(LIBS) $(INCLUDES) -o $@

clean:
	        rm -f *.o msisensor-ct
		rm CMakeCache.txt
		rm -r CMakeFiles
		$(MAKE) -C ${SAMTOOLS_ROOT} clean
		$(MAKE) -C ${XGBOOST_ROOT} clean_all

