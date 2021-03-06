#############################################################################
# Makefile for building: libsrc.a
# Generated by qmake (2.01a) (Qt 4.4.1) on: Thu Aug 27 21:26:13 2009
# Project:  src.pro
# Template: lib
# Command: /usr/bin/qmake -spec /usr/local/Qt4.4/mkspecs/darwin-g++ -macx -o Makefile src.pro
#############################################################################

# simplified by rlf 

####### Compiler, tools and options

CC            = cc
CXX           = g++
DEFINES       = -D__USE_WS_X11__
CFLAGS        = -pipe -O2 -fPIC $(DEFINES)
CXXFLAGS      = -pipe -O2 -fPIC $(DEFINES) -arch x86_64
INCPATH       = -I/usr/local/Qt4.4/mkspecs/darwin-g++ -I. -I.
AR            = ar cq
RANLIB        = ranlib -s
QMAKE         = /usr/bin/qmake
TAR           = tar -cf
COMPRESS      = gzip -9f
COPY          = cp -f
SED           = sed
COPY_FILE     = cp -f
COPY_DIR      = cp -f -r
INSTALL_FILE  = $(COPY_FILE)
INSTALL_DIR   = $(COPY_DIR)
INSTALL_PROGRAM = $(COPY_FILE)
DEL_FILE      = rm -f
SYMLINK       = ln -sf
DEL_DIR       = rmdir
MOVE          = mv -f
CHK_DIR_EXISTS= test -d
MKDIR         = mkdir -p

####### Output directory

OBJECTS_DIR   = obj/

####### Files

SOURCES       = alleletrim.cpp \
		bg.cpp \
		bgmatrix.cpp \
		BNodeList.cpp \
		chromosome.cpp \
		data.cpp \
		datanode.cpp \
		dblock.cpp \
		doublematrix.cpp \
		field.cpp \
		fieldstruct.cpp \
		fpmm.cpp \
		geneticdist.cpp \
		genome.cpp \
		ginverse1.cpp \
		glmm1.cpp \
		glmm2.cpp \
		gnodederived.cpp \
		gnodestuff.cpp \
		hashtable.cpp \
		histogram.cpp \
		IBDMatrix.cpp \
		individual.cpp \
		maternalPaternal.cpp \
		mim.cpp \
		min_nr_powell.cpp \
		min_praxis.cpp \
		model.cpp \
		model_blup.cpp \
		model_gibbs.cpp \
		model_hypo.cpp \
		model_peeling.cpp \
		model_vce.cpp \
		mqtl.cpp \
		nonsymm_eigen.cpp \
		nufamily.cpp \
		parmMap.cpp \
		pedigree.cpp \
		plotter.cpp \
		pop_bayes.cpp \
		pop_gibbs.cpp \
		pop_graph.cpp \
		pop_mblup.cpp \
		pop_peeling.cpp \
		popMIM.cpp \
		population.cpp \
		population1.cpp \
		rpedigree.cpp \
		rsamplerparms.cpp \
		rutil.cpp \
		session.cpp \
		sparsebgmatrix.cpp \
		sparsematrix.cpp \
		stat1.cpp \
		stat2.cpp \
		statdist.cpp \
		symm_eigen.cpp \
		termlist.cpp \
		util.cpp 
OBJECTS       = obj/alleletrim.o \
		obj/bg.o \
		obj/bgmatrix.o \
		obj/BNodeList.o \
		obj/chromosome.o \
		obj/data.o \
		obj/datanode.o \
		obj/dblock.o \
		obj/doublematrix.o \
		obj/field.o \
		obj/fieldstruct.o \
		obj/fpmm.o \
		obj/geneticdist.o \
		obj/genome.o \
		obj/ginverse1.o \
		obj/glmm1.o \
		obj/glmm2.o \
		obj/gnodederived.o \
		obj/gnodestuff.o \
		obj/hashtable.o \
		obj/histogram.o \
		obj/IBDMatrix.o \
		obj/individual.o \
		obj/maternalPaternal.o \
		obj/mim.o \
		obj/min_nr_powell.o \
		obj/min_praxis.o \
		obj/model.o \
		obj/model_blup.o \
		obj/model_gibbs.o \
		obj/model_hypo.o \
		obj/model_peeling.o \
		obj/model_vce.o \
		obj/mqtl.o \
		obj/nonsymm_eigen.o \
		obj/nufamily.o \
		obj/parmMap.o \
		obj/pedigree.o \
		obj/plotter.o \
		obj/pop_bayes.o \
		obj/pop_gibbs.o \
		obj/pop_graph.o \
		obj/pop_mblup.o \
		obj/pop_peeling.o \
		obj/popMIM.o \
		obj/population.o \
		obj/population1.o \
		obj/rpedigree.o \
		obj/rsamplerparms.o \
		obj/rutil.o \
		obj/session.o \
		obj/sparsebgmatrix.o \
		obj/sparsematrix.o \
		obj/stat1.o \
		obj/stat2.o \
		obj/statdist.o \
		obj/symm_eigen.o \
		obj/termlist.o \
		obj/util.o

TARGET        = libmatvec.a

####### Build rules


$(TARGET):  $(OBJECTS) $(OBJCOMP) 
	-$(DEL_FILE) $(TARGET)
	$(AR) $(TARGET) $(OBJECTS)
	$(RANLIB) $(TARGET)

clean: 
	-$(DEL_FILE) $(OBJECTS) $(TARGET)
	-$(DEL_FILE) *~ core *.core

####### Compile

obj/alleletrim.o: alleletrim.cpp safe_vectors.h \
		session.h \
		MersenneTwister.h \
		vector.h \
		exception.h \
		util.h \
		stat.h \
		gnodederived.h \
		gnodestuff.h \
		geneticdist.h \
		doublematrix.h \
		matrix.h \
		gnodesetderived.h \
		individual.h \
		datanode.h \
		genome.h \
		chromosome.h \
		dblock.h \
		bg.h \
		population.h \
		hashtable.h \
		pedigree.h \
		rpedigree.h \
		rutil.h \
		data.h \
		field.h \
		fieldstruct.h \
		mim.h \
		parmMap.h \
		mqtl.h \
		IBDMatrix.h \
		BNodeList.h \
		nufamily.h \
		sparsematrix.h \
		fpmm.h \
		model.h \
		termlist.h \
		minimizer.h \
		rsamplerparms.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/alleletrim.o alleletrim.cpp

obj/bg.o: bg.cpp bg.h \
		doublematrix.h \
		session.h \
		MersenneTwister.h \
		matrix.h \
		vector.h \
		exception.h \
		util.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/bg.o bg.cpp

obj/bgmatrix.o: bgmatrix.cpp session.h \
		MersenneTwister.h \
		util.h \
		bgmatrix.h \
		matrix.h \
		vector.h \
		exception.h \
		bg.h \
		doublematrix.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/bgmatrix.o bgmatrix.cpp

obj/BNodeList.o: BNodeList.cpp stat.h \
		gnodederived.h \
		exception.h \
		gnodestuff.h \
		geneticdist.h \
		util.h \
		doublematrix.h \
		session.h \
		MersenneTwister.h \
		matrix.h \
		vector.h \
		safe_vectors.h \
		gnodesetderived.h \
		individual.h \
		datanode.h \
		genome.h \
		chromosome.h \
		dblock.h \
		bg.h \
		population.h \
		hashtable.h \
		pedigree.h \
		rpedigree.h \
		rutil.h \
		data.h \
		field.h \
		fieldstruct.h \
		mim.h \
		parmMap.h \
		mqtl.h \
		IBDMatrix.h \
		BNodeList.h \
		nufamily.h \
		sparsematrix.h \
		fpmm.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/BNodeList.o BNodeList.cpp

obj/chromosome.o: chromosome.cpp util.h \
		chromosome.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/chromosome.o chromosome.cpp

obj/data.o: data.cpp session.h \
		MersenneTwister.h \
		util.h \
		doublematrix.h \
		matrix.h \
		vector.h \
		exception.h \
		data.h \
		datanode.h \
		field.h \
		fieldstruct.h \
		hashtable.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/data.o data.cpp

obj/datanode.o: datanode.cpp datanode.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/datanode.o datanode.cpp

obj/dblock.o: dblock.cpp dblock.h \
		doublematrix.h \
		session.h \
		MersenneTwister.h \
		matrix.h \
		vector.h \
		exception.h \
		util.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/dblock.o dblock.cpp

obj/doublematrix.o: doublematrix.cpp session.h \
		MersenneTwister.h \
		util.h \
		doublematrix.h \
		matrix.h \
		vector.h \
		exception.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/doublematrix.o doublematrix.cpp

obj/field.o: field.cpp field.h \
		util.h \
		vector.h \
		session.h \
		MersenneTwister.h \
		exception.h \
		datanode.h \
		fieldstruct.h \
		hashtable.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/field.o field.cpp

obj/fieldstruct.o: fieldstruct.cpp fieldstruct.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/fieldstruct.o fieldstruct.cpp

obj/fpmm.o: fpmm.cpp fpmm.h \
		dblock.h \
		doublematrix.h \
		session.h \
		MersenneTwister.h \
		matrix.h \
		vector.h \
		exception.h \
		util.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/fpmm.o fpmm.cpp

obj/geneticdist.o: geneticdist.cpp geneticdist.h \
		util.h \
		doublematrix.h \
		session.h \
		MersenneTwister.h \
		matrix.h \
		vector.h \
		exception.h \
		safe_vectors.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/geneticdist.o geneticdist.cpp

obj/genome.o: genome.cpp genome.h \
		geneticdist.h \
		util.h \
		doublematrix.h \
		session.h \
		MersenneTwister.h \
		matrix.h \
		vector.h \
		exception.h \
		safe_vectors.h \
		chromosome.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/genome.o genome.cpp

obj/ginverse1.o: ginverse1.cpp exception.h \
		bg.h \
		doublematrix.h \
		session.h \
		MersenneTwister.h \
		matrix.h \
		vector.h \
		util.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/ginverse1.o ginverse1.cpp

obj/glmm1.o: glmm1.cpp session.h \
		MersenneTwister.h \
		util.h \
		doublematrix.h \
		matrix.h \
		vector.h \
		exception.h \
		model.h \
		data.h \
		datanode.h \
		field.h \
		fieldstruct.h \
		hashtable.h \
		geneticdist.h \
		safe_vectors.h \
		termlist.h \
		minimizer.h \
		fpmm.h \
		dblock.h \
		sparsematrix.h \
		rsamplerparms.h \
		population.h \
		pedigree.h \
		rpedigree.h \
		rutil.h \
		individual.h \
		genome.h \
		chromosome.h \
		gnodederived.h \
		gnodestuff.h \
		bg.h \
		mim.h \
		parmMap.h \
		mqtl.h \
		IBDMatrix.h \
		BNodeList.h \
		nufamily.h \
		gnodesetderived.h \
		statdist.h \
		stat.h \
		statdistbase.h \
		glmm.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/glmm1.o glmm1.cpp

obj/glmm2.o: glmm2.cpp session.h \
		MersenneTwister.h \
		util.h \
		doublematrix.h \
		matrix.h \
		vector.h \
		exception.h \
		model.h \
		data.h \
		datanode.h \
		field.h \
		fieldstruct.h \
		hashtable.h \
		geneticdist.h \
		safe_vectors.h \
		termlist.h \
		minimizer.h \
		fpmm.h \
		dblock.h \
		sparsematrix.h \
		rsamplerparms.h \
		population.h \
		pedigree.h \
		rpedigree.h \
		rutil.h \
		individual.h \
		genome.h \
		chromosome.h \
		gnodederived.h \
		gnodestuff.h \
		bg.h \
		mim.h \
		parmMap.h \
		mqtl.h \
		IBDMatrix.h \
		BNodeList.h \
		nufamily.h \
		gnodesetderived.h \
		statdist.h \
		stat.h \
		statdistbase.h \
		glmm.h \
		pt.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/glmm2.o glmm2.cpp

obj/gnodederived.o: gnodederived.cpp safe_vectors.h \
		session.h \
		MersenneTwister.h \
		vector.h \
		exception.h \
		util.h \
		gnodederived.h \
		gnodestuff.h \
		geneticdist.h \
		doublematrix.h \
		matrix.h \
		gnodesetderived.h \
		individual.h \
		datanode.h \
		genome.h \
		chromosome.h \
		dblock.h \
		bg.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/gnodederived.o gnodederived.cpp

obj/gnodestuff.o: gnodestuff.cpp safe_vectors.h \
		session.h \
		MersenneTwister.h \
		vector.h \
		exception.h \
		util.h \
		stat.h \
		population.h \
		hashtable.h \
		pedigree.h \
		doublematrix.h \
		matrix.h \
		rpedigree.h \
		rutil.h \
		data.h \
		datanode.h \
		field.h \
		fieldstruct.h \
		individual.h \
		genome.h \
		geneticdist.h \
		chromosome.h \
		gnodederived.h \
		gnodestuff.h \
		dblock.h \
		bg.h \
		mim.h \
		parmMap.h \
		mqtl.h \
		IBDMatrix.h \
		BNodeList.h \
		nufamily.h \
		gnodesetderived.h \
		sparsematrix.h \
		fpmm.h \
		model.h \
		termlist.h \
		minimizer.h \
		rsamplerparms.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/gnodestuff.o gnodestuff.cpp

obj/hashtable.o: hashtable.cpp exception.h \
		util.h \
		hashtable.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/hashtable.o hashtable.cpp

obj/histogram.o: histogram.cpp histogram.h \
		safe_vectors.h \
		session.h \
		MersenneTwister.h \
		vector.h \
		exception.h \
		util.h \
		plotter.h \
		doublematrix.h \
		matrix.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/histogram.o histogram.cpp

obj/IBDMatrix.o: IBDMatrix.cpp IBDMatrix.h \
		safe_vectors.h \
		session.h \
		MersenneTwister.h \
		vector.h \
		exception.h \
		util.h \
		matrix.h \
		population.h \
		hashtable.h \
		pedigree.h \
		doublematrix.h \
		rpedigree.h \
		rutil.h \
		data.h \
		datanode.h \
		field.h \
		fieldstruct.h \
		individual.h \
		genome.h \
		geneticdist.h \
		chromosome.h \
		gnodederived.h \
		gnodestuff.h \
		dblock.h \
		bg.h \
		mim.h \
		parmMap.h \
		mqtl.h \
		BNodeList.h \
		nufamily.h \
		gnodesetderived.h \
		sparsematrix.h \
		fpmm.h \
		stat.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/IBDMatrix.o IBDMatrix.cpp

obj/individual.o: individual.cpp population.h \
		hashtable.h \
		pedigree.h \
		doublematrix.h \
		session.h \
		MersenneTwister.h \
		matrix.h \
		vector.h \
		exception.h \
		util.h \
		rpedigree.h \
		rutil.h \
		safe_vectors.h \
		data.h \
		datanode.h \
		field.h \
		fieldstruct.h \
		individual.h \
		genome.h \
		geneticdist.h \
		chromosome.h \
		gnodederived.h \
		gnodestuff.h \
		dblock.h \
		bg.h \
		mim.h \
		parmMap.h \
		mqtl.h \
		IBDMatrix.h \
		BNodeList.h \
		nufamily.h \
		gnodesetderived.h \
		sparsematrix.h \
		fpmm.h \
		model.h \
		termlist.h \
		minimizer.h \
		rsamplerparms.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/individual.o individual.cpp

obj/maternalPaternal.o: maternalPaternal.cpp gnodederived.h \
		exception.h \
		gnodestuff.h \
		geneticdist.h \
		util.h \
		doublematrix.h \
		session.h \
		MersenneTwister.h \
		matrix.h \
		vector.h \
		safe_vectors.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/maternalPaternal.o maternalPaternal.cpp

obj/mim.o: mim.cpp mim.h \
		doublematrix.h \
		session.h \
		MersenneTwister.h \
		matrix.h \
		vector.h \
		exception.h \
		util.h \
		geneticdist.h \
		safe_vectors.h \
		parmMap.h \
		rutil.h \
		rpedigree.h \
		mqtl.h \
		IBDMatrix.h \
		individual.h \
		datanode.h \
		genome.h \
		chromosome.h \
		gnodederived.h \
		gnodestuff.h \
		dblock.h \
		bg.h \
		BNodeList.h \
		rsamplerparms.h \
		population.h \
		hashtable.h \
		pedigree.h \
		data.h \
		field.h \
		fieldstruct.h \
		nufamily.h \
		gnodesetderived.h \
		sparsematrix.h \
		fpmm.h \
		stat.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/mim.o mim.cpp

obj/min_nr_powell.o: min_nr_powell.cpp model.h \
		session.h \
		MersenneTwister.h \
		data.h \
		datanode.h \
		field.h \
		util.h \
		vector.h \
		exception.h \
		fieldstruct.h \
		hashtable.h \
		geneticdist.h \
		doublematrix.h \
		matrix.h \
		safe_vectors.h \
		termlist.h \
		minimizer.h \
		fpmm.h \
		dblock.h \
		sparsematrix.h \
		rsamplerparms.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/min_nr_powell.o min_nr_powell.cpp

obj/min_praxis.o: min_praxis.cpp minimizer.h \
		matrix.h \
		session.h \
		MersenneTwister.h \
		vector.h \
		exception.h \
		util.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/min_praxis.o min_praxis.cpp

obj/model.o: model.cpp histogram.h \
		safe_vectors.h \
		session.h \
		MersenneTwister.h \
		vector.h \
		exception.h \
		util.h \
		plotter.h \
		doublematrix.h \
		matrix.h \
		stat.h \
		model.h \
		data.h \
		datanode.h \
		field.h \
		fieldstruct.h \
		hashtable.h \
		geneticdist.h \
		termlist.h \
		minimizer.h \
		fpmm.h \
		dblock.h \
		sparsematrix.h \
		rsamplerparms.h \
		population.h \
		pedigree.h \
		rpedigree.h \
		rutil.h \
		individual.h \
		genome.h \
		chromosome.h \
		gnodederived.h \
		gnodestuff.h \
		bg.h \
		mim.h \
		parmMap.h \
		mqtl.h \
		IBDMatrix.h \
		BNodeList.h \
		nufamily.h \
		gnodesetderived.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/model.o model.cpp

obj/model_blup.o: model_blup.cpp session.h \
		MersenneTwister.h \
		model.h \
		data.h \
		datanode.h \
		field.h \
		util.h \
		vector.h \
		exception.h \
		fieldstruct.h \
		hashtable.h \
		geneticdist.h \
		doublematrix.h \
		matrix.h \
		safe_vectors.h \
		termlist.h \
		minimizer.h \
		fpmm.h \
		dblock.h \
		sparsematrix.h \
		rsamplerparms.h \
		population.h \
		pedigree.h \
		rpedigree.h \
		rutil.h \
		individual.h \
		genome.h \
		chromosome.h \
		gnodederived.h \
		gnodestuff.h \
		bg.h \
		mim.h \
		parmMap.h \
		mqtl.h \
		IBDMatrix.h \
		BNodeList.h \
		nufamily.h \
		gnodesetderived.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/model_blup.o model_blup.cpp

obj/model_gibbs.o: model_gibbs.cpp model.h \
		session.h \
		MersenneTwister.h \
		data.h \
		datanode.h \
		field.h \
		util.h \
		vector.h \
		exception.h \
		fieldstruct.h \
		hashtable.h \
		geneticdist.h \
		doublematrix.h \
		matrix.h \
		safe_vectors.h \
		termlist.h \
		minimizer.h \
		fpmm.h \
		dblock.h \
		sparsematrix.h \
		rsamplerparms.h \
		stat.h \
		population.h \
		pedigree.h \
		rpedigree.h \
		rutil.h \
		individual.h \
		genome.h \
		chromosome.h \
		gnodederived.h \
		gnodestuff.h \
		bg.h \
		mim.h \
		parmMap.h \
		mqtl.h \
		IBDMatrix.h \
		BNodeList.h \
		nufamily.h \
		gnodesetderived.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/model_gibbs.o model_gibbs.cpp

obj/model_hypo.o: model_hypo.cpp session.h \
		MersenneTwister.h \
		util.h \
		model.h \
		data.h \
		datanode.h \
		field.h \
		vector.h \
		exception.h \
		fieldstruct.h \
		hashtable.h \
		geneticdist.h \
		doublematrix.h \
		matrix.h \
		safe_vectors.h \
		termlist.h \
		minimizer.h \
		fpmm.h \
		dblock.h \
		sparsematrix.h \
		rsamplerparms.h \
		stat.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/model_hypo.o model_hypo.cpp

obj/model_peeling.o: model_peeling.cpp model.h \
		session.h \
		MersenneTwister.h \
		data.h \
		datanode.h \
		field.h \
		util.h \
		vector.h \
		exception.h \
		fieldstruct.h \
		hashtable.h \
		geneticdist.h \
		doublematrix.h \
		matrix.h \
		safe_vectors.h \
		termlist.h \
		minimizer.h \
		fpmm.h \
		dblock.h \
		sparsematrix.h \
		rsamplerparms.h \
		population.h \
		pedigree.h \
		rpedigree.h \
		rutil.h \
		individual.h \
		genome.h \
		chromosome.h \
		gnodederived.h \
		gnodestuff.h \
		bg.h \
		mim.h \
		parmMap.h \
		mqtl.h \
		IBDMatrix.h \
		BNodeList.h \
		nufamily.h \
		gnodesetderived.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/model_peeling.o model_peeling.cpp

obj/model_vce.o: model_vce.cpp session.h \
		MersenneTwister.h \
		model.h \
		data.h \
		datanode.h \
		field.h \
		util.h \
		vector.h \
		exception.h \
		fieldstruct.h \
		hashtable.h \
		geneticdist.h \
		doublematrix.h \
		matrix.h \
		safe_vectors.h \
		termlist.h \
		minimizer.h \
		fpmm.h \
		dblock.h \
		sparsematrix.h \
		rsamplerparms.h \
		stat.h \
		population.h \
		pedigree.h \
		rpedigree.h \
		rutil.h \
		individual.h \
		genome.h \
		chromosome.h \
		gnodederived.h \
		gnodestuff.h \
		bg.h \
		mim.h \
		parmMap.h \
		mqtl.h \
		IBDMatrix.h \
		BNodeList.h \
		nufamily.h \
		gnodesetderived.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/model_vce.o model_vce.cpp

obj/mqtl.o: mqtl.cpp mqtl.h \
		rpedigree.h \
		rutil.h \
		safe_vectors.h \
		session.h \
		MersenneTwister.h \
		vector.h \
		exception.h \
		util.h \
		doublematrix.h \
		matrix.h \
		mim.h \
		geneticdist.h \
		parmMap.h \
		IBDMatrix.h \
		individual.h \
		datanode.h \
		genome.h \
		chromosome.h \
		gnodederived.h \
		gnodestuff.h \
		dblock.h \
		bg.h \
		BNodeList.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/mqtl.o mqtl.cpp

obj/nonsymm_eigen.o: nonsymm_eigen.cpp util.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/nonsymm_eigen.o nonsymm_eigen.cpp

obj/nufamily.o: nufamily.cpp nufamily.h \
		individual.h \
		doublematrix.h \
		session.h \
		MersenneTwister.h \
		matrix.h \
		vector.h \
		exception.h \
		util.h \
		datanode.h \
		genome.h \
		geneticdist.h \
		safe_vectors.h \
		chromosome.h \
		gnodederived.h \
		gnodestuff.h \
		dblock.h \
		bg.h \
		population.h \
		hashtable.h \
		pedigree.h \
		rpedigree.h \
		rutil.h \
		data.h \
		field.h \
		fieldstruct.h \
		mim.h \
		parmMap.h \
		mqtl.h \
		IBDMatrix.h \
		BNodeList.h \
		gnodesetderived.h \
		sparsematrix.h \
		fpmm.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/nufamily.o nufamily.cpp

obj/parmMap.o: parmMap.cpp parmMap.h \
		exception.h \
		rutil.h \
		safe_vectors.h \
		session.h \
		MersenneTwister.h \
		vector.h \
		util.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/parmMap.o parmMap.cpp

obj/pedigree.o: pedigree.cpp session.h \
		MersenneTwister.h \
		util.h \
		stat.h \
		hashtable.h \
		pedigree.h \
		doublematrix.h \
		matrix.h \
		vector.h \
		exception.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/pedigree.o pedigree.cpp

obj/plotter.o: plotter.cpp session.h \
		MersenneTwister.h \
		util.h \
		plotter.h \
		doublematrix.h \
		matrix.h \
		vector.h \
		exception.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/plotter.o plotter.cpp

obj/pop_bayes.o: pop_bayes.cpp model.h \
		session.h \
		MersenneTwister.h \
		data.h \
		datanode.h \
		field.h \
		util.h \
		vector.h \
		exception.h \
		fieldstruct.h \
		hashtable.h \
		geneticdist.h \
		doublematrix.h \
		matrix.h \
		safe_vectors.h \
		termlist.h \
		minimizer.h \
		fpmm.h \
		dblock.h \
		sparsematrix.h \
		rsamplerparms.h \
		population.h \
		pedigree.h \
		rpedigree.h \
		rutil.h \
		individual.h \
		genome.h \
		chromosome.h \
		gnodederived.h \
		gnodestuff.h \
		bg.h \
		mim.h \
		parmMap.h \
		mqtl.h \
		IBDMatrix.h \
		BNodeList.h \
		nufamily.h \
		gnodesetderived.h \
		statdist.h \
		stat.h \
		statdistbase.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/pop_bayes.o pop_bayes.cpp

obj/pop_gibbs.o: pop_gibbs.cpp population.h \
		hashtable.h \
		pedigree.h \
		doublematrix.h \
		session.h \
		MersenneTwister.h \
		matrix.h \
		vector.h \
		exception.h \
		util.h \
		rpedigree.h \
		rutil.h \
		safe_vectors.h \
		data.h \
		datanode.h \
		field.h \
		fieldstruct.h \
		individual.h \
		genome.h \
		geneticdist.h \
		chromosome.h \
		gnodederived.h \
		gnodestuff.h \
		dblock.h \
		bg.h \
		mim.h \
		parmMap.h \
		mqtl.h \
		IBDMatrix.h \
		BNodeList.h \
		nufamily.h \
		gnodesetderived.h \
		sparsematrix.h \
		fpmm.h \
		stat.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/pop_gibbs.o pop_gibbs.cpp

obj/pop_graph.o: pop_graph.cpp session.h \
		MersenneTwister.h \
		population.h \
		hashtable.h \
		pedigree.h \
		doublematrix.h \
		matrix.h \
		vector.h \
		exception.h \
		util.h \
		rpedigree.h \
		rutil.h \
		safe_vectors.h \
		data.h \
		datanode.h \
		field.h \
		fieldstruct.h \
		individual.h \
		genome.h \
		geneticdist.h \
		chromosome.h \
		gnodederived.h \
		gnodestuff.h \
		dblock.h \
		bg.h \
		mim.h \
		parmMap.h \
		mqtl.h \
		IBDMatrix.h \
		BNodeList.h \
		nufamily.h \
		gnodesetderived.h \
		sparsematrix.h \
		fpmm.h \
		stat.h \
		model.h \
		termlist.h \
		minimizer.h \
		rsamplerparms.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/pop_graph.o pop_graph.cpp

obj/pop_mblup.o: pop_mblup.cpp session.h \
		MersenneTwister.h \
		util.h \
		population.h \
		hashtable.h \
		pedigree.h \
		doublematrix.h \
		matrix.h \
		vector.h \
		exception.h \
		rpedigree.h \
		rutil.h \
		safe_vectors.h \
		data.h \
		datanode.h \
		field.h \
		fieldstruct.h \
		individual.h \
		genome.h \
		geneticdist.h \
		chromosome.h \
		gnodederived.h \
		gnodestuff.h \
		dblock.h \
		bg.h \
		mim.h \
		parmMap.h \
		mqtl.h \
		IBDMatrix.h \
		BNodeList.h \
		nufamily.h \
		gnodesetderived.h \
		sparsematrix.h \
		fpmm.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/pop_mblup.o pop_mblup.cpp

obj/pop_peeling.o: pop_peeling.cpp session.h \
		MersenneTwister.h \
		population.h \
		hashtable.h \
		pedigree.h \
		doublematrix.h \
		matrix.h \
		vector.h \
		exception.h \
		util.h \
		rpedigree.h \
		rutil.h \
		safe_vectors.h \
		data.h \
		datanode.h \
		field.h \
		fieldstruct.h \
		individual.h \
		genome.h \
		geneticdist.h \
		chromosome.h \
		gnodederived.h \
		gnodestuff.h \
		dblock.h \
		bg.h \
		mim.h \
		parmMap.h \
		mqtl.h \
		IBDMatrix.h \
		BNodeList.h \
		nufamily.h \
		gnodesetderived.h \
		sparsematrix.h \
		fpmm.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/pop_peeling.o pop_peeling.cpp

obj/popMIM.o: popMIM.cpp model.h \
		session.h \
		MersenneTwister.h \
		data.h \
		datanode.h \
		field.h \
		util.h \
		vector.h \
		exception.h \
		fieldstruct.h \
		hashtable.h \
		geneticdist.h \
		doublematrix.h \
		matrix.h \
		safe_vectors.h \
		termlist.h \
		minimizer.h \
		fpmm.h \
		dblock.h \
		sparsematrix.h \
		rsamplerparms.h \
		population.h \
		pedigree.h \
		rpedigree.h \
		rutil.h \
		individual.h \
		genome.h \
		chromosome.h \
		gnodederived.h \
		gnodestuff.h \
		bg.h \
		mim.h \
		parmMap.h \
		mqtl.h \
		IBDMatrix.h \
		BNodeList.h \
		nufamily.h \
		gnodesetderived.h \
		stat.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/popMIM.o popMIM.cpp

obj/population.o: population.cpp model.h \
		session.h \
		MersenneTwister.h \
		data.h \
		datanode.h \
		field.h \
		util.h \
		vector.h \
		exception.h \
		fieldstruct.h \
		hashtable.h \
		geneticdist.h \
		doublematrix.h \
		matrix.h \
		safe_vectors.h \
		termlist.h \
		minimizer.h \
		fpmm.h \
		dblock.h \
		sparsematrix.h \
		rsamplerparms.h \
		population.h \
		pedigree.h \
		rpedigree.h \
		rutil.h \
		individual.h \
		genome.h \
		chromosome.h \
		gnodederived.h \
		gnodestuff.h \
		bg.h \
		mim.h \
		parmMap.h \
		mqtl.h \
		IBDMatrix.h \
		BNodeList.h \
		nufamily.h \
		gnodesetderived.h \
		stat.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/population.o population.cpp

obj/population1.o: population1.cpp population.h \
		hashtable.h \
		pedigree.h \
		doublematrix.h \
		session.h \
		MersenneTwister.h \
		matrix.h \
		vector.h \
		exception.h \
		util.h \
		rpedigree.h \
		rutil.h \
		safe_vectors.h \
		data.h \
		datanode.h \
		field.h \
		fieldstruct.h \
		individual.h \
		genome.h \
		geneticdist.h \
		chromosome.h \
		gnodederived.h \
		gnodestuff.h \
		dblock.h \
		bg.h \
		mim.h \
		parmMap.h \
		mqtl.h \
		IBDMatrix.h \
		BNodeList.h \
		nufamily.h \
		gnodesetderived.h \
		sparsematrix.h \
		fpmm.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/population1.o population1.cpp

obj/rpedigree.o: rpedigree.cpp rpedigree.h \
		rutil.h \
		safe_vectors.h \
		session.h \
		MersenneTwister.h \
		vector.h \
		exception.h \
		util.h \
		doublematrix.h \
		matrix.h \
		mim.h \
		geneticdist.h \
		parmMap.h \
		mqtl.h \
		IBDMatrix.h \
		individual.h \
		datanode.h \
		genome.h \
		chromosome.h \
		gnodederived.h \
		gnodestuff.h \
		dblock.h \
		bg.h \
		BNodeList.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/rpedigree.o rpedigree.cpp

obj/rsamplerparms.o: rsamplerparms.cpp rsamplerparms.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/rsamplerparms.o rsamplerparms.cpp

obj/rutil.o: rutil.cpp rutil.h \
		safe_vectors.h \
		session.h \
		MersenneTwister.h \
		vector.h \
		exception.h \
		util.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/rutil.o rutil.cpp

obj/session.o: session.cpp exception.h \
		session.h \
		MersenneTwister.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/session.o session.cpp

obj/sparsebgmatrix.o: sparsebgmatrix.cpp session.h \
		MersenneTwister.h \
		util.h \
		sparsebgmatrix.h \
		doublematrix.h \
		matrix.h \
		vector.h \
		exception.h \
		dblock.h \
		bgmatrix.h \
		bg.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/sparsebgmatrix.o sparsebgmatrix.cpp

obj/sparsematrix.o: sparsematrix.cpp session.h \
		MersenneTwister.h \
		util.h \
		sparsematrix.h \
		doublematrix.h \
		matrix.h \
		vector.h \
		exception.h \
		dblock.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/sparsematrix.o sparsematrix.cpp

obj/stat1.o: stat1.cpp exception.h \
		session.h \
		MersenneTwister.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/stat1.o stat1.cpp

obj/stat2.o: stat2.cpp exception.h \
		session.h \
		MersenneTwister.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/stat2.o stat2.cpp

obj/statdist.o: statdist.cpp statdist.h \
		doublematrix.h \
		session.h \
		MersenneTwister.h \
		matrix.h \
		vector.h \
		exception.h \
		util.h \
		stat.h \
		statdistbase.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/statdist.o statdist.cpp

obj/symm_eigen.o: symm_eigen.cpp vector.h \
		session.h \
		MersenneTwister.h \
		exception.h \
		util.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/symm_eigen.o symm_eigen.cpp

obj/termlist.o: termlist.cpp util.h \
		termlist.h \
		geneticdist.h \
		doublematrix.h \
		session.h \
		MersenneTwister.h \
		matrix.h \
		vector.h \
		exception.h \
		safe_vectors.h \
		fieldstruct.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/termlist.o termlist.cpp

obj/util.o: util.cpp util.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o obj/util.o util.cpp

