mf: mfMain.o mf.o subjectMerge.o generateAlign.o parse.o \
	convert.o classify.o globalSeg.o lenStatistics.o accuracyStatistics.o \
	util.o reads.o queryIndex.o map.o misass.o indel.o sv.o misinfo.o correct.o \
	circosData.o misReg.o
	gcc -g -o $@ $^ -lm -lpthread # $@: targets; $^: dependences
	rm *.o
.c.o:  # .o relys on the corresponding .c.
	gcc -g -c $<

clean:
	rm *.o
	rm *~
	rm *.bak
