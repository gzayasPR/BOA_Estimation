#CFLAGS = -O3
CFLAGS = -g -O3


CC = g++ -g -O3
#CC = gcc

all: triolanc unolanc haplanc gettriophase unolanc2way random.o unolanc5way trainhmm

random.o: random.cpp
	$(CC) $(CFLAGS) -c random.cpp 
unolanc: unolanc.cpp random.o hmm_phase.h
	$(CC) $(CFLAGS) random.o unolanc.cpp -lm -o unolanc
trainhmm: trainhmm.cpp random.o hmm_phase.h
	$(CC) $(CFLAGS) random.o trainhmm.cpp -lm -o trainhmm
unolanc2way: unolanc2way.cpp random.o hmm_phase.h
	$(CC) $(CFLAGS) random.o unolanc2way.cpp -lm -o unolanc2way
unolanc5way: unolanc5way.cpp random.o hmm_phase.h
	$(CC) $(CFLAGS) random.o unolanc5way.cpp -lm -o unolanc5way
haplanc: haplanc.cpp random.o hmm_phase.h
	$(CC) $(CFLAGS) random.o haplanc.cpp -lm -o haplanc
triolanc: triolanc.cpp random.o hmm_phase.h
	$(CC) $(CFLAGS) triolanc.cpp random.o -lm -o triolanc
gettriophase: gettriophase.cpp random.o hmm_phase.h
	$(CC) $(CFLAGS) random.o gettriophase.cpp -lm -o gettriophase
install:
	mkdir bin
	mv triolanc bin/ ; mv unolanc bin/ ; mv haplanc bin/ ; mv gettriophase bin/ ; mv unolanc2way bin/; mv unolanc5way bin/ ; mv trainhmm bin/
clean:	
	rm -f -v *.o
	rm -f haplanc unolanc triolanc gettriophase random.o unolanc2way unolanc5way
	rm -r bin
