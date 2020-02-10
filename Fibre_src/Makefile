CC = g++
CFLAG = -std=c++11 -fopenmp 
all : main analyze
.PHONY : all
main : main.cpp *.h
	$(CC) $(CFLAG) -o main main.cpp
analyze : analyze.cpp *.h
	$(CC) $(CFLAG) -o analyze analyze.cpp
clean :
	rm main analyze
install : main analyze
	mv main analyze ./tmp