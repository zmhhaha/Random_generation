CC = g++
CFLAG = -std=c++11 -fopenmp 
all : main analyze
.PHONY : all
main : src/main.cpp src/*.h
	$(CC) $(CFLAG) -o main src/main.cpp 
analyze : src/analyze.cpp src/*.h
	$(CC) $(CFLAG) -o analyze src/analyze.cpp
clean :
	rm main analyze
install : main analyze
	mv main analyze ./bin