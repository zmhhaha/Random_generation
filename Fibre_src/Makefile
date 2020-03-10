CC = g++
CXXFLAG = -std=c++11 -fopenmp
all : main analyze
.PHONY : all
main : src/main.cpp src/*.h
	$(CC) $(CXXFLAG) -g -o main src/main.cpp 
analyze : src/analyze.cpp src/*.h
	$(CC) $(CXXFLAG) -g -o analyze src/analyze.cpp
clean :
	rm main analyze
install : main analyze
	cp main analyze /mnt/sdb1
	mv main analyze ./bin