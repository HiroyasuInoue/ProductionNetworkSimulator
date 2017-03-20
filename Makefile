CC = x86_64-w64-mingw32-g++
TGT = $(SRC:%.cpp=%)
CXX_DEBUG_FLAGS = -g
CXX_RELEASE_FLAGS = -O3

all: clean release

release: proNet.exe

proNet.exe: proNet.cpp
	$(CC) -std=c++11 $(CXX_RELEASE_FLAGS) -static -o proNet.exe proNet.cpp

debug: proNet.exe
	$(CC) -std=c++11 $(CXX_DEBUG_FLAGS) -o proNet.exe proNet.cpp

clean:
	rm -f proNet.exe

