all: shake128

shake128: shake128.cpp
	g++ -O3 -o shake128 shake128.cpp