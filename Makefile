all:
	g++ src_rm/main.cc -Wall  -std=c++0x -O3 src_rm/allocator.cc src_rm/utils.cc src_rm/TimGraph.cc src_rm/anyoption.cc src_rm/sfmt/SFMT.c  -o main_rm
