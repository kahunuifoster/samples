all: randmst 0d

%: %.cpp
	g++-4.8 --std=c++11 $< -o $@