.PHONY: clean

all: test

main: 	compress.o rlcsa/rlcsa.a
	g++ -g $^ -lm -o compress

rlcsa/rlcsa.a:
	$(MAKE) -C rlcsa
	
test: main 
	./run
clean:
	-rm compress
