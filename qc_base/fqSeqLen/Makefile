all : fqSeqLen

fqSeqLen : fqSeqLen.o
	@gcc -o fqSeqLen fqSeqLen.o
fqSeqLen.o : fqSeqLen.c
	@gcc -c fqSeqLen.c
clean :
	@rm fqSeqLen fqSeqLen.o
