settle:  settle.o
	g++ -fPIC -c *.cc *.c
	g++ -shared -fPIC -Wl,-soname,libsettle.so  -o libsettle.so *.o
	#objdump -T libsettle.so | grep main

clean:
	rm -fv *.o

cleaner: clean
	rm -fv *.so *.pyc
