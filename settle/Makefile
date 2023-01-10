# settle:  settle.o
# 	g++ -fPIC -c *.cc *.c
# 	g++ -shared -fPIC -Wl,-install_name,libsettle.so  -o libsettle.so *.o
# 	#objdump -T libsettle.so | grep main

mac: settle.o
#	g++ -Ofast -fPIC -c *.cc *.c
#	g++ -Ofast -shared -fPIC -Wl,-install_name,libsettle.so  -o libsettle.so *.o
	g++ -fPIC -c *.cc *.c
	g++ -shared -fPIC -Wl,-install_name,libsettle.so  -o libsettle.so *.o
	#objdump -T libsettle.so | grep main

linux: settle.o 
	gcc -g -Wall -O0 -fPIC -c *.c
	g++ -g -Wall -O0 -fPIC -c *.cc
	g++ -O0 -shared -fPIC -Wl,-soname,libsettle.so  -o libsettle.so *.o
	#objdump -T libsettle.so | grep main

clean:
	rm -fv *.o

cleaner: clean
	rm -fv *.so *.pyc
