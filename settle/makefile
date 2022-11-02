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
	g++ -Ofast -fPIC -c *.cc *.c
	g++ -Ofast -shared -fPIC -Wl,-soname,libsettle.so  -o libsettle.so *.o
	#objdump -T libsettle.so | grep main

clean:
	rm -fv *.o

cleaner: clean
	rm -fv *.so *.pyc
