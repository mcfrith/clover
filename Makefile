CXX = g++
CXXFLAGS = -Wall -O3 -DNDEBUG

clover: *.cpp *.hpp
	$(CXX) $(CXXFLAGS) -o $@ *.cpp

clean:
	rm -f clover

distdir = clover-`date +%Y-%m-%d`

dist:
	mkdir $(distdir)
	cp *.?pp Makefile $(distdir)
	zip -qrm $(distdir) $(distdir)
