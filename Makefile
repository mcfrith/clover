CXX = g++
CXXFLAGS = -Wall -O3 -DNDEBUG

clover: *.cpp *.hpp
	$(CXX) $(CXXFLAGS) -o $@ *.cpp

clean:
	rm -f clover

log:
	git2cl > ChangeLog.txt

distdir = clover-`git log --date=short --pretty=format:%cd -1`

dist: log
	mkdir $(distdir)
	cp *.?pp Makefile *.txt $(distdir)
	tar -cf $(distdir).tar $(distdir)
	rm -r $(distdir)
	gzip $(distdir).tar
