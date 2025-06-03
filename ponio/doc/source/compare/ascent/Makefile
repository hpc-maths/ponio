GIT_repository = https://github.com/AnyarInc/Ascent.git
ASCENT_include = Ascent

all: lorenz.txt transport.txt arenstorf.txt


%.txt: %.out
	./$<

%.out: %.cpp $(ASCENT_include)
	${CXX} -std=c++20 $< -I $(ASCENT_include)/include -o $@

$(ASCENT_include):
	git clone $(GIT_repository)

clean:
	rm -rf *.txt *.out $(ASCENT_include)
