ROOT_DIR=$(shell pwd)
ODIR  = $(ROOT_DIR)/obj
SDIR  = $(ROOT_DIR)/src

CXX   = g++
CFLAG = -std=c++11
CLIB  = -lvec
 
DEPS  = $(shell ls $(SDIR)/*.h)
SRC   = $(shell ls $(SDIR)/*.cpp)
OBJ   = $(patsubst $(SDIR)/%.cpp,$(ODIR)/%.o,$(SRC))

bvvv.x : $(OBJ)
	$(CXX) -o $@ $^ $(CFLAG) $(CLIB)

$(ODIR)/%.o : $(SDIR)/%.cpp $(DEPS) | $(ODIR)/.
	$(CXX) -c -o $@ $< $(CFLAG) $(CLIB)

%/. : 
	mkdir -p $(patsubst %/.,%,$@)
	
.PRECIOUS: %/.
.PHONY: clean
clean:
	rm -rf *.x $(ODIR)
