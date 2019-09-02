ROOT_DIR=$(shell pwd)
ODIR  = $(ROOT_DIR)/obj
SDIR  = $(ROOT_DIR)/src/main
TDIR  = $(ROOT_DIR)/test

CXX   = g++
CFLAG = -std=c++11
CLIB  = -lvec
 
DEPS  = $(shell ls $(SDIR)/*.h)
SRC   = $(shell ls $(SDIR)/*.cpp)
OBJ   = $(patsubst $(SDIR)/%.cpp,$(ODIR)/%.o,$(SRC))

.PHONY: all
all :
	make bvvv.x
	make mc_aux.x

bvvv.x : $(OBJ)
	$(CXX) -o $@ $^ $(CFLAG) $(CLIB)
	cp $@ $(TDIR)

mc_aux.x : $(ROOT_DIR)/src/aux/mc_convertor.cpp
	$(CXX) -o $@ $^ $(CFLAG) $(CLIB)
	cp $@ $(TDIR)


$(ODIR)/%.o : $(SDIR)/%.cpp $(DEPS) | $(ODIR)/.
	$(CXX) -c -o $@ $< $(CFLAG) $(CLIB)

%/. : 
	mkdir -p $(patsubst %/.,%,$@)

.PRECIOUS: %/.
.PHONY: clean
clean:
	rm -rf *.x test/*.x $(ODIR)
