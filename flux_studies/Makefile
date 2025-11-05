ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS   := $(shell root-config --libs)
CXXFLAGS := -O2 -Wall -fPIC $(ROOTCFLAGS) -I. -I/exp/minerva/app/users/qvuong/MAT_AL9/opt/include
LDFLAGS    := $(ROOTLIBS)

DICT_NAME := FluxDict
DICT_HDR := LinkDef.h
DICT_SRC := ${DICT_NAME}.cxx
DICT_OBJ := ${DICT_NAME}.o

TARGET := libFluxLoop.so
OBJS := FluxCalculatorLoop.o $(DICT_OBJ)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) -shared -o $@ $^ $(LDFLAGS)

%.o: %.cxx
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(DICT_SRC): FluxCalculatorLoop.h $(DICT_HDR)
	rootcling -f $(DICT_SRC) -c $^

clean:
	rm -f *.o *.so *.pcm $(DICT_SRC)
