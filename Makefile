############################# Makefile ##########################

LIB = lib/
INCLUDE	=	headers/
SOURCES	=	sources/

SRC1	=	SystemOfParticles.cpp
SRC2	=	Force.cpp
SRC3    =   etime.cpp
MAIN  = main.cpp

SRCS	=	$(SRC3) $(SRC2) $(SRC1) $(MAIN)

INC1	=	$(SRC1:.cpp=.hpp)
INC2	=	$(SRC2:.cpp=.hpp)
INC3    =   $(SRC3:.cpp=.hpp)

OBJS	=	$(addprefix $(LIB),$(SRCS:.cpp=.o))

CXX	=   g++

CPPFLAGS    =	-std=c++11 -I $(INCLUDE)

EXEC = rmoldyn

build: $(EXEC)

$(LIB)$(SRC1:.cpp=.o): $(SOURCES)$(SRC1) $(INCLUDE)$(INC1) $(INCLUDE)$(INC3)
	$(CXX) $(CPPFLAGS)	-c $< -o $@
$(LIB)$(SRC2:.cpp=.o): $(SOURCES)$(SRC2) $(INCLUDE)$(INC2)
	$(CXX) $(CPPFLAGS)	-c $< -o $@
$(LIB)$(SRC3:.cpp=.o): $(SOURCES)$(SRC3) $(INCLUDE)$(INC3)
	$(CXX) $(CPPFLAGS)	-c $< -o $@
$(LIB)$(MAIN:.cpp=.o): $(MAIN)
	$(CXX) $(CPPFLAGS)	-c $< -o $@
	
$(EXEC): $(OBJS)
	@echo $(OBJS)
	$(CXX) -o $@ $^

	@echo 	"Program built."
	
.PHONY: clean clean_all rebuild

clean:
	@echo "Object files removed!"
	@rm -f $(LIB)*.o
clean_all:
	@echo "Object files and executable removed!"
	@rm -f $(LIB)*.o $(EXEC) *~ $(SOURCES)/*~ $(INCLUDE)/*~
rebuild:
	@make clean_all
	@make
	
run:
	@./$(EXEC)
