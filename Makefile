CC=gcc
CFLAGS=-W -Wall -ansi -pedantic -std=c11 -g
INC=-I include/
SRC=src/
TST=test/
EXE=exe/
LIBS= -lm
OBJ= $(SRC)raisin.o $(SRC)m.o $(SRC)fm.o $(SRC)dlist.o $(SRC)uf.o $(SRC)rbh.o $(SRC)crbh.o $(SRC)ipart.o $(SRC)io.o $(SRC)a.o $(SRC)matrix.o $(SRC)vector.o $(SRC)pqueue.o $(SRC)list.o
OBJTST= $(TST)test_rbh.o $(SRC)rbh.o $(SRC)crbh.o $(SRC)io.o $(SRC)a.o $(SRC)matrix.o $(SRC)vector.o $(SRC)pqueue.o $(SRC)list.o
OBJTSTC= $(TST)test_crbh.o $(SRC)rbh.o $(SRC)crbh.o $(SRC)io.o $(SRC)a.o $(SRC)matrix.o $(SRC)vector.o $(SRC)pqueue.o $(SRC)list.o
OBJTSTP= $(TST)test_ipart.o $(SRC)uf.o $(SRC)rbh.o $(SRC)crbh.o $(SRC)ipart.o $(SRC)io.o $(SRC)a.o $(SRC)matrix.o $(SRC)vector.o $(SRC)pqueue.o $(SRC)list.o
OBJTSTFM= $(TST)test_fm.o $(SRC)fm.o $(SRC)dlist.o $(SRC)uf.o $(SRC)rbh.o $(SRC)crbh.o $(SRC)ipart.o $(SRC)io.o $(SRC)a.o $(SRC)matrix.o $(SRC)vector.o $(SRC)pqueue.o $(SRC)list.o
OBJTSTM= $(TST)test_m.o $(SRC)m.o $(SRC)fm.o $(SRC)dlist.o $(SRC)uf.o $(SRC)rbh.o $(SRC)crbh.o $(SRC)ipart.o $(SRC)io.o $(SRC)a.o $(SRC)matrix.o $(SRC)vector.o $(SRC)pqueue.o $(SRC)list.o
EXEC=raisin
TEST=test_rbh

all: $(EXEC)

test_m: $(OBJTSTM)
	$(CC) $(INC) -o $(EXE)$@ $^ $(CFLAGS) $(LIBS)

test_fm: $(OBJTSTFM)
	$(CC) $(INC) -o $(EXE)$@ $^ $(CFLAGS) $(LIBS)

test_ipart: $(OBJTSTP)
	$(CC) $(INC) -o $(EXE)$@ $^ $(CFLAGS) $(LIBS)

test_crbh: $(OBJTSTC)
	$(CC) $(INC) -o $(EXE)$@ $^ $(CFLAGS)

test_rbh: $(OBJTST)
	$(CC) $(INC) -o $(EXE)$@ $^ $(CFLAGS)

raisin: $(OBJ)
	$(CC) $(INC) -o $(EXE)$@ $^ $(CFLAGS) $(LIBS)

$(SRC)%.o : $(SRC)%.c
	$(CC) $(INC) -o $@ -c $< $(CFLAGS)

clean:
	rm -rf $(SRC)*.o $(TST)*.o $(EXE)/*
