CC              = gcc
SRC             = Ovarian.c 
TARGET          = project
OBJ             = $(SRC:.c=.o)
DDDEBUG         =
#------------------------------------------------------

$(TARGET): $(OBJ)
	$(CC)  -o $@ $(OBJ)

$(OBJ): $(SRC)
	$(CC) -g -c $(SRC) $(DDDEBUG)

clean:
	rm -f *.o *~ $(TARGET)
	rm -f Result.txt
	rm -f ResultAUC.txt
