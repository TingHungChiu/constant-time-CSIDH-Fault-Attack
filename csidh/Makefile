CC=gcc

# List of object files
OBJS = constants.o fp.o randombytes.o uint.o csidh.o mont.o speed.o

# Executable target
test: $(OBJS)
	$(CC) $(OBJS) -o test

# Compile rule for .c -> .o  
%.o: %.c
	$(CC) -c $< -o $@

# Clean rule
clean: 
	rm -f $(OBJS) test