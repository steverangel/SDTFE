TARGET = dtfe
LIB = libdtfe.a

QHULLLIBDIR = ../qhull/lib 
QHULLINCDIR = -I../qhull/src/libqhull

TIFFLIBDIR = ../tiff-4.0.9/install/lib
TIFFINCDIR = -I../tiff-4.0.9/libtiff

CC = cc -arch x86_64
CFLAGS = -O2 $(QHULLINCDIR) $(TIFFINCDIR)
LINKER = ld -arch x86_64 -o

ARCHIVER = ar cvq 

SRCDIR = src
OBJDIR = obj
BINDIR = bin
LIBDIR = bin/lib

SOURCES  := $(wildcard $(SRCDIR)/*.c)
INCLUDES := $(wildcard $(SRCDIR)/*.h)
OBJECTS  := $(SOURCES:$(SRCDIR)/%.c=$(OBJDIR)/%.o)

rm = rm -f

all: $(BINDIR)/$(TARGET) $(LIBDIR)/$(LIB)

$(BINDIR)/$(TARGET): $(OBJECTS) ../tiff-4.0.9/libtiff/libtiff.la
	@$(LINKER) $@ $(OBJECTS) -L$(QHULLLIBDIR) -lqhullstatic -L$(TIFFLIBDIR) -ltiff -lm
	@echo "Linking complete!"

$(LIBDIR)/$(LIB): $(OBJECTS)
	@$(ARCHIVER) $@ $(OBJECTS)
	ar d $(LIBDIR)/$(LIB) main.o
	@echo "Static library created!"

$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.c
	@$(CC) $(CFLAGS) -c $< -o $@
	@echo "Compiled "$<" successfully!"

.PHONEY: clean
clean:
	@$(rm) $(OBJECTS)
	@echo "Cleanup complete!"

.PHONEY: remove
remove: clean
	@$(rm) $(BINDIR)/$(TARGET) $(LIBDIR)/$(LIB)
	@echo "Executable removed!"

# end of Makefile
