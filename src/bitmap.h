#ifndef _BITMAP_h
#define _BITMAP_h

#include <stdint.h>   // for uint32_t

typedef uint32_t word_t;    

// bitmap function prototypes
void set_bit(word_t *, int); 

void clear_bit(word_t *, int );

int get_bit(word_t *, int );

word_t* create_bitmap(int);

void delete_bitmap(word_t *);

#endif
