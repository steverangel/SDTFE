/*
 * A simple set of functions for creating and using bitmaps
 * */

#include <stdlib.h>
#include "bitmap.h"

enum { BITS_PER_WORD = sizeof(word_t) * sizeof(char) };

#define WORD_OFFSET(b) ((b) / BITS_PER_WORD)
#define BIT_OFFSET(b)  ((b) % BITS_PER_WORD)

void set_bit(word_t *words, int n) { 
    words[WORD_OFFSET(n)] |= (1 << BIT_OFFSET(n));
}

void clear_bit(word_t *words, int n) {
    words[WORD_OFFSET(n)] &= ~(1 << BIT_OFFSET(n)); 
}

int get_bit(word_t *words, int n) {
    word_t bit = words[WORD_OFFSET(n)] & (1 << BIT_OFFSET(n));
    return bit != 0; 
}

word_t* create_bitmap(int n) {
    return (word_t*)calloc(WORD_OFFSET(n) + 1, sizeof(word_t));
}

void delete_bitmap(word_t *words) {
    free(words);
}

