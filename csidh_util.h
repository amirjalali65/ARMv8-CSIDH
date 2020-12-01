#ifndef CSIDH_UTIL_H
#define CSIDH_UTIL_H
#include "csidh_api.h"
#include "arith.h"

#define __STDC_FORMAT_MACROS
#include <inttypes.h>

#define VERSION 0.1

void pprint_pk(public_key *x);

void pprint_sk(private_key *y);

void pprint_ss(uint64_t *z);

void save_file(char *file, void *buf, size_t len);

void error_exit(char *str);

int read_file(const char *file, uint8_t *buf, size_t len);

int read_stdin(uint8_t *buf, int len);
#endif 
