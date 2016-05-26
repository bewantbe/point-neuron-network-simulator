#ifndef HEADER_LEGACY_LIB
#define HEADER_LEGACY_LIB

#include <string>

int ReadOneLongCmdPara(int argc, char *argv[], int pp, double *vec, int size);
int CheckDirAndCreate(const std::string &filepath);

#endif