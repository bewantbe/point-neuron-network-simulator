#ifndef HEADER_LEGACY_LIB
#define HEADER_LEGACY_LIB

int ReadOneLongCmdPara(int argc, char *argv[], int pp, double *vec, int size);
int CheckDirAndCreate(const std::string &filepath);
unsigned int GetSeedFromTime();

#endif