//Taken from stackoverflow for use in parallel computing applications
//https://stackoverflow.com/questions/2513505/how-to-get-available-memory-c-g
#ifndef MEMHELPER_H
#define MEMHELPER_H

#include <unistd.h>

unsigned long long getTotalSystemMemory()
{
    long pages = sysconf(_SC_PHYS_PAGES);
    long page_size = sysconf(_SC_PAGE_SIZE);
    return pages * page_size;
}

#endif
