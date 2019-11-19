#ifndef SRC_TILE_H_
#define SRC_TILE_H_

#include <stddef.h>
#if defined(__APPLE__)

#include <sys/sysctl.h>
size_t CacheSize() {
	size_t lineSize = 0;
	size_t sizeOfLineSize = sizeof(lineSize);
	sysctlbyname("hw.cachelinesize", &lineSize, &sizeOfLineSize, 0, 0);
	return lineSize;
}

#elif defined(_WIN32)

#include <stdlib.h>
#include <windows.h>
size_t CacheSize() {
	size_t lineSize = 0;
	DWORD bufferSize = 0;
	DWORD i = 0;
	SYSTEM_LOGICAL_PROCESSOR_INFORMATION * buffer = 0;

	GetLogicalProcessorInformation(0, &bufferSize);
	buffer = (SYSTEM_LOGICAL_PROCESSOR_INFORMATION *) malloc(bufferSize);
	GetLogicalProcessorInformation(&buffer[0], &bufferSize);

	for (i = 0; i != bufferSize / sizeof(SYSTEM_LOGICAL_PROCESSOR_INFORMATION); ++i) {
		if (buffer[i].Relationship == RelationCache && buffer[i].Cache.Level == 2) {
			lineSize = buffer[i].Cache.Size;
			break;
		}
	}

	free(buffer);
	return lineSize;
}

#elif defined(__linux__)

#include <stdio.h>
size_t CacheSize() {
	FILE * p = 0;
	const char *cs = "/sys/devices/system/cpu/cpu0/cache/index2/size";
	p = fopen(cs, "r");
	unsigned int lineSize = 0;
	if (p) {
		int ret = fscanf(p, "%d", &lineSize);
		if (ret == EOF) { /* do not use it */ fprintf(stderr, "Warning: could not read %s to get cache sizes, using defaults\n", cs); }
		fclose(p);
	}
	return lineSize;
}

#else

size_t CacheSize() {
	return 0;	
}

#endif

#endif /* SRC_TILE_H_ */
