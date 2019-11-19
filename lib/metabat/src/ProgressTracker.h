#ifndef PROGRESS_TRACKER_H_
#define PROGRESS_TRACKER_H_

#include <sys/time.h>

#include "OpenMP.h"

class ProgressTracker {
public:
	ProgressTracker(long size, long steps = 100) {
		reset(size, steps);
	}
	~ProgressTracker() {
	}

	inline void track(int count = 1) {
#pragma omp atomic
		_progress += count;
	}
	inline void reset(long size, long steps = 100) {
		_progress = 0;
		_size = size <= 0 ? 1 : size;
		_steps = steps <= 0 ? 1 : steps;
		_nextStep = getIncrement();
		gettimeofday(&_start, NULL);
	}
	inline void setProgress(int progress) {
		_progress = progress;
	}
	inline bool isStepMarker() const {
		return _progress >= _nextStep;
	}
	inline long getIncrement() const {
		return (_size + _steps - 1) / _steps;
	}
	inline float getElapsed() const {
		timeval now;
		gettimeofday(&now, NULL);
		float elapsed = (now.tv_sec - _start.tv_sec) + (now.tv_usec - _start.tv_usec) / 1000000.0; //seconds
		return elapsed;
	}

	const char *getProgress() {
		long _point = _progress;
		float fraction = 1.0 * _point / _size;
		float seconds = getElapsed();
		long secondsLeft = (long) (seconds / fraction - seconds + 0.5);
		sprintf(msg, "%0.1f%% (%ld of %ld), ETA %ld:%02ld:%02ld    ", 100.0 * fraction, _point, _size, secondsLeft / 3600, (secondsLeft % 3600) / 60,
				secondsLeft % 60);
		long increment = getIncrement();
		_nextStep = (_point / increment + 1) * increment;
		return msg;
	}

private:
	long _progress, _size, _steps, _nextStep;
	timeval _start;
	char msg[80]; // 80 chars should fit easily: "100.0% (verylargenumber of verylargenumber) xxxxx.yy min"
};

#endif
