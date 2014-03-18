#pragma once
#include <vector>
#include "Point.h"
double BellRand(double width, double center)
{
	// Generates a random number around a given center and a distribution max length
	double randSum = 0;
	for (int k = 0; k < 100; k++)
	{
		randSum += rand() % 100;
	}
	randSum = randSum / (double)100.0; //now have a normally distrubuted number between 0 and 100,000
	randSum = randSum / (double)100.0; // now its between 0-1
	randSum *= width;
	randSum -= (width / 2);
	randSum += center;
	return randSum;
}
inline double TimeDifference(SYSTEMTIME A, SYSTEMTIME B)
{
	double ret = (A.wHour - B.wHour) * 3600 + (A.wMinute - B.wMinute) * 60 + A.wSecond - B.wSecond + (double)(A.wMilliseconds - B.wMilliseconds) / 1000;
	return ret;
}
class MouseMovement
{
public:
	void LinearMove(Point start, Point end, double time_to_move)
	{
		if (start.x == end.x && start.y == end.y)
		{
			return;
		}
		SYSTEMTIME start_time;
		GetSystemTime(&start_time);
		Point parametricvector(end.x - start.x, end.y - start.y);
		SYSTEMTIME now_time;
		double time_percent = 0;
		while (time_percent < 1)
		{
			GetSystemTime(&now_time);
			time_percent = (dif_frac_seconds(now_time, start_time)) / time_to_move;
			SetCursorPos(start.x + time_percent*(parametricvector.x), start.y + time_percent*(parametricvector.y));
		}
	}
private:
	Point generateDest(Point center, double radius)
	{
		// The goal of this function is to give our move target a valid location
		// that follows a normal distribution. radius=max radius
		double randRadius = genNormalDist(radius*(double)2.0, 0);
		if (randRadius == 0)
		{
			return center;
		}
		double theta = rand();
		double retX = sin(theta) * randRadius + center.x;
		double retY = cos(theta) * randRadius + center.y;
		return Point(retX, retY);
	}
	// Resolution
	int pointspersecond;
	// Holds all the structs.
	std::vector<Point> Storage;
};