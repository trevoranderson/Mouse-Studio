#pragma once
#include <vector>
#include "Point.h"
#include "Mat2.h"
#define DISTLENGTH 1000
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
			time_percent = (TimeDifference(now_time, start_time)) / time_to_move;
			SetCursorPos(start.x + time_percent*(parametricvector.x), start.y + time_percent*(parametricvector.y));
		}
	}
	void Click()
	{
		POINT * curserpos = new POINT();
		GetCursorPos(curserpos);
		INPUT input;
		input.type = INPUT_MOUSE;
		input.mi.dx = curserpos->x;
		input.mi.dy = curserpos->y;
		input.mi.dwFlags = (MOUSEEVENTF_ABSOLUTE |/*MOUSEEVENTF_MOVE|*/MOUSEEVENTF_LEFTDOWN | MOUSEEVENTF_LEFTUP);
		input.mi.mouseData = 0;
		input.mi.dwExtraInfo = NULL;
		input.mi.time = 0;
		SendInput(1, &input, sizeof(INPUT));
	}
	void PressKey(int key = 0x75, double waitTime = 0, double waitRadius =0)
	{
		// Calculate wait Time inbetween key down and key up
		double sleepDuration = BellRand(waitTime, waitRadius);

		INPUT ip;
		ip.type = INPUT_KEYBOARD;
		ip.ki.wScan = 0; // hardware scan code for key
		ip.ki.time = 0;
		ip.ki.dwExtraInfo = 0;

		// Press the key
		ip.ki.wVk = key; // virtual-key code for the f6 key -> replace later
		ip.ki.dwFlags = 0; // 0 for key press
		SendInput(1, &ip, sizeof(INPUT));

		// Sleep for specified time
		Sleep(sleepDuration);

		// Release the key
		ip.ki.dwFlags = KEYEVENTF_KEYUP; // KEYEVENTF_KEYUP for key release
		SendInput(1, &ip, sizeof(INPUT));
	}
	void PlayCurrentMovement(double time_to_move)
	{
		SYSTEMTIME start_time;
		GetSystemTime(&start_time);
		SYSTEMTIME now_time;
		double elapsedTime = 0;
		int bucket;
		double numPoints = Storage.size();
		while (elapsedTime<1)
		{
			GetSystemTime(&now_time);
			elapsedTime = (TimeDifference(now_time, start_time)) / time_to_move;
			//select the proper bucket for smooth movement
			bucket = (double)numPoints*elapsedTime;
			if (bucket >= numPoints)
			{
				bucket = numPoints - 1;
			}
			SetCursorPos(Storage[bucket].x, Storage[bucket].y);
		}
	}
	void PlayCurrentBetweenPoints(Point begin, Point end, double time_to_move)
	{
			//First scale -> rotate -> translate such that you cover the wanted path
			//Scale:
			Point dirVec(end.x - begin.x, end.y - begin.y);
			double displacement = sqrt((dirVec.x) * (dirVec.x) + (dirVec.y) * (dirVec.y));
			double scaleFactor = displacement / (double)DISTLENGTH;
			int numPoints = Storage.size();
			for (int cnt = 0; cnt<numPoints; cnt++)
			{
				Storage[cnt] = vecmatmult(matscale(scaleFactor, scaleFactor), Storage[cnt]);
			}
			//Rotate:
			double theta = atan2(dirVec.y, dirVec.x);
			for (int cnt = 0; cnt<numPoints; cnt++)
			{
				Storage[cnt] = vecmatmult(matrot(theta), Storage[cnt]);
			}
			// Translate:
			for (int cnt = 0; cnt<numPoints; cnt++)
			{
				Storage[cnt].x += begin.x;
				Storage[cnt].y += begin.y;
			}
			//PlayMouseMovement(time_to_move);
			/////////////////////////////////////////////////////////////////////////
			SYSTEMTIME start_time;
			GetSystemTime(&start_time);
			SYSTEMTIME now_time;
			double elapsedTime = 0;
			int bucket;
			double largenum;
			while (elapsedTime<1)
			{
				GetSystemTime(&now_time);
				elapsedTime = (TimeDifference(now_time, start_time)) / time_to_move;
				//select the proper bucket for smooth movement

				largenum = (double)numPoints*elapsedTime;
				if (largenum > INT_MAX)
				{
					bucket = numPoints - 1;
				}
				else
				{
					bucket = largenum;
				}


				if (bucket >= numPoints)
				{
					bucket = numPoints - 1;
				}
				SetCursorPos(Storage[bucket].x, Storage[bucket].y);
			}
			ToStorageForm();
		}
	void Record(double time_to_move)
	{
		// This algorithm works by filling a large vector of Points
		// according to the time elapsed, then filling the the gaps by 
		// linear interpolation



		//Allocate enough Storage that the vector to have the correct resolution
		int sizeToUse = pointspersecond * time_to_move;
		Storage.reserve(sizeToUse);
		while (Storage.size() < sizeToUse)
		{
			Storage.push_back(Point());
		}
		int numPoints = Storage.size();

		POINT Q;
		POINT * curserposinit = &Q;
		SYSTEMTIME start_time;
		GetSystemTime(&start_time);
		SYSTEMTIME now_time;
		double elapsedTime = 0;
		while (elapsedTime<1.0)
		{
			GetSystemTime(&now_time);
			elapsedTime = (TimeDifference(now_time, start_time)) / time_to_move;
			GetCursorPos(curserposinit);
			int bucket = (int)numPoints*elapsedTime;
			if (bucket < numPoints)
			{
				Storage[bucket] = Point(curserposinit->x, curserposinit->y);
			}
			else
			{
				Storage[numPoints -1 ] = Point(curserposinit->x, curserposinit->y);
			}
			
		}
		// Interpolate
		cleanMM();
	}
	void ToStorageForm()
	{
		//Translates to origin-> rotates theta to zero -> scales to 500 total pixels in X-axis
		// first translate all points to origin
		Point first(Storage[0].x, Storage[0].y);
		int numPoints = Storage.size();
		for (int cnt = 0; cnt<numPoints; cnt++)
		{
			Storage[cnt].x -= first.x;
			Storage[cnt].y -= first.y;
		}
		//find the angle between the ending points and set it equal to zero.
		double theta = atan2((Storage[numPoints - 1].y), (Storage[numPoints - 1].x));
		for (int cnt = 0; cnt<numPoints; cnt++)
		{
			Storage[cnt] = vecmatmult(matrot(-theta), Storage[cnt]);
		}
		//finally, scale the length to a predetermined amount
		double scaleFactor = (double)DISTLENGTH / (Storage[numPoints - 1].x);
		for (int cnt = 0; cnt<numPoints; cnt++)
		{
			Storage[cnt] = vecmatmult(matscale(scaleFactor, scaleFactor), Storage[cnt]);
		}
	}
private:
	void cleanMM()
	{
		//Here we want to take our array of mouse positions and interpolate between measurements
		Point Previous(0, 0);
		Point Current(0, 0);
		int firstUnfilled = 0;
		int lastUnfilled = 0;
		int totUnfilled = 0;
		int numPoints = Storage.size();
		if (Storage[0].x <0 && Storage[0].y <0)
		{
			for (int k = 1; k<numPoints - 1; k++)
			{
				if (Storage[k].x >= 0 && Storage[k].y >= 0)
				{
					Storage[0] = Storage[k];
					cleanMM();
					return;
				}
			}
		}
		for (int cnt = 0; cnt<numPoints- 1; cnt++)
		{
			if (Storage[cnt].x<0 && Storage[cnt].y<0)
			{
				firstUnfilled = cnt;
				lastUnfilled = cnt;
				Previous.x = Storage[cnt - 1].x;
				Previous.y = Storage[cnt - 1].y;
				for (cnt = cnt + 1; Storage[cnt].x<0 && Storage[cnt].y<0; cnt++)
				{
					lastUnfilled++;
				}
				if (cnt >= numPoints - 1)
				{
					break;
				}
				//Now we have Previous, first unfilled, last unfilled
				Current.x = Storage[cnt].x;
				Current.y = Storage[cnt].y;

				//now loop through unfilled buckets, interpolating the values
				double divider = lastUnfilled - firstUnfilled + 2; // if we have 3,5, we need to divide the gap into 4 spaces
				Point dirVec = Point(Current.x - Previous.x, Current.y - Previous.y);
				for (int p = 0; p<divider - 1; p++)
				{
					Storage[p + firstUnfilled].x = dirVec.x*((p + 1) / divider) + Previous.x;
					Storage[p + firstUnfilled].y = dirVec.y*((p + 1) / divider) + Previous.y;

				}
				totUnfilled += lastUnfilled - firstUnfilled + 1;

			}
		}
		//if the last value is garbage, we need to put a value there
		// and loop back through the array in order to fill things.
		if (Storage[numPoints - 1].x<0 && Storage[numPoints - 1].y<0)
		{
			Storage[numPoints - 1].x = Previous.x;
			Storage[numPoints - 1].y = Previous.y;
			cleanMM();
		}
	}
	Point generateDest(Point center, double radius)
	{
		// The goal of this function is to give our move target a valid location
		// that follows a normal distribution. radius=max radius
		double randRadius = BellRand(radius*(double)2.0, 0);
		if (randRadius == 0)
		{
			return center;
		}
		double theta = rand();
		double retX = sin(theta) * randRadius + center.x;
		double retY = cos(theta) * randRadius + center.y;
		return Point(retX, retY);
	}
	// Turn Storage into the valid storage form (goes from (0,0) to (500,0))

	// Resolution
	int pointspersecond = 1000;
	// Holds all the structs.
	std::vector<Point> Storage;
};