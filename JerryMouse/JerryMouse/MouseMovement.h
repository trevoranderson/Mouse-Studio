#pragma once
#include <vector>
#include <fstream>
#include <stdio.h>
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
	void Click(double waitTime = 0, double waitRadius = 0)
	{
		if (waitTime == 0 && waitRadius == 0)
		{
			POINT p;
			POINT * curserpos = &p;
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
		else
		{
			//Calculate sleep time
			double sleepDuration = BellRand(waitTime, waitRadius);

			POINT p;
			POINT * curserpos = &p;
			GetCursorPos(curserpos);
			INPUT input;
			input.type = INPUT_MOUSE;
			input.mi.dx = curserpos->x;
			input.mi.dy = curserpos->y;
			input.mi.dwFlags = (MOUSEEVENTF_ABSOLUTE |/*MOUSEEVENTF_MOVE|*/MOUSEEVENTF_LEFTDOWN);
			input.mi.mouseData = 0;
			input.mi.dwExtraInfo = NULL;
			input.mi.time = 0;
			SendInput(1, &input, sizeof(INPUT));

			//wait
			Sleep(sleepDuration);
			input.mi.dwFlags = (MOUSEEVENTF_ABSOLUTE | MOUSEEVENTF_LEFTUP);
			SendInput(1, &input, sizeof(INPUT));
		}
	}
	void PressKey(int key = 0x75, double waitTime = 0, double waitRadius = 0)
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
		int numPoints = Storage.size();
		int bucket = 0;
		while (bucket < numPoints)
		{
			// Apply the current bucket's action
			SetCursorPos(Storage[bucket].x, Storage[bucket].y);
			// Prepare for next
			GetSystemTime(&now_time);
			elapsedTime = (TimeDifference(now_time, start_time)) / time_to_move;
			//select the proper bucket for smooth movement
			bucket = (double)numPoints*elapsedTime;
		}
	}
	void PlayCurrentBetweenPoints(Point begin, Point end, double time_to_move)
	{
		// Since we are doing transformations, build a tmp storage vector from the old
		std::vector<Point> toPlay = Storage;
		//scale -> rotate -> translate such that you cover the wanted path
		//Scale:
		Point dirVec(end.x - begin.x, end.y - begin.y);
		double displacement = sqrt((dirVec.x) * (dirVec.x) + (dirVec.y) * (dirVec.y));
		double scaleFactor = 1.0;
		int numPoints = toPlay.size();
		for (int cnt = 0; cnt < numPoints; cnt++)
		{
			toPlay[cnt] = vecmatmult(matscale(scaleFactor, scaleFactor), toPlay[cnt]);
		}
		//Rotate:
		double theta = atan2(dirVec.y, dirVec.x);
		for (int cnt = 0; cnt < numPoints; cnt++)
		{
			toPlay[cnt] = vecmatmult(matrot(theta), toPlay[cnt]);
		}
		// Translate:
		for (int cnt = 0; cnt < numPoints; cnt++)
		{
			toPlay[cnt].x += begin.x;
			toPlay[cnt].y += begin.y;
		}
		//PlayMouseMovement(time_to_move);
		/////////////////////////////////////////////////////////////////////////
		SYSTEMTIME start_time;
		GetSystemTime(&start_time);
		SYSTEMTIME now_time;
		double elapsedTime = 0;
		int bucket = 0;
		while (bucket < numPoints)
		{
			// Apply the current bucket's action
			SetCursorPos(toPlay[bucket].x, toPlay[bucket].y);
			// Prepare for next
			GetSystemTime(&now_time);
			elapsedTime = (TimeDifference(now_time, start_time)) / time_to_move;
			//select the proper bucket for smooth movement
			bucket = (double)numPoints*elapsedTime;
		}
	}
	void Record(double time_to_move, int resolutionpps = 1000)
	{
		// This algorithm works by filling a large vector of Points
		// according to the time elapsed, then filling the the gaps by 
		// linear interpolation

		ClearStorage();
		//Allocate enough Storage that the vector to have the correct resolution
		int sizeToUse = resolutionpps * time_to_move;
		Storage.resize(sizeToUse);
		int numPoints = Storage.size();

		POINT Q;
		POINT * curserposinit = &Q;
		SYSTEMTIME start_time;
		GetSystemTime(&start_time);
		SYSTEMTIME now_time;
		double elapsedTime = 0;
		// This structure guarantees that Storage[0] is nonempty
		// and we never attempt to put something out of range
		int bucket = 0;
		while (bucket < numPoints)
		{
			// Put something in the bucket
			GetCursorPos(curserposinit);
			Storage[bucket] = Point(curserposinit->x, curserposinit->y);
			//Calculate next bucket
			GetSystemTime(&now_time);
			elapsedTime = (TimeDifference(now_time, start_time)) / time_to_move;
			bucket = (int)numPoints*elapsedTime;
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
		for (int cnt = 0; cnt < numPoints; cnt++)
		{
			Storage[cnt].x -= first.x;
			Storage[cnt].y -= first.y;
		}
		//find the angle between the ending points and set it equal to zero.
		double theta = atan2((Storage[numPoints - 1].y), (Storage[numPoints - 1].x));
		for (int cnt = 0; cnt < numPoints; cnt++)
		{
			Storage[cnt] = vecmatmult(matrot(-theta), Storage[cnt]);
		}
		//finally, scale the length to a predetermined amount
		double scaleFactor = (double)DISTLENGTH / (Storage[numPoints - 1].x);
		for (int cnt = 0; cnt < numPoints; cnt++)
		{
			Storage[cnt] = vecmatmult(matscale(scaleFactor, scaleFactor), Storage[cnt]);
		}
	}
	void OutputStorage(std::string path)
	{
		// copy Storage into a large char array.
		std::ofstream file(path, std::ios::binary);
		int StorageSize = Storage.size();
		int bytesize = sizeof(Point)*(Storage.size());
		char * bytearr = new char[bytesize];
		memcpy(bytearr, &Storage[0], sizeof(Point)* StorageSize);
		file.write(bytearr, bytesize);
		file.close();
		delete bytearr;
	}
	void ReadToStorage(std::string path)
	{
		FILE * file = NULL;
		if ((file = fopen(path.c_str(), "rb")) == NULL)
		{
			std::cout << "Could not open specified file" << std::endl;
		}
		long lCurPos, lEndPos;
		lCurPos = ftell(file);
		fseek(file, 0, 2);
		lEndPos = ftell(file);
		fseek(file, lCurPos, 0);
		
		int fileSize = lEndPos;
		void * fileBuff = new char[fileSize];
		fread(fileBuff, fileSize, 1, file);

		//Cast it to Points and fill Storage
		int StorageSize = fileSize / (sizeof(Point));
		Storage.resize(StorageSize);
		Point * ArrToCpy = (Point*)fileBuff;
		for (int k = 0; k < StorageSize; k++)
		{
			Storage[k].x = ArrToCpy[k].x;
			Storage[k].y = ArrToCpy[k].y;
		}
		fclose(file);
		delete fileBuff;
	}
private:
	void ClearStorage()
	{
		for (auto &i : Storage)
		{
			i.x = -1;
			i.y = -1;
		}
	}
	void cleanMM()
	{
		//Here we want to take our array of mouse positions and interpolate between measurements
		Point Previous(0, 0);
		Point Current(0, 0);
		int firstUnfilled = 0;
		int lastUnfilled = 0;
		int totUnfilled = 0;
		int numPoints = Storage.size();
		// The bucket at index zero is guaranteed to have something, if it was recorded using our Record
		for (int cnt = 0; cnt < numPoints - 1; cnt++)
		{
			if (Storage[cnt].x < 0 && Storage[cnt].y < 0)
			{
				firstUnfilled = cnt;
				lastUnfilled = cnt;
				Previous.x = Storage[cnt - 1].x;
				Previous.y = Storage[cnt - 1].y;
				for (cnt = cnt + 1; Storage[cnt].x < 0 && Storage[cnt].y < 0; cnt++)
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
				for (int p = 0; p < divider - 1; p++)
				{
					Storage[p + firstUnfilled].x = dirVec.x*((p + 1) / divider) + Previous.x;
					Storage[p + firstUnfilled].y = dirVec.y*((p + 1) / divider) + Previous.y;

				}
				totUnfilled += lastUnfilled - firstUnfilled + 1;

			}
		}
		//if the last value is garbage, we need to put a value there
		// and loop back through the array in order to fill things.
		if (Storage[numPoints - 1].x < 0 && Storage[numPoints - 1].y < 0)
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
	// Holds all the structs.
	std::vector<Point> Storage;
};