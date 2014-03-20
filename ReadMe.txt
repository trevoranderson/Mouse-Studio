simple script:

void test(std::string filename)
{
	auto testObj = new MouseMovement();
	testObj->Record(2.0);
	testObj->PlayCurrentMovement(2.0);
	testObj->OutputStorage(filename);
}

This records 2 seconds of mouse movement, then plays it back at the same speed, 
outputting it to the filename specified. If you wanted to play it back in the future,
call the ReadToStorage() function with the same filename.
