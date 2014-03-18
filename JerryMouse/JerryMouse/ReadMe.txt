simple script:

	POINT PI, PE;
	GetCursorPos(&PI);
	recordMouseMovement(2.0);
	GetCursorPos(&PE);
	MMtoOrigin(DinoNugget);
	PlayHumanMouse(MMCoord(PI.x, PI.y), MMCoord(PE.x, PE.y), 2.0, DinoNugget);

This records 2 seconds of mouse movement, then plays it back at the same speed.