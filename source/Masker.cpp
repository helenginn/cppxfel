//
//  Masker.cpp
//  cppxfel_local
//
//  Created by Helen Ginn on 03/01/2018.
//  Copyright © 2018 Division of Structural Biology Oxford. All rights reserved.
//

#include "Masker.h"

Masker Masker::masker;

Masker::Masker()
{
	
}

bool Masker::isMasked(int x, int y)
{
	for (int i = 0; i < _rectangles.size(); i++)
	{
		Rectangle *rect = &_rectangles[i];

		if ((x >= rect->x0 && x <= rect->x1) &&
			(y >= rect->y0 && y <= rect->y1))
		{
			return true;
		}
	}

	return false;
}

void Masker::addRectangle(int x0, int y0, int x1, int y1)
{
	Rectangle rect;
	rect.x0 = x0;
	rect.x1 = x1;
	rect.y0 = y0;
	rect.y1 = y1;
	_rectangles.push_back(rect);
}
