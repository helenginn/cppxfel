//
//  Masker.hpp
//  cppxfel_local
//
//  Created by Helen Ginn on 03/01/2018.
//  Copyright © 2018 Division of Structural Biology Oxford. All rights reserved.
//

#ifndef Masker_hpp
#define Masker_hpp

#include <stdio.h>
#include <vector>

typedef struct
{
	int x0;
	int x1;
	int y0;
	int y1;
} Rectangle;

class Masker
{
public:
	Masker();

	static Masker *getMasker()
	{
		return &masker;
	}
	
	bool isMasked(int x, int y);
	void addRectangle(int x0, int y0, int x1, int y1);
private:


	std::vector<Rectangle> _rectangles;
	static Masker masker;
};

#endif /* Masker_hpp */
