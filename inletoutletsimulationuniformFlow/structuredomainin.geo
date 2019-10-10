	// base corner points
	Point(1) = {0,0,0};
        Point(2) = {0,4,0};
        Point(3) = {0,8,0};
	Point(4) = {0,12,0};
	Point(5) = {4,12,0};
        Point(6) = {8,12,0};
        Point(7) = {12,12,0};

	Point(8) = {12,8,0};
	Point(9) = {12,4,0};
	Point(10) = {12,0,0};
	Point(11) = {8,0,0};
	Point(12) = {4,0,0};
	
	// base lines
	Line(1) = {1, 2};
	Line(2) = {2, 3};
	Line(3) = {3, 4};
	Line(4) = {4, 5};
	Line(5) = {5, 6};
	Line(6) = {6, 7};
	Line(7) = {7, 8};
	Line(8) = {8, 9};
        Line(9) = {9, 10};
	Line(10) ={10,11};
	Line(11) = {11,12};
        Line(12) = {12,1};


	//To add surface
	Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8,9,10,11,12};
	
	Plane Surface(1) = {1};
	
	// extra work for creating structured mesh

	  // Transfinite Line{1, 9} = 33;
      //   Transfinite Line{2, 8} = 33;
           //Transfinite Line{3, 7} = 33;
	   //Transfinite Line{4, 6} = 33;
	  // Transfinite Line{5, 12} = 33;
	  // Transfinite Line{6, 11} = 33;  
	  // Transfinite Line{7, 10} = 33;
        
         

           Transfinite Curve{1, 2, 3, 4, 5, 6, 7, 8, 9, 10,11, 12} = 33;

	   Transfinite Surface{1}= {1,10,7,4} ;
	     Recombine Surface{1};


	// extrude base surface in z direction for create 


	// dummy 3d geometry for OpenFoam
	Extrude {0, 0, 0.5} {
	   Surface{1}; Layers{1}; Recombine;
	}

        // Physical Surface
	
	Physical Surface("inlet") = {69};
	Physical Surface("outlet") = {57};
	Physical Surface("frontAndBack") = {74, 1};
	Physical Surface("leftWall") = {73, 29, 33, 37, 41, 45, 49, 53};
	Physical Surface("rightWall") = {61, 65};
        Physical Volume("body") = {1};
